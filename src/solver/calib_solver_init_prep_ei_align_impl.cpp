// iKalibr: Unified Targetless Spatiotemporal Calibration Framework
// Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China
// https://github.com/Unsigned-Long/iKalibr.git
//
// Author: Shuolong Chen (shlchen@whu.edu.cn)
// GitHub: https://github.com/Unsigned-Long
//  ORCID: 0000-0002-5283-9057
//
// Purpose: See .h/.hpp file.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * The names of its contributors can not be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "solver/calib_solver.h"
#include "util/tqdm.h"
#include "spdlog/spdlog.h"
#include "calib/calib_data_manager.h"
#include "calib/calib_param_manager.h"
#include "viewer/viewer.h"
#include "util/status.hpp"
#include "core/tracked_event_feature.h"
#include "core/event_trace_sac.h"
#include "calib/estimator.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepEventInertialAlign() const {
    if (!Configor::IsEventIntegrated()) {
        return;
    }
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    /**
     * we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are
     * poor
     */
    const double st = std::max(so3Spline.MinTime(), scaleSpline.MinTime()) +  // the max as start
                      Configor::Prior::TimeOffsetPadding;
    const double et = std::min(so3Spline.MaxTime(), scaleSpline.MaxTime()) -  // the min as end
                      Configor::Prior::TimeOffsetPadding;

    /**
     * we first perform event-based feature tracking.
     * We first dedistort the events, then export them to files and use third-party software for
     * feature tracking.
     */
    constexpr double TRACKING_LEN_PERCENT_THD = 0.3;
    constexpr double TRACKING_FIT_SAC_THD = 5.0;
    constexpr double TRACKING_AGE_PERCENT_THD = 0.3;
    constexpr double TRACKING_FREQ_PERCENT_THD = 0.3;
    // topic, batch index, tracked features
    std::map<std::string, std::map<int, EventFeatTrackingBatch>> eventFeatTrackingRes;
    for (const auto &[topic, eventMes] : _dataMagr->GetEventMeasurements()) {
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        // create a workspace for event-based feature tracking
        const std::string hasteWorkspace =
            Configor::DataStream::OutputPath + "/events/" + topic + "/haste_ws";
        if (!std::filesystem::exists(hasteWorkspace)) {
            if (!std::filesystem::create_directories(hasteWorkspace)) {
                throw Status(Status::CRITICAL,
                             "can not create output directory '{}' for event camera '{}'!!!",
                             hasteWorkspace, topic);
            }
        }

        if (auto eventsInfo = HASTEDataIO::TryLoadEventsInfo(hasteWorkspace);
            eventsInfo != std::nullopt) {
            spdlog::info("try to load feature tracking results from haste for camera '{}'...",
                         topic);
            // the output tracking results from HASTE should be distortion-free?
            auto tracking =
                HASTEDataIO::TryLoadHASTEResults(*eventsInfo, _dataMagr->GetRawStartTimestamp());
            if (tracking != std::nullopt) {
                auto bar = std::make_shared<tqdm>();
                int barIndex = 0;
                spdlog::info("rep-process loaded event tracking by haste for '{}'...", topic);
                for (auto iter = tracking->begin(); iter != tracking->end(); ++iter) {
                    bar->progress(barIndex++, static_cast<int>(tracking->size()));
                    auto &[index, batch] = *iter;
                    // aligned time (start and end)
                    const auto &batchInfo = eventsInfo->batches.at(index);
                    const auto &batchSTime = batchInfo.start_time + eventsInfo->raw_start_time -
                                             _dataMagr->GetRawStartTimestamp();
                    const auto &batchETime = batchInfo.end_time + eventsInfo->raw_start_time -
                                             _dataMagr->GetRawStartTimestamp();

                    if ((batchSTime < st && batchETime < st) ||
                        (batchSTime > et && batchETime > et)) {
                        iter = tracking->erase(iter);
                        continue;
                    }

                    // const auto oldSize = batch.size();
                    EventTrackingFilter::FilterByTrackingLength(batch, TRACKING_LEN_PERCENT_THD);
                    EventTrackingFilter::FilterByTraceFittingSAC(batch, TRACKING_FIT_SAC_THD);
                    EventTrackingFilter::FilterByTrackingAge(batch, TRACKING_AGE_PERCENT_THD);
                    EventTrackingFilter::FilterByTrackingFreq(batch, TRACKING_FREQ_PERCENT_THD);
                    // spdlog::info(
                    //     "size before filtering: {}, size after filtering: {}, filtered: {}",
                    //     oldSize, batch.size(), oldSize - batch.size());

                    // draw
                    _viewer->ClearViewer(Viewer::VIEW_MAP);
                    _viewer->AddEventFeatTracking(batch, intri, static_cast<float>(batchSTime),
                                                  static_cast<float>(batchETime), Viewer::VIEW_MAP,
                                                  0.01, 20);
                    // auto iters = _dataMagr->ExtractEventDataPiece(topic, batchSTime, batchETime);
                    // _viewer->AddEventData(iters.first, iters.second, batchSTime,
                    // Viewer::VIEW_MAP, 0.01, 20);
                }
                bar->finish();
                // save tracking results
                eventFeatTrackingRes[topic] = *tracking;
                continue;
            }
        }
        // if tracking is not performed, we output raw event data for haste-powered feature
        // tracking
        /**
         * |--> 'outputSIter1'
         * |            |<- BATCH_TIME_WIN_THD ->|<- BATCH_TIME_WIN_THD ->|
         * ------------------------------------------------------------------
         * |<- BATCH_TIME_WIN_THD ->|<- BATCH_TIME_WIN_THD ->|
         * | data in this windown would be output for event-based feature tracking
         * |--> 'outputSIter2'
         */
        spdlog::info("saving event data of camera '{}' for haste-based feature tracking...", topic);
        const std::size_t EVENT_FRAME_NUM_THD = intri->imgHeight * intri->imgWidth / 5;
        SaveEventDataForFeatureTracking(topic, hasteWorkspace, 0.2, EVENT_FRAME_NUM_THD, 200);
    }
    cv::destroyAllWindows();
    _viewer->ClearViewer(Viewer::VIEW_MAP);
    if (eventFeatTrackingRes.size() != Configor::DataStream::EventTopics.size()) {
        throw Status(Status::FINE,
                     "files for haste-based feature tracking have been output to '{}', run "
                     "corresponding commands shell files to generate tracking results!",
                     Configor::DataStream::OutputPath + "/events/");
    }

    /**
     * Based on the results given by HASTE, we perform consistency detection based on the visual
     * model to estimate the camera rotation and remove bad tracking.
     */
    constexpr double REL_ROT_TIME_INTERVAL = 0.03 /* about 30 Hz */;
    constexpr double FMAT_THRESHOLD = 1.0;
    constexpr double ROT_ONLY_RANSAC_THRESHOLD = 1.0;
    for (const auto &[topic, tracking] : eventFeatTrackingRes) {
        // traces of all features
        std::map<FeatureTrackingTrace::Ptr, EventFeatTrackingVec> traceVec;
        for (const auto &[id, batch] : tracking) {
            for (const auto &[fId, featVec] : batch) {
                auto trace = FeatureTrackingTrace::CreateFrom(featVec);
                if (trace != nullptr) {
                    traceVec[trace] = featVec;
                }
            }
        }
        auto traceVecOldSize = traceVec.size();
        auto FindInRangeTrace = [&traceVec](double time) {
            std::map<FeatureTrackingTrace::Ptr, Eigen::Vector2d> inRangeTraceVec;
            for (const auto &[trace, featVec] : traceVec) {
                if (auto pos = trace->PositionAt(time); pos != std::nullopt) {
                    inRangeTraceVec[trace] = *pos;
                }
            }
            return inRangeTraceVec;
        };
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        // time from {t1} to {t2}, relative rotation from {t2} to {t1}
        RotationEstimator::RelRotationSequence relRotations;
        spdlog::info("perform rotation-only relative rotation estimation for '{}'...", topic);
        auto bar = std::make_shared<tqdm>();
        const int totalSize = static_cast<int>((et - st) / REL_ROT_TIME_INTERVAL) - 1;
        int barIndex = 0;
        for (double time = st; time < et;) {
            bar->progress(barIndex++, totalSize);
            const double t1 = time, t2 = time + REL_ROT_TIME_INTERVAL;
            time += REL_ROT_TIME_INTERVAL;

            const auto traceVec1 = FindInRangeTrace(t1), traceVec2 = FindInRangeTrace(t2);

            // trace, pos at t1, pos at t2
            std::map<FeatureTrackingTrace::Ptr, std::pair<Eigen::Vector2d, Eigen::Vector2d>>
                matchedTraceVec;
            for (const auto &[trace, pos1] : traceVec1) {
                auto iter = traceVec2.find(trace);
                if (iter != traceVec2.cend()) {
                    matchedTraceVec[trace] = {pos1, iter->second};
                }
            }
            // spdlog::info("matched feature count from '{:.3f}' to '{:.3f}': {}", t1, t2,
            //              matchedTraceVec.size());

            // using rotation-only relative rotation solver to estimate the relative rotation
            auto size = matchedTraceVec.size();
            std::vector<Eigen::Vector2d> featUndisto1, featUndisto2;
            featUndisto1.reserve(size), featUndisto2.reserve(size);
            std::map<int, FeatureTrackingTrace::Ptr> featIdxMap;
            int index = 0;
            for (const auto &[trace, posPair] : matchedTraceVec) {
                featUndisto1.push_back(posPair.first);
                featUndisto2.push_back(posPair.second);
                featIdxMap[index++] = trace;
            }
            auto resFMat =
                RotOnlyVisualOdometer::RejectUsingFMat(featUndisto1, featUndisto2, FMAT_THRESHOLD);
            if (resFMat.first) {
                // solving successful
                const auto &status = resFMat.second;
                std::vector<Eigen::Vector2d> featUndisto1Temp, featUndisto2Temp;
                std::map<int, FeatureTrackingTrace::Ptr> featIdxMapTemp;
                int inlierIdx = 0;
                for (int i = 0; i < static_cast<int>(status.size()); ++i) {
                    if (status.at(i)) {
                        // is an inlier
                        featUndisto1Temp.push_back(featUndisto1.at(i));
                        featUndisto2Temp.push_back(featUndisto2.at(i));
                        featIdxMapTemp[inlierIdx++] = featIdxMap.at(i);
                    } else {
                        // is an outlier
                        traceVec.erase(featIdxMap.at(i));
                    }
                }
                // spdlog::info("'FMat Rejection': inliers count: {}, outliers count: {}, total:
                // {}", featIdxMapTemp.size(), featIdxMap.size() - featIdxMapTemp.size(),
                //              featIdxMap.size());

                // draw match figure
                // cv::Mat matchImg(static_cast<int>(intri->imgHeight),
                //                  static_cast<int>(intri->imgWidth) * 2, CV_8UC3,
                //                  cv::Scalar(255, 255, 255));
                // matchImg.colRange(0.0, static_cast<int>(intri->imgWidth))
                //     .setTo(cv::Scalar(200, 200, 200));
                // for (int i = 0; i < static_cast<int>(featUndisto1.size()); ++i) {
                //     const Eigen::Vector2d &p1 = featUndisto1.at(i);
                //     Eigen::Vector2d p2 = featUndisto2.at(i) + Eigen::Vector2d(intri->imgWidth,
                //     0.0);
                //     // feature pair
                //     DrawKeypointOnCVMat(matchImg, p1, true, cv::Scalar(0, 0, 255));
                //     DrawKeypointOnCVMat(matchImg, p2, true, cv::Scalar(0, 0, 255));
                //     // line
                //     DrawLineOnCVMat(matchImg, p1, p2, cv::Scalar(0, 0, 255));
                // }
                // for (int i = 0; i < static_cast<int>(featUndisto1Temp.size()); ++i) {
                //     const Eigen::Vector2d &p1 = featUndisto1Temp.at(i);
                //     Eigen::Vector2d p2 =
                //         featUndisto2Temp.at(i) + Eigen::Vector2d(intri->imgWidth, 0.0);
                //     // feature pair
                //     DrawKeypointOnCVMat(matchImg, p1, true, cv::Scalar(0, 255, 0));
                //     DrawKeypointOnCVMat(matchImg, p2, true, cv::Scalar(0, 255, 0));
                //     // line
                //     DrawLineOnCVMat(matchImg, p1, p2, cv::Scalar(0, 255, 0));
                // }
                // cv::imshow("FMat Rejection", matchImg);
                // cv::waitKey(0);

                // assign
                featIdxMap = featIdxMapTemp;
                featUndisto1 = featUndisto1Temp;
                featUndisto2 = featUndisto2Temp;
            } else {
                continue;
            }
            auto resRotOnly = RotOnlyVisualOdometer::RelRotationRecovery(
                featUndisto1,                // features at t1
                featUndisto2,                // corresponding features at t2
                intri,                       // the camera intrinsics
                ROT_ONLY_RANSAC_THRESHOLD);  // the threshold
            if (!resRotOnly.second.empty()) {
                // solving successful
                relRotations.emplace_back(
                    t1, t2, Sophus::SO3d(Sophus::makeRotationMatrix(resRotOnly.first)));
                // find outliers
                std::set<int> inlierIdxSet;
                for (int idx : resRotOnly.second) {
                    inlierIdxSet.insert(idx);
                }
                for (const auto &[idx, trace] : featIdxMap) {
                    if (inlierIdxSet.find(idx) == inlierIdxSet.cend()) {
                        // is an outlier, we remove this trace as it's poor
                        traceVec.erase(trace);
                    }
                }
                // spdlog::info(
                //     "'Rotation-Only RANSAC': inliers count: {}, outliers count: {}, total: {}",
                //     inlierIdxSet.size(), featIdxMap.size() - inlierIdxSet.size(),
                //     featIdxMap.size());

                // draw match figure
                // cv::Mat matchImg(static_cast<int>(intri->imgHeight),
                //                  static_cast<int>(intri->imgWidth) * 2, CV_8UC3,
                //                  cv::Scalar(255, 255, 255));
                // matchImg.colRange(0.0, static_cast<int>(intri->imgWidth))
                //     .setTo(cv::Scalar(200, 200, 200));
                // for (int i = 0; i < static_cast<int>(featUndisto1.size()); ++i) {
                //     const Eigen::Vector2d &p1 = featUndisto1.at(i);
                //     Eigen::Vector2d p2 = featUndisto2.at(i) + Eigen::Vector2d(intri->imgWidth,
                //     0.0);
                //     // feature pair
                //     DrawKeypointOnCVMat(matchImg, p1, true, cv::Scalar(0, 0, 255));
                //     DrawKeypointOnCVMat(matchImg, p2, true, cv::Scalar(0, 0, 255));
                //     // line
                //     DrawLineOnCVMat(matchImg, p1, p2, cv::Scalar(0, 0, 255));
                // }
                // for (int idx : inlierIdxSet) {
                //     const Eigen::Vector2d &p1 = featUndisto1.at(idx);
                //     Eigen::Vector2d p2 =
                //         featUndisto2.at(idx) + Eigen::Vector2d(intri->imgWidth, 0.0);
                //     // feature pair
                //     DrawKeypointOnCVMat(matchImg, p1, true, cv::Scalar(0, 255, 0));
                //     DrawKeypointOnCVMat(matchImg, p2, true, cv::Scalar(0, 255, 0));
                //     // line
                //     DrawLineOnCVMat(matchImg, p1, p2, cv::Scalar(0, 255, 0));
                // }
                // cv::imshow("Rotation-Only RANSAC", matchImg);
                // cv::waitKey(0);
            }
        }
        bar->finish();
        auto traceVecNewSize = traceVec.size();
        spdlog::info("event trace count before rejection: {}, count after rejection: {}",
                     traceVecOldSize, traceVecNewSize);

        auto rotEstimator = RotationEstimator::Create();
        rotEstimator->Estimate(so3Spline, relRotations);
        if (rotEstimator->SolveStatus()) {
            // assign the estimated extrinsic rotation
            _parMagr->EXTRI.SO3_EsToBr.at(topic) = rotEstimator->GetSO3SensorToSpline();
        } else {
            throw Status(Status::ERROR,
                         "initialize rotation 'SO3_EsToBr' failed, this may be related to "
                         "insufficiently excited motion or bad event data streams.");
        }
        /**
         * once the extrinsic rotation is recovered, if time offset is also required, we
         * continue to recover it and refine extrineic rotation using
         * continuous-time-based alignment
         */
        if (Configor::Prior::OptTemporalParams) {
            auto estimator = Estimator::Create(_splines, _parMagr);

            auto optOption = OptOption::OPT_SO3_EsToBr | OptOption::OPT_TO_EsToBr;
            double TO_EsToBr = _parMagr->TEMPORAL.TO_EsToBr.at(topic);
            double weight = Configor::DataStream::EventTopics.at(topic).Weight;

            for (const auto &[lastTime, curTime, SO3_CurToLast] : relRotations) {
                // we throw the head and tail data as the rotations from the fitted SO3
                // Spline in that range are poor
                if (lastTime + TO_EsToBr < st || curTime + TO_EsToBr > et) {
                    continue;
                }
                estimator->AddHandEyeRotationAlignmentForEvent(
                    topic,          // the ros topic
                    lastTime,       // the time of start rotation stamped by the camera
                    curTime,        // the time of end rotation stamped by the camera
                    SO3_CurToLast,  // the relative rotation
                    optOption,      // the optimization option
                    weight          // the weight
                );
            }

            // we don't want to output the solving information
            auto optWithoutOutput =
                Estimator::DefaultSolverOptions(Configor::Preference::AvailableThreads(),
                                                false,  // do not output the solving information
                                                Configor::Preference::UseCudaInSolving);

            estimator->Solve(optWithoutOutput, _priori);
        }
        _viewer->UpdateSensorViewer();
    }
    _parMagr->ShowParamStatus();
    spdlog::warn("developing!!!");
    std::cin.get();
}
}  // namespace ns_ikalibr