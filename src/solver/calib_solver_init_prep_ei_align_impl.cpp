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
#include "core/feature_tracking.h"
#include "core/event_trace_sac.h"
#include "calib/estimator.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepEventInertialAlign() const {
    throw Status(
        Status::WARNING,
        "Although point-based optical flow event-inertial calibration has been developed "
        "in iKalibr, the poor results of event-based sparse optical flow estimation have led to "
        "unsatisfactory calibration outcomes, leading to the temporary suspension of "
        "support for this method. Accurate event-based sparse optical flow estimation is "
        "all it needs!");
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
    constexpr double TRACKING_FIT_SAC_THD = 3.0;
    constexpr double TRACKING_AGE_PERCENT_THD = 0.3;
    constexpr double TRACKING_FREQ_PERCENT_THD = 0.3;
    constexpr double BATCH_TIME_WIN_THD = 0.2;
    constexpr int HASTE_SEED_COUNT = 50;
    // topic, batch index, tracked features
    std::map<std::string, std::map<int, FeatureVecMap>> eventFeatTrackingRes;
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
            // todo: the output tracking results from HASTE should be distortion-free?
            auto tracking = HASTEDataIO::TryLoadHASTEResultsFromBinary(
                *eventsInfo, intri, _dataMagr->GetRawStartTimestamp());
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
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                    _viewer->ClearViewer(Viewer::VIEW_MAP);
                    _viewer->AddEventFeatTracking(batch, intri, static_cast<float>(batchSTime),
                                                  static_cast<float>(batchETime), Viewer::VIEW_MAP);
                    // auto iters = _dataMagr->ExtractEventDataPiece(topic, batchSTime, batchETime);
                    // _viewer->AddEventData(iters.first, iters.second, batchSTime,
                    // Viewer::VIEW_MAP, 0.01, 20);
                }
                bar->finish();
                // save tracking results
                eventFeatTrackingRes[topic] = *tracking;
                _viewer->ClearViewer(Viewer::VIEW_MAP);
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
        SaveEventDataForFeatureTracking(topic, hasteWorkspace, BATCH_TIME_WIN_THD,
                                        HASTE_SEED_COUNT);
    }
    cv::destroyAllWindows();
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
    constexpr double DISCRETE_TIME_INTERVAL = 0.03 /* about 30 Hz */;
    constexpr double FMAT_THRESHOLD = 1.0;
    constexpr double ROT_ONLY_RANSAC_THD = 1.0;
    std::map<std::string, std::set<FeatureTrackingCurve::Ptr>> eventTraceMap;
    for (const auto &[topic, tracking] : eventFeatTrackingRes) {
        // traces of all features
        auto &traceVec = eventTraceMap[topic];
        for (const auto &[id, batch] : tracking) {
            for (const auto &[fId, featVec] : batch) {
                auto trace = FeatureTrackingCurve::CreateFrom(featVec, true);
                if (trace != nullptr) {
                    traceVec.insert(trace);
                }
            }
        }
        auto traceVecOldSize = traceVec.size();
        auto FindInRangeTrace = [&traceVec](double time) {
            std::map<FeatureTrackingCurve::Ptr, Eigen::Vector2d> inRangeTraceVec;
            for (const auto &trace : traceVec) {
                if (auto pos = trace->PositionAt(time); pos != std::nullopt) {
                    inRangeTraceVec[trace] = *pos;
                }
            }
            return inRangeTraceVec;
        };
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        // time from {t1} to {t2}, relative rotation from {t2} to {t1}
        RotationEstimator::RelRotationSequence relRotations;
        auto rotEstimator = RotationEstimator::Create();
        spdlog::info("perform rotation-only relative rotation estimation for '{}'...", topic);
        auto bar = std::make_shared<tqdm>();
        const int totalSize = static_cast<int>((et - st) / DISCRETE_TIME_INTERVAL) - 1;
        int barIndex = 0;
        for (double time = st; time < et - DISCRETE_TIME_INTERVAL;) {
            bar->progress(barIndex++, totalSize);
            const double t1 = time, t2 = time + DISCRETE_TIME_INTERVAL;
            time += DISCRETE_TIME_INTERVAL;

            const auto traceVec1 = FindInRangeTrace(t1), traceVec2 = FindInRangeTrace(t2);

            // trace, pos at t1, pos at t2
            std::map<FeatureTrackingCurve::Ptr, std::pair<Eigen::Vector2d, Eigen::Vector2d>>
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
            std::map<int, FeatureTrackingCurve::Ptr> featIdxMap;
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
                std::map<int, FeatureTrackingCurve::Ptr> featIdxMapTemp;
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
                featUndisto1,          // features at t1
                featUndisto2,          // corresponding features at t2
                intri,                 // the camera intrinsics
                ROT_ONLY_RANSAC_THD);  // the threshold
            if (!resRotOnly.second.empty()) {
                // solving successful, find outliers
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

                relRotations.emplace_back(
                    t1, t2, Sophus::SO3d(Sophus::makeRotationMatrix(resRotOnly.first)));
                if (rotEstimator->SolveStatus() || (relRotations.size() < 50) ||
                    (relRotations.size() % 5 != 0)) {
                    continue;
                }
                rotEstimator->Estimate(so3Spline, relRotations);
                if (!rotEstimator->SolveStatus()) {
                    continue;
                }
                // assign the estimated extrinsic rotation
                _parMagr->EXTRI.SO3_EsToBr.at(topic) = rotEstimator->GetSO3SensorToSpline();

                /**
                 * once the extrinsic rotation is recovered, if time offset is also
                 * required, we continue to recover it and refine extrineic rotation using
                 * continuous-time-based alignment
                 */
                if (Configor::Prior::OptTemporalParams) {
                    auto estimator = Estimator::Create(_splines, _parMagr);

                    auto optOption = OptOption::OPT_SO3_EsToBr | OptOption::OPT_TO_EsToBr;
                    double TO_EsToBr = _parMagr->TEMPORAL.TO_EsToBr.at(topic);
                    double weight = Configor::DataStream::EventTopics.at(topic).Weight;

                    for (const auto &[lastTime, curTime, SO3_CurToLast] : relRotations) {
                        // we throw the head and tail data as the rotations from the fitted
                        // SO3 Spline in that range are poor
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
                    auto optWithoutOutput = Estimator::DefaultSolverOptions(
                        Configor::Preference::AvailableThreads(),
                        false,  // do not output the solving information
                        Configor::Preference::UseCudaInSolving);

                    estimator->Solve(optWithoutOutput, _priori);
                }
                _viewer->UpdateSensorViewer();

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

        if (!rotEstimator->SolveStatus()) {
            throw Status(Status::ERROR,
                         "initialize rotation 'SO3_EsToBr' failed, this may be related to "
                         "insufficiently excited motion or bad event data streams.");
        }
    }

    // save eventTraceMap to data manager
    for (const auto &[topic, traceSet] : eventTraceMap) {
        std::vector<FeatureTrackingCurve::Ptr> traceVec;
        traceVec.reserve(traceSet.size());
        for (const auto &trace : traceSet) {
            traceVec.push_back(trace);
        }
        _dataMagr->SetVisualFeatureTrackingCurve(topic, traceVec);
    }

#define USE_NEW_CAM_VEL_ESTIMATE 1
#if USE_NEW_CAM_VEL_ESTIMATE
    for (const auto &[topic, eventTrace] : eventTraceMap) {
        auto FindInRangeTrackedTripleFeatures = [&eventTrace](double time, double interval) {
            std::vector<std::array<Eigen::Vector2d, 3>> triplePosVec;
            const double t1 = time - interval;
            const double t2 = time;
            const double t3 = time + interval;

            for (const auto &trace : eventTrace) {
                auto pos1 = trace->PositionAt(t1);
                if (pos1 == std::nullopt) {
                    continue;
                }
                auto pos2 = trace->PositionAt(t2);
                if (pos2 == std::nullopt) {
                    continue;
                }
                auto pos3 = trace->PositionAt(t3);
                if (pos3 == std::nullopt) {
                    continue;
                }
                triplePosVec.emplace_back(std::array<Eigen::Vector2d, 3>{*pos1, *pos2, *pos3});
            }
            return triplePosVec;
        };
        std::vector<OpticalFlowTripleTrace::Ptr> traceVec;
        std::list<CameraFrame::Ptr> frameTripleList;
        for (double time = st + DISCRETE_TIME_INTERVAL; time < et - DISCRETE_TIME_INTERVAL;) {
            const double t0 = time - DISCRETE_TIME_INTERVAL;
            const double t1 = time;
            const double t2 = time + DISCRETE_TIME_INTERVAL;
            time += DISCRETE_TIME_INTERVAL;

            if (t0 < st || t2 > et) {
                continue;
            }

            auto triplePosVec = FindInRangeTrackedTripleFeatures(t1, DISCRETE_TIME_INTERVAL);

            if (frameTripleList.empty()) {
                frameTripleList.push_back(CameraFrame::Create(t0, cv::Mat(), cv::Mat(), 0));
                frameTripleList.push_back(CameraFrame::Create(t1, cv::Mat(), cv::Mat(), 1));
                frameTripleList.push_back(CameraFrame::Create(t2, cv::Mat(), cv::Mat(), 2));
            } else {
                frameTripleList.pop_front();
                auto id = frameTripleList.back()->GetId() + 1;
                frameTripleList.push_back(CameraFrame::Create(t2, cv::Mat(), cv::Mat(), id));
            }
            auto f0 = frameTripleList.front();
            auto f1 = *std::next(frameTripleList.cbegin());
            auto f2 = frameTripleList.back();
            for (const auto &ary : triplePosVec) {
                traceVec.push_back(OpticalFlowTripleTrace::Create({std::pair{f0, ary.at(0)},
                                                                   std::pair{f1, ary.at(1)},
                                                                   std::pair{f2, ary.at(2)}}));
            }
        }
        _dataMagr->SetVisualOpticalFlowTrace(topic, traceVec);
        spdlog::info("create optical flow triple trace count for '{}': {}", topic, traceVec.size());
    }

    /**
     * resort the optical flow correspondences based on the camera frame index
     */
    std::map<std::string, std::map<CameraFrame::Ptr, std::vector<OpticalFlowCorr::Ptr>>>
        opticalFlowInFrame;
    for (const auto &[topic, _] : Configor::DataStream::EventTopics) {
        const auto &traceVec = _dataMagr->GetVisualOpticalFlowTrace(topic);
        auto &curOpticalFlowInFrame = opticalFlowInFrame[topic];

        for (const auto &trace : traceVec) {
            auto midCamFrame = trace->GetMidCameraFrame();
            // the depth information stored in 'OpticalFlowCorr' is invalid currently
            const auto &corr = trace->CreateEventOpticalFlowCorr();
            curOpticalFlowInFrame[midCamFrame].emplace_back(corr);

    #define VISUALIZE_OPTICAL_FLOW_TRACE 0
    #if VISUALIZE_OPTICAL_FLOW_TRACE
            const double readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);
            Eigen::Vector2d vel = corr->MidPointVel(readout);  // the pixel velocity
            if (vel.norm() > 5.0 * Configor::Prior::LossForOpticalFlowFactor) {
                // just for visualization
                const auto &intri = _parMagr->INTRI.Camera.at(topic);
                auto h = static_cast<int>(intri->imgHeight), w = static_cast<int>(intri->imgWidth);
                auto movement = trace->GetTrace();
                std::array<cv::Scalar, 3> color = {
                    cv::Scalar{200, 200, 200}, {255, 255, 255}, {200, 200, 200}};
                for (int i = 0; i != movement.size(); ++i) {
                    auto &mv = movement.at(i);
                    auto cImg = cv::Mat(h, w, CV_8UC3, color.at(i));
                    mv.first = CameraFrame::Create(mv.first->GetTimestamp(), cv::Mat(), cImg,
                                                   mv.first->GetId());
                }
                // show the visual pixel trace image (features, mid-point pixel velocity)
                auto img = OpticalFlowTripleTrace::Create(movement)->CreateOpticalFlowMat(
                    _parMagr->INTRI.Camera.at(topic), corr->MidPointVel(readout));
                Eigen::Vector2d mp = corr->MidPoint();
                const auto filename = fmt::format(
                    "{}/{}-{}-{}.png", Configor::DataStream::DebugPath, midCamFrame->GetId(),
                    static_cast<int>(mp(0)), static_cast<int>(mp(1)));
                // save this image to disk
                // cv::imwrite(filename, img);
                cv::imshow("img", img);
                cv::waitKey(0);
            }
    #endif
    #undef VISUALIZE_OPTICAL_FLOW_TRACE
        }
    }

    /**
     * for each event camera (virtual) frame of each topic, we estimate the linear velocity, and
     * store them in 'eventBodyFrameVelDirs' [topic, camera frame, body-frame velocity]
     */
    auto &eventBodyFrameVelDirs = _initAsset->eventBodyFrameVelDirs;
    for (const auto &[topic, _] : Configor::DataStream::EventTopics) {
        spdlog::info("estimate event-derived linear velocities for '{}'...", topic);
        const auto &readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);
        const double TO_EsToBr = _parMagr->TEMPORAL.TO_EsToBr.at(topic);

        auto bar = std::make_shared<tqdm>();
        const auto &curOpticalFlowInFrame = opticalFlowInFrame.at(topic);
        auto totalSize = static_cast<int>(curOpticalFlowInFrame.size());
        int curIdx = 0;
        for (const auto &[frame, ofVec] : curOpticalFlowInFrame) {
            bar->progress(curIdx++, totalSize);
            const double timeByBr = frame->GetTimestamp() + TO_EsToBr;
            // at least two measurements are required, here we up the ante
            if (timeByBr < st || timeByBr > et || ofVec.size() < 5) {
                continue;
            }
            // if the camera is moving too slow, we do not estimate its velocity direction
            double avgPixelVel = 0.0;
            for (const auto &corr : ofVec) {
                avgPixelVel += corr->MidPointVel(readout).norm();
            }
            avgPixelVel /= static_cast<double>(ofVec.size());
            if (avgPixelVel < Configor::Prior::LossForOpticalFlowFactor) {
                continue;
            }
            auto estimator = Estimator::Create(_splines, _parMagr);
            Eigen::Vector3d velDir(0.0, 0.0, 1.0f);
            for (auto &corr : ofVec) {
                // we initialize the depth as 1.0, this value would be optimized in estimator
                corr->depth = 1.0;
                estimator->AddVisualVelocityDepthFactorForEvent(
                    &velDir,  // the direction of linear velocity to be estimated
                    corr,     // the optical flow correspondence
                    topic,    // rostopic
                    1.0,      // weight
                    true,     // estimate the depth information
                    true);    // only estimate the direction of the linear velocity
            }
            // we don't want to output the solving information
            auto optWithoutOutput =
                Estimator::DefaultSolverOptions(Configor::Preference::AvailableThreads(),
                                                false,   // do not output the solving information
                                                false);  // we do not use cuda solving here
            auto sum = estimator->Solve(optWithoutOutput, this->_priori);

            eventBodyFrameVelDirs[topic].emplace_back(frame->GetTimestamp(), velDir);
        }
        bar->finish();
        /**
         * the obtained camera-frame velocities may not time-ordered, thus we sort these quantities
         * based on their timestamps
         */
        auto &curVelDirs = eventBodyFrameVelDirs.at(topic);
        std::sort(curVelDirs.begin(), curVelDirs.end(),
                  [](const auto &p1, const auto &p2) { return p1.first < p2.first; });
    }
#else
    /**
     * create discrete optical flow trace using fitted event trace
     * |    [---][---][---][---][---]   [---][---]-  [---][---][---][---][---]-  |
     * |[---][---]-   [---][---]--   [---][---][---][---]  [---][---][---][---]- |
     * |  [---][---][---]       [---][---][---]  [---][---][---][---][---]--     |
     */
    // constexpr double TRACKING_FINAL_FIT_SAC_THD = 2.0;
    // constexpr double TRACKING_FINAL_AGE_THD = DISCRETE_TIME_INTERVAL * 3 /*sed*/;
    // std::map<std::string, std::map<FeatureTrackingTrace::Ptr, EventFeatTrackingVec>>
    //     eventTraceMapFiltered;
    for (const auto &[topic, eventTrace] : eventTraceMap) {
        /**
         * At first we wanted to filter the obtained event trace again, but it is not necessary.
         */
        // auto &eventTraceFiltered = eventTraceMapFiltered[topic];
        // for (const auto &[trace, trackList] : eventTrace) {
        //     auto res = EventTrackingTraceSacProblem::EventTrackingTraceSac(
        //         trackList, TRACKING_FINAL_FIT_SAC_THD);
        //     if (res.first != nullptr) {
        //         double age = res.first->eTime - res.first->sTime;
        //         if (age > TRACKING_FINAL_AGE_THD) {
        //             eventTraceFiltered[res.first] = res.second;
        //         }
        //     }
        // }
        // spdlog::info("event trace count before filtered for '{}': {}, count after filtered: {}",
        //              topic, eventTrace.size(), eventTraceFiltered.size());

        RotOnlyVisualOdometer::FeatTrackingInfo trackInfoList;
        ns_veta::IndexT index = 0;
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        for (const auto &trace : eventTrace) {
            std::list<std::pair<CameraFramePtr, Feature>> featList;
            for (auto t = trace->sTime; t < trace->eTime;) {
                if (auto up = trace->PositionAt(t); up != std::nullopt) {
                    Eigen::Vector2d rp = intri->GetDistoPixel(*up);
                    Feature feat(cv::Point2d(rp(0), rp(1)), cv::Point2d((*up)(0), (*up)(1)));
                    // we create a fake camera frame here to store the timestamp
                    featList.emplace_back(CameraFrame::Create(t), feat);
                }
                t += DISCRETE_TIME_INTERVAL;
            }
            trackInfoList[index++] = featList;
        }
        // store
        _dataMagr->SetVisualOpticalFlowTrace(
            // ros topic of this camera
            topic,
            // create the optical flow trace
            CreateOpticalFlowTrace({trackInfoList}, 3));
        spdlog::info("create 'OpticalFlowTripleTrace' count for event camera '{}': {}", topic,
                     _dataMagr->GetVisualOpticalFlowTrace(topic).size());
    }

    /**
     * create discrete 'OpticalFlowTripleTrace' to estimate discrete camera veloicty
     * |   -|----|----|    |----|----|----|-   | ---|----|--- |    | ---|----|----|    |
     * |----|--- |  --|----|----|----|-   |----|----|----|-   |----|----|----|   -|--- |
     * |  --|----|--- |   -|----|----|----|-   |----|----|----|----|  --|----|----|    |
     */
    auto &eventBodyFrameVelDirs = _initAsset->eventBodyFrameVelDirs;
    for (const auto &[topic, trackList] : eventTraceMap) {
        // depth, position, velocity
        using DepthPosVelTuple = std::tuple<double, Eigen::Vector2d, Eigen::Vector2d>;
        std::list<std::pair<double, std::vector<DepthPosVelTuple>>> ofsPerStampList;
        for (double time = st; time < et;) {
            std::vector<DepthPosVelTuple> inRangePosVelVec;
            for (const auto &trace : trackList) {
                auto pos = trace->PositionAt(time);
                auto vel = trace->VelocityAt(time);
                if (pos != std::nullopt && vel != std::nullopt) {
                    inRangePosVelVec.emplace_back(-1.0, *pos, *vel);
                }
    #define VISUALIZE_OPTICAL_FLOW_TRACE 0
    #if VISUALIZE_OPTICAL_FLOW_TRACE
                auto timeBefore = time - DISCRETE_TIME_INTERVAL;
                auto posBefore = trace->PositionAt(timeBefore);
                auto timeAfter = time + DISCRETE_TIME_INTERVAL;
                auto posBack = trace->PositionAt(timeAfter);
                if (posBefore != std::nullopt && pos != std::nullopt && vel != std::nullopt &&
                    posBack != std::nullopt &&
                    vel->norm() > 5.0 * Configor::Prior::LossForOpticalFlowFactor) {
                    const auto &intri = _parMagr->INTRI.Camera.at(topic);
                    auto h = static_cast<int>(intri->imgHeight);
                    auto w = static_cast<int>(intri->imgWidth);
                    cv::Mat img(h, w, CV_8UC3, cv::Scalar(255, 255, 255));
                    auto frameBefore = CameraFrame::Create(timeBefore, cv::Mat(), img.clone(), 0);
                    auto frame = CameraFrame::Create(time, cv::Mat(), img.clone(), 0);
                    auto frameBack = CameraFrame::Create(timeAfter, cv::Mat(), img.clone(), 0);
                    std::array<std::pair<CameraFrame::Ptr, Eigen::Vector2d>, 3> movement;
                    movement.at(0) = std::pair{frameBefore, *posBefore};
                    movement.at(1) = std::pair{frame, *pos};
                    movement.at(2) = std::pair{frameBack, *posBack};
                    auto of = OpticalFlowTripleTrace::Create(movement);
                    auto mat = of->CreateOpticalFlowMat(intri, *vel);
                    cv::imshow("mat", mat);
                    cv::waitKey(0);
                }
    #endif
    #undef VISUALIZE_OPTICAL_FLOW_TRACE
            }
            ofsPerStampList.emplace_back(time, inRangePosVelVec);
            time += DISCRETE_TIME_INTERVAL;
        }

        spdlog::info("estimate event-derived linear velocities for '{}'...", topic);
        auto bar = std::make_shared<tqdm>();
        auto totalSize = static_cast<int>(ofsPerStampList.size());
        int curIdx = 0;
        const double TO_EsToBr = _parMagr->TEMPORAL.TO_EsToBr.at(topic);
        for (auto &[timeByCam, dpvTupleVec] : ofsPerStampList) {
            bar->progress(curIdx++, totalSize);
            const double timeByBr = timeByCam + TO_EsToBr;
            // at least two measurements are required, here we up the ante
            if (timeByBr < st || timeByBr > et || dpvTupleVec.size() < 5) {
                continue;
            }
            // if the camera is moving too slow, we do not estimate its velocity direction
            double avgPixelVel = 0.0;
            for (const auto &val : dpvTupleVec) {
                avgPixelVel += std::get<2>(val).norm();
            }
            avgPixelVel /= static_cast<double>(dpvTupleVec.size());
            if (avgPixelVel < Configor::Prior::LossForOpticalFlowFactor) {
                continue;
            }
            auto estimator = Estimator::Create(_splines, _parMagr);
            Eigen::Vector3d velDir(0.0, 0.0, 1.0f);

            for (auto &[depth, position, velocity] : dpvTupleVec) {
                // we initialize the depth as 1.0, this value would be optimized in estimator
                depth = 1.0;
                estimator->AddVisualVelocityDepthFactorForEvent(
                    topic,      // ros topic for event camera
                    &velDir,    // the direction of linear velocity to be estimated
                    timeByCam,  // time stamped by event camera, we don't consider readout time here
                    position,   // the 2d pixel position of this feature
                    velocity,   // the 2d pixel velocity of this feature
                    &depth,     // the optical flow correspondence
                    1.0,        // weight
                    true,       // estimate the depth information
                    true);      // only estimate the direction of the linear velocity
            }
            // we don't want to output the solving information
            auto optWithoutOutput =
                Estimator::DefaultSolverOptions(Configor::Preference::AvailableThreads(),
                                                false,   // do not output the solving information
                                                false);  // we do not use cuda solving here
            auto sum = estimator->Solve(optWithoutOutput, this->_priori);

            eventBodyFrameVelDirs[topic].emplace_back(timeByCam, velDir);
        }
        bar->finish();
    }
#endif
}
}  // namespace ns_ikalibr