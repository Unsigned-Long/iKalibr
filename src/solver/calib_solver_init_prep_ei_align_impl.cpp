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
                for (auto &[index, batch] : *tracking) {
                    // aligned time (start and end)
                    const auto &batchInfo = eventsInfo->batches.at(index);
                    const auto &batchSTime = batchInfo.start_time + eventsInfo->raw_start_time -
                                             _dataMagr->GetRawStartTimestamp();
                    const auto &batchETime = batchInfo.end_time + eventsInfo->raw_start_time -
                                             _dataMagr->GetRawStartTimestamp();

                    // const auto oldSize = batch.size();
                    EventTrackingFilter::FilterByTrackingLength(batch, 0.1 /*percent*/);
                    EventTrackingFilter::FilterByTraceFittingSAC(batch, 5.0 /*pixel*/);
                    EventTrackingFilter::FilterByTrackingAge(batch, 0.1 /*percent*/);
                    EventTrackingFilter::FilterByTrackingFreq(batch, 0.1 /*percent*/);
                    // spdlog::info(
                    //     "size before filtering: {}, size after filtering: {}, filtered: {}",
                    //     oldSize, batch.size(), oldSize - batch.size());

                    // draw
                    _viewer->ClearViewer(Viewer::VIEW_MAP);
                    _viewer->AddEventFeatTracking(batch, intri, batchSTime, batchETime,
                                                  Viewer::VIEW_MAP, 0.01, 20);
                    // auto iters = _dataMagr->ExtractEventDataPiece(topic, batchSTime, batchETime);
                    // _viewer->AddEventData(iters.first, iters.second, batchSTime,
                    // Viewer::VIEW_MAP, 0.01, 20);
                }
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
        auto FindInRangeTrace = [&traceVec](double time) {
            std::map<FeatureTrackingTrace::Ptr, Eigen::Vector2d> inRangeTraceVec;
            for (const auto &[trace, featVec] : traceVec) {
                if (auto pos = trace->PositionAt(time); pos != std::nullopt) {
                    inRangeTraceVec[trace] = *pos;
                }
            }
            return inRangeTraceVec;
        };
        for (double time = st; time < et;) {
            const double t1 = time, t2 = time + REL_ROT_TIME_INTERVAL;
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
            spdlog::info("matched feature count from '{:.3f}' to '{:.3f}': {}", t1, t2,
                         matchedTraceVec.size());

            // todo: using rotation-only relative rotation solver to estimate the relative rotation
            std::cin.get();

            time += REL_ROT_TIME_INTERVAL;
        }
    }

    spdlog::warn("developing!!!");
    std::cin.get();
}
}  // namespace ns_ikalibr