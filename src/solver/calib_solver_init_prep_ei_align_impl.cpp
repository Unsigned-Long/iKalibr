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

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepEventInertialAlign() const {
    if (!Configor::IsEventIntegrated()) {
        return;
    }
    // const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    // const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    /**
     * we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are
     * poor
     */
    // const double st = std::max(so3Spline.MinTime(), scaleSpline.MinTime()) +  // the max as start
    //                   Configor::Prior::TimeOffsetPadding;
    // const double et = std::min(so3Spline.MaxTime(), scaleSpline.MaxTime()) -  // the min as end
    //                   Configor::Prior::TimeOffsetPadding;

    /**
     * we first perform event-based feature tracking.
     * We first dedistort the events, then export them to files and use third-party software for
     * feature tracking.
     */
    int needFeatureTrackingCount = 0;
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

                    // draw
                    {
                        _viewer->ClearViewer(Viewer::VIEW_MAP);
                        _viewer->AddHASTETracking(batch, intri, batchSTime, batchETime,
                                                  Viewer::VIEW_MAP, 0.01, 20);
                        auto iters =
                            _dataMagr->ExtractEventDataPiece(topic, batchSTime, batchETime);
                        _viewer->AddEventData(iters.first, iters.second, batchSTime,
                                              Viewer::VIEW_MAP, 0.01, 20);
                        std::cin.get();
                    }

                    // todo: filter raw tracking results, select good ones
                    HASTEDataIO::FilterResultsByTrackingLength(batch, 0.1);

                    // draw

                    _viewer->ClearViewer(Viewer::VIEW_MAP);
                    _viewer->AddHASTETracking(batch, intri, batchSTime, batchETime,
                                              Viewer::VIEW_MAP, 0.01, 20);
                    auto iters = _dataMagr->ExtractEventDataPiece(topic, batchSTime, batchETime);
                    _viewer->AddEventData(iters.first, iters.second, batchSTime, Viewer::VIEW_MAP,
                                          0.01, 20);
                    std::cin.get();
                }
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
        SaveEventDataForFeatureTracking(topic, hasteWorkspace, 0.2 /*sed*/);
        ++needFeatureTrackingCount;
    }
    cv::destroyAllWindows();
    if (needFeatureTrackingCount != 0) {
        throw Status(Status::FINE,
                     "files for haste-based feature tracking have been output to '{}', run "
                     "corresponding commands shell files to generate tracking results!",
                     Configor::DataStream::OutputPath + "/events/");
    }
    spdlog::warn("developing!!!");
    std::cin.get();
}
}  // namespace ns_ikalibr