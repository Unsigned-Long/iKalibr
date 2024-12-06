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
#include "core/ev_circle_tracking.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepEventInertialAlignCircleBased() const {
    if (!Configor::IsEventIntegrated()) {
        return;
    }

    /**
     * estimate norm flows and recover extrinsic rotations and time offsets of cameras under the
     * pure rotation motion assumption
     */
    constexpr double TIME_SURFACE_DECAY_TIME = 0.01;
    for (const auto &[topic, eventMes] : _dataMagr->GetEventMeasurements()) {
        spdlog::info("perform norm flow estimation for event camera '{}'...", topic);
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        auto saeCreator = ActiveEventSurface::Create(intri, 0.01);
        double lastNfEventTime = eventMes.front()->GetTimestamp();
        auto bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(eventMes.size()); i++) {
            bar->progress(i, eventMes.size());

            const auto &eventAry = eventMes.at(i);
            saeCreator->GrabEvent(eventAry);

            if (saeCreator->GetTimeLatest() - eventMes.front()->GetTimestamp() < 0.05 ||
                saeCreator->GetTimeLatest() - lastNfEventTime < TIME_SURFACE_DECAY_TIME) {
                continue;
            }
            // estimate norm flows
            auto nfCreator = EventNormFlow(saeCreator);
            auto res = nfCreator.ExtractNormFlows(
                TIME_SURFACE_DECAY_TIME,  // decay seconds for time surface
                1,                        // window size to fit local planes
                1,                        // distance between neighbor norm flows
                0.8,                      // the ratio, for ransac and in-range candidates
                false,                    // whether undistort images
                2E-3,  // the point to plane threshold in temporal domain, unit (s)
                3);    // ransac iteration count
            lastNfEventTime = res->timestamp;

            if (res->nfs.empty()) {
                continue;
            }

            EventCircleTracking::Create()->Process(res);

            cv::imshow("Time Surface & Norm Flow", res->Visualization(TIME_SURFACE_DECAY_TIME));
            cv::waitKey(0);
        }
        bar->finish();
    }
    cv::destroyAllWindows();

    _parMagr->ShowParamStatus();
    std::cin.get();
}
}  // namespace ns_ikalibr