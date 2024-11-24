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
#include "calib/estimator.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepEventInertialAlignLineBased() const {
    if (!Configor::IsEventIntegrated()) {
        return;
    }

    /**
     * estimate norm flows and recover extrinsic rotations and time offsets of cameras under the
     * pure rotation motion assumption
     */
    constexpr double TIME_SURFACE_DECAY_TIME = 0.03;
    std::map<std::string, std::list<std::list<NormFlow::Ptr>>> nfsForEventCams;
    for (const auto &[topic, eventMes] : _dataMagr->GetEventMeasurements()) {
        spdlog::info("perform norm flow estimation for event camera '{}'...", topic);
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        auto saeCreator = ActiveEventSurface::Create(intri, 0.01);
        double lastNfEventTime = eventMes.front()->GetTimestamp();
        auto &nfsCurCam = nfsForEventCams[topic];
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
                2,                        // window size to fit local planes
                4,                        // distance between neighbor norm flows
                0.9,                      // the ratio, for ransac and in-range candidates
                2E-3,  // the point to plane threshold in temporal domain, unit (s)
                2);    // ransac iteration count
            lastNfEventTime = res.timestamp;

            if (res.nfs.empty()) {
                continue;
            }

            nfsCurCam.push_back(res.nfs);

            cv::imshow("Time Surface & Norm Flow", res.Visualization(0.02));
            // _viewer->AddEventData(res.ActiveEvents(0.02), res.timestamp, Viewer::VIEW_MAP,
            //                       {0.01, 100});
            // _viewer->AddEventData(res.NormFlowEvents(), res.timestamp, Viewer::VIEW_MAP,
            //                       {0.01, 100}, ns_viewer::Colour::Green());
            // _viewer->ClearViewer(Viewer::VIEW_MAP);
            cv::waitKey(1);
        }
        bar->finish();
    }
    cv::destroyAllWindows();

    for (const auto &[topic, nfsList] : nfsForEventCams) {
        spdlog::info("init extrinsic rotations and time offsets for event camera '{}'...", topic);
        auto estimator = Estimator::Create(_splines, _parMagr);
        auto opt = OptOption::OPT_SO3_EsToBr;
        /**
         * Although the 'AddEventNormFlowRotConstraint' supports time delay estimation, it assumes
         * pure rotational motion. When the linear velocity is non-zero (i.e., additional
         * translation exists), model errors may be introduced, which could potentially be absorbed
         * into the time delay. Therefore, time delay is not optimized in this case.
         */
        // if (Configor::Prior::OptTemporalParams) {
        //     opt |= OptOption::OPT_TO_EsToBr;
        // }
        for (const auto &nfs : nfsList) {
            for (const auto &nf : nfs) {
                estimator->AddEventNormFlowRotConstraint(nf, topic, opt, 1.0);
            }
        }
        auto sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
        _viewer->UpdateSensorViewer();
    }
    _parMagr->ShowParamStatus();
    std::cin.get();
}
}  // namespace ns_ikalibr