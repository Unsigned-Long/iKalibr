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

#include "calib/ceres_callback.h"
#include "viewer/viewer.h"
#include "calib/calib_param_manager.h"
#include "spdlog/spdlog.h"

namespace ns_ikalibr{

// ------------------
// CeresDebugCallBack
// ------------------

CeresDebugCallBack::CeresDebugCallBack(CalibParamManager::Ptr calibParamManager)
    : _parMagr(std::move(calibParamManager)),
      _outputDir(Configor::DataStream::OutputPath + "/iteration/epoch"),
      _idx(0) {
    if (std::filesystem::exists(_outputDir)) {
        std::filesystem::remove_all(_outputDir);
    }

    if (!std::filesystem::create_directories(_outputDir)) {
        spdlog::warn("create directory failed: '{}'", _outputDir);
    } else {
        _iterInfoFile = std::ofstream(_outputDir + "/epoch_info.csv", std::ios::out);
        _iterInfoFile << "cost,gradient,tr_radius(1/lambda)" << std::endl;
    }
}

ceres::CallbackReturnType CeresDebugCallBack::operator()(const ceres::IterationSummary &summary) {
    if (std::filesystem::exists(_outputDir)) {
        // save param
        const std::string paramFilename = _outputDir + "/ikalibr_param_" + std::to_string(_idx) +
                                          ns_ikalibr::Configor::GetFormatExtension();
        _parMagr->Save(paramFilename, ns_ikalibr::Configor::Preference::OutputDataFormat);

        // save iter info
        _iterInfoFile << _idx << ',' << summary.cost << ',' << summary.gradient_norm << ','
                      << summary.trust_region_radius << std::endl;

        ++_idx;
    }
    return ceres::SOLVER_CONTINUE;
}

CeresDebugCallBack::~CeresDebugCallBack() { _iterInfoFile.close(); }

// -------------------
// CeresViewerCallBack
// -------------------
CeresViewerCallBack::CeresViewerCallBack(Viewer::Ptr viewer)
    : _viewer(std::move(viewer)) {}

ceres::CallbackReturnType CeresViewerCallBack::operator()(const ceres::IterationSummary &summary) {
    _viewer->UpdateSensorViewer().UpdateSplineViewer();
    return ceres::CallbackReturnType::SOLVER_CONTINUE;
}

}
