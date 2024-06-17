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

#include "core/radar_velocity_sac.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

int RadarVelocitySacProblem::getSampleSize() const { return 3; }

bool RadarVelocitySacProblem::computeModelCoefficients(
    const std::vector<int> &indices, RadarVelocitySacProblem::model_t &outModel) const {
    std::vector<RadarTarget::Ptr> selectedTargets(indices.size());
    const auto &allTargets = _data->GetTargets();
    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        selectedTargets.at(i) = allTargets.at(indices.at(i));
    }
    outModel = RadarTargetArray(_data->GetTimestamp(), selectedTargets)
                   .RadarVelocityFromStaticTargetArray();
    return true;
}

void RadarVelocitySacProblem::getSelectedDistancesToModel(
    const RadarVelocitySacProblem::model_t &model,
    const std::vector<int> &indices,
    std::vector<double> &scores) const {
    scores.resize(indices.size());
    const auto &allTargets = _data->GetTargets();
    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        const auto &curTar = allTargets.at(indices.at(i));
        const Eigen::Vector3d &tPos = curTar->GetTargetXYZ();
        scores.at(i) = curTar->GetRadialVelocity() + tPos.dot(model) / tPos.norm();
    }
}

void RadarVelocitySacProblem::optimizeModelCoefficients(
    const std::vector<int> &inliers,
    const RadarVelocitySacProblem::model_t &model,
    RadarVelocitySacProblem::model_t &optimized_model) {
    computeModelCoefficients(inliers, optimized_model);
}
}  // namespace ns_ikalibr