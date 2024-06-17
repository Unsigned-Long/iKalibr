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

#include "sensor/radar.h"
#include "ctraj/utils/utils.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
// access
[[nodiscard]] double RadarTarget::GetTimestamp() const { return _timestamp; }

void RadarTarget::SetTimestamp(double timestamp) { _timestamp = timestamp; }

const Eigen::Vector3d &RadarTarget::GetTargetXYZ() const { return _target; }

Eigen::Vector3d RadarTarget::GetTargetRTP() const { return ns_ctraj::XYZtoRTP<double>(_target); }

double RadarTarget::GetRadialVelocity() const { return _radialVel; }

RadarTarget::Ptr RadarTarget::Create(double timestamp, const Eigen::Vector4d &rawMes) {
    return std::make_shared<RadarTarget>(timestamp, rawMes);
}

RadarTarget::RadarTarget(double timestamp, const Eigen::Vector4d &rawMes)
    : _timestamp(timestamp),
      _target(ns_ctraj::RTPtoXYZ<double>({rawMes(0), rawMes(1), rawMes(2)})),
      _radialVel(rawMes(3)),
      _range(rawMes(0)),
      _invRange(1.0 / rawMes(0)) {}

RadarTarget::Ptr RadarTarget::Create(double timestamp,
                                     const Eigen::Vector3d &target,
                                     double radialVel) {
    return std::make_shared<RadarTarget>(timestamp, target, radialVel);
}

RadarTarget::RadarTarget(double timestamp, Eigen::Vector3d target, double radialVel)
    : _timestamp(timestamp),
      _target(std::move(target)),
      _radialVel(radialVel),
      _range(_target.norm()),
      _invRange(1.0 / _target.norm()) {}

double RadarTarget::GetRange() const { return _range; }

double RadarTarget::GetInvRange() const { return _invRange; }

RadarTargetArray::RadarTargetArray(double timestamp, const std::vector<RadarTarget::Ptr> &targets)
    : _timestamp(timestamp),
      _targets(targets) {}

RadarTargetArray::Ptr RadarTargetArray::Create(double timestamp,
                                               const std::vector<RadarTarget::Ptr> &targets) {
    return std::make_shared<RadarTargetArray>(timestamp, targets);
}

double RadarTargetArray::GetTimestamp() const { return _timestamp; }

const std::vector<RadarTarget::Ptr> &RadarTargetArray::GetTargets() const { return _targets; }

void RadarTargetArray::SetTimestamp(double timestamp) { _timestamp = timestamp; }

bool RadarTargetArray::SaveTargetArraysToDisk(const std::string &filename,
                                              const std::vector<RadarTargetArray::Ptr> &arrays,
                                              int precision) {
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(precision);
    cereal::JSONOutputArchive ar(file);
    ar(cereal::make_nvp("radar_arrays", arrays));
    return true;
}

Eigen::Vector3d RadarTargetArray::RadarVelocityFromStaticTargetArray() {
    Eigen::VectorXd lVec(_targets.size());
    Eigen::MatrixXd BMat(_targets.size(), 3);
    for (int i = 0; i < static_cast<int>(_targets.size()); ++i) {
        const auto &tar = _targets.at(i);
        lVec(i) = tar->GetRadialVelocity() * tar->GetTargetXYZ().norm();
        BMat.block<1, 3>(i, 0) = -tar->GetTargetXYZ().transpose();
    }
    Eigen::Vector3d xVec = (BMat.transpose() * BMat).inverse() * BMat.transpose() * lVec;
    return xVec;
}
}  // namespace ns_ikalibr
