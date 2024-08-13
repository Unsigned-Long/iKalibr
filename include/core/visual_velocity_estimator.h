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

#ifndef IKALIBR_VISUAL_VELOCITY_ESTIMATOR_H
#define IKALIBR_VISUAL_VELOCITY_ESTIMATOR_H

#include "config/configor.h"
#include "ctraj/core/spline_bundle.h"
#include "ctraj/core/pose.hpp"
#include "veta/camera/pinhole.h"
#include "opencv2/core.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;

class VisualVelocityEstimator {
public:
    using Ptr = std::shared_ptr<VisualVelocityEstimator>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;
    using So3SplineType = SplineBundleType::So3SplineType;
    using RotationSequence = std::vector<std::pair<double, Sophus::SO3d>>;

protected:
    // pixel, velocity, depth
    std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>> _dynamics;
    ns_veta::PinholeIntrinsic::Ptr _intri;

public:
    explicit VisualVelocityEstimator(
        const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>> &dynamics,
        ns_veta::PinholeIntrinsic::Ptr intri);

    static Ptr Create(
        const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>> &dynamics,
        const ns_veta::PinholeIntrinsic::Ptr &intri);

    [[nodiscard]] std::optional<Eigen::Vector3d> Estimate(double timeByBr,
                                                          const So3SplineType &spline,
                                                          const Sophus::SO3d &SO3_DnToBr) const;

    static cv::Mat DrawVisualVelocityMat(
        const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>> &dynamics,
        const ns_veta::PinholeIntrinsic::Ptr &intri,
        double timeByBr,
        const So3SplineType &spline,
        const Sophus::SO3d &SO3_DnToBr,
        const Eigen::Vector3d &LIN_VEL_DnToWInDn,
        const CameraFramePtr &frame,
        double factor);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_VISUAL_VELOCITY_ESTIMATOR_H
