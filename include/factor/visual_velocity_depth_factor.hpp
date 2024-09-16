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

#ifndef IKALIBR_VISUAL_VELOCITY_DEPTH_FACTOR_HPP
#define IKALIBR_VISUAL_VELOCITY_DEPTH_FACTOR_HPP

#include "util/utils.h"
#include "factor/data_correspondence.h"
#include "ceres/dynamic_autodiff_cost_function.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct VisualVelocityDepthFactor {
private:
    Eigen::Vector2d point;
    Eigen::Vector2d pointVel;
    const Eigen::Vector3d ANG_VEL_CmToWInCm;
    double fx, fy;
    double cx, cy;

    double weight;

public:
    explicit VisualVelocityDepthFactor(Eigen::Vector2d point,
                                       Eigen::Vector2d pointVel,
                                       Eigen::Vector3d camFrameAngVel,
                                       const ns_veta::PinholeIntrinsic::Ptr &intri,
                                       double weight)

        : point(std::move(point)),
          pointVel(std::move(pointVel)),
          ANG_VEL_CmToWInCm(std::move(camFrameAngVel)),
          weight(weight) {
        fx = intri->FocalX();
        fy = intri->FocalY();
        cx = intri->PrincipalPoint()(0);
        cy = intri->PrincipalPoint()(1);
    }

    static auto Create(const Eigen::Vector2d &point,
                       const Eigen::Vector2d &pointVel,
                       const Eigen::Vector3d &camFrameAngVel,
                       const ns_veta::PinholeIntrinsic::Ptr &intri,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<VisualVelocityDepthFactor>(
            new VisualVelocityDepthFactor(point, pointVel, camFrameAngVel, intri, weight));
    }

    static std::size_t TypeHashCode() { return typeid(VisualVelocityDepthFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ LIN_VEL_CmToWInCm | DEPTH ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        // get value
        Eigen::Map<const Eigen::Vector3<T>> LIN_VEL_CmToWInCm(sKnots[0]);
        T DEPTH = sKnots[1][0];

        Eigen::Matrix<double, 2, 3> subAMat, subBMat;
        OpticalFlowCorr::SubMats<double>(&fx, &fy, &cx, &cy, point, &subAMat, &subBMat);

        Eigen::Vector2<T> pred =
            (1.0 / DEPTH) * subAMat * LIN_VEL_CmToWInCm + subBMat * ANG_VEL_CmToWInCm;

        Eigen::Map<Eigen::Vector2<T>> residuals(sResiduals);
        residuals = T(weight) * (pred - pointVel);

        return true;
    }
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_VISUAL_VELOCITY_DEPTH_FACTOR_HPP
