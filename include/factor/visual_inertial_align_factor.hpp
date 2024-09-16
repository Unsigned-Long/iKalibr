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

#ifndef IKALIBR_VISUAL_INERTIAL_ALIGN_FACTOR_HPP
#define IKALIBR_VISUAL_INERTIAL_ALIGN_FACTOR_HPP

#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/core/spline_bundle.h"
#include "ctraj/core/pose.hpp"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"
#include "config/configor.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

template <int Order>
struct VisualInertialAlignHelper {
public:
    using Ptr = std::shared_ptr<VisualInertialAlignHelper>;
    using SplineBundleType = ns_ctraj::SplineBundle<Order>;
    using So3SplineType = typename SplineBundleType::So3SplineType;

public:
    double dt, dt2;
    Eigen::Vector3d velVec, posVec;
    Eigen::Matrix3d velMat, posMat;

    Eigen::Matrix3d sCMat, eCMat, diffCMat;

    Eigen::Vector3d diffDMat;

    Eigen::Matrix3d diffEMat;

    Sophus::SO3d mtSo3BrToBr0;

public:
    VisualInertialAlignHelper(const So3SplineType &so3Spline,
                              const ns_ctraj::Posed &sPose,
                              const ns_ctraj::Posed &ePose,
                              double mapTime,
                              double TO_CmToBr,
                              const std::pair<Eigen::Vector3d, Eigen::Matrix3d> &velVecMat,
                              const std::pair<Eigen::Vector3d, Eigen::Matrix3d> &posVecMat) {
        dt = ePose.timeStamp - sPose.timeStamp;
        dt2 = dt * dt;

        velVec = velVecMat.first;
        velMat = velVecMat.second;

        posVec = posVecMat.first;
        posMat = posVecMat.second;

        sCMat = CMat(so3Spline, sPose.timeStamp + TO_CmToBr);
        eCMat = CMat(so3Spline, ePose.timeStamp + TO_CmToBr);
        diffCMat = eCMat - sCMat;

        diffDMat = ePose.t - sPose.t;

        diffEMat = so3Spline.Evaluate(ePose.timeStamp + TO_CmToBr).matrix() -
                   so3Spline.Evaluate(sPose.timeStamp + TO_CmToBr).matrix();

        mtSo3BrToBr0 = so3Spline.Evaluate(mapTime + TO_CmToBr);
    }

protected:
    static Eigen::Matrix3d CMat(const So3SplineType &so3Spline, double t) {
        auto so3 = so3Spline.Evaluate(t);
        auto angVelInW = so3 * so3Spline.VelocityBody(t);
        return Sophus::SO3d::hat(angVelInW) * so3.matrix();
    }
};

extern template struct VisualInertialAlignHelper<Configor::Prior::SplineOrder>;

template <int Order>
struct VisualInertialAlignFactor {
private:
    VisualInertialAlignHelper<Order> helper;
    double weight;

public:
    VisualInertialAlignFactor(const VisualInertialAlignHelper<Order> &helper, double weight)
        : helper(helper),
          weight(weight) {}

    static auto Create(const VisualInertialAlignHelper<Order> &helper, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<VisualInertialAlignFactor>(
            new VisualInertialAlignFactor(helper, weight));
    }

    static std::size_t TypeHashCode() { return typeid(VisualInertialAlignFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3_CmToBr | POS_CmInBr | POS_BiInBr | S_VEL | E_VEL | GRAVITY | SCALE ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        Eigen::Map<const Sophus::SO3<T>> const SO3_CmToBr(sKnots[0]);
        Eigen::Map<const Eigen::Vector3<T>> const POS_CmInBr(sKnots[1]);
        Eigen::Map<const Eigen::Vector3<T>> const POS_BiInBr(sKnots[2]);
        Eigen::Map<const Eigen::Vector3<T>> const S_VEL(sKnots[3]);
        Eigen::Map<const Eigen::Vector3<T>> const E_VEL(sKnots[4]);
        Eigen::Map<const Eigen::Vector3<T>> const GRAVITY(sKnots[5]);
        T SCALE = sKnots[6][0];

        Eigen::Vector3<T> velLeft =
            helper.velVec.template cast<T>() - helper.velMat.template cast<T>() * POS_BiInBr;
        Eigen::Vector3<T> velRight =
            E_VEL - S_VEL - helper.diffCMat.template cast<T>() * POS_CmInBr - GRAVITY * helper.dt;

        Eigen::Vector3<T> posLeft =
            helper.posVec.template cast<T>() - helper.posMat.template cast<T>() * POS_BiInBr;
        Eigen::Vector3<T> posRight =
            helper.mtSo3BrToBr0.template cast<T>() * SO3_CmToBr *
                (SCALE * helper.diffDMat.template cast<T>()) -
            helper.diffEMat.template cast<T>() * POS_CmInBr -
            (S_VEL - helper.sCMat.template cast<T>() * POS_CmInBr) * helper.dt -
            0.5 * GRAVITY * helper.dt2;

        Eigen::Map<Eigen::Vector6<T>> residuals(sResiduals);

        residuals.template block<3, 1>(0, 0) = velLeft - velRight;
        residuals.template block<3, 1>(3, 0) = posLeft - posRight;

        residuals = T(weight) * residuals;

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

extern template struct VisualInertialAlignFactor<Configor::Prior::SplineOrder>;
}  // namespace ns_ikalibr

#endif  // IKALIBR_VISUAL_INERTIAL_ALIGN_FACTOR_HPP
