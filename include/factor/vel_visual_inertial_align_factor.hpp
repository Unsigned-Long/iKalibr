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

#ifndef VEL_VISUAL_INERTIAL_ALIGN_FACTOR_HPP
#define VEL_VISUAL_INERTIAL_ALIGN_FACTOR_HPP

#include "factor/rgbd_inertial_align_factor.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
template <int Order>
using VelVisualInertialAlignHelper = RGBDInertialAlignHelper<Order>;

template <int Order>
struct VelVisualInertialAlignFactor {
private:
    VelVisualInertialAlignHelper<Order> _helper;

    double _weight;

public:
    explicit VelVisualInertialAlignFactor(const VelVisualInertialAlignHelper<Order> &helper,
                                          double weight)
        : _helper(helper),
          _weight(weight) {}

    static auto Create(const VelVisualInertialAlignHelper<Order> &helper, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<VelVisualInertialAlignFactor>(
            new VelVisualInertialAlignFactor(helper, weight));
    }

    static std::size_t TypeHashCode() { return typeid(VelVisualInertialAlignFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ POS_BiInBr | SO3_CmToBr | POS_CmInBr | GRAVITY | VEL_SCALE_S | VEL_SCALE_E ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        Eigen::Map<const Eigen::Vector3<T>> POS_BiInBr(sKnots[0]);
        Eigen::Map<Sophus::SO3<T> const> const SO3_CmToBr(sKnots[1]);
        Eigen::Map<const Eigen::Vector3<T>> POS_CmInBr(sKnots[2]);
        Eigen::Map<const Eigen::Vector3<T>> GRAVITY(sKnots[3]);
        T sVelScale = sKnots[4][0];
        T eVelScale = sKnots[5][0];

        Eigen::Vector3<T> sLIN_VEL_BrToBr0InBr0 =
            _helper.sSO3_BrToBr0 * SO3_CmToBr * (sVelScale * _helper.sDVec) +
            Sophus::SO3<T>::hat(_helper.sSO3_BrToBr0 * POS_CmInBr) * _helper.sANG_VEL_BcToBc0;

        Eigen::Vector3<T> eLIN_VEL_BrToBr0InBr0 =
            _helper.eSO3_BrToBr0 * SO3_CmToBr * (eVelScale * _helper.eDVec) +
            Sophus::SO3<T>::hat(_helper.eSO3_BrToBr0 * POS_CmInBr) * _helper.eANG_VEL_BcToBc0;

        Eigen::Vector3<T> pred =
            eLIN_VEL_BrToBr0InBr0 - sLIN_VEL_BrToBr0InBr0 - GRAVITY * _helper.dt;
        Eigen::Vector3<T> mes =
            _helper.velVec.template cast<T>() - _helper.velMat.template cast<T>() * POS_BiInBr;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals.template block<3, 1>(0, 0) = T(_weight) * (pred - mes);

        return true;
    }
};

extern template struct VelVisualInertialAlignFactor<Configor::Prior::SplineOrder>;
}  // namespace ns_ikalibr

#endif  // VEL_VISUAL_INERTIAL_ALIGN_FACTOR_HPP
