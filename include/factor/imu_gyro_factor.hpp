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

#ifndef IKALIBR_IMU_GYRO_FACTOR_HPP
#define IKALIBR_IMU_GYRO_FACTOR_HPP

#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "sensor/imu.h"
#include "util/utils.h"
#include "config/configor.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
template <int Order>
struct IMUGyroFactor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta;
    IMUFrame::Ptr _frame{};

    double _so3DtInv;
    double _weight;

public:
    explicit IMUGyroFactor(ns_ctraj::SplineMeta<Order> so3Meta, IMUFrame::Ptr frame, double weight)
        : _so3Meta(std::move(so3Meta)),
          _frame(std::move(frame)),
          _so3DtInv(1.0 / _so3Meta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &so3Meta,
                       const IMUFrame::Ptr &frame,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<IMUGyroFactor>(
            new IMUGyroFactor(so3Meta, frame, weight));
    }

    static std::size_t TypeHashCode() { return typeid(IMUGyroFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 | GYRO_BIAS | GYRO_MAP_COEFF | SO3_AtoG | SO3_BiToBr | TO_BiToBr ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        // array offset
        std::size_t SO3_OFFSET;
        std::size_t GYRO_BIAS_OFFSET = _so3Meta.NumParameters();
        std::size_t GYRO_MAP_COEFF_OFFSET = GYRO_BIAS_OFFSET + 1;
        std::size_t SO3_AtoG_OFFSET = GYRO_MAP_COEFF_OFFSET + 1;
        std::size_t SO3_BiToBr_OFFSET = SO3_AtoG_OFFSET + 1;
        std::size_t TO_BiToBr_OFFSET = SO3_BiToBr_OFFSET + 1;

        T TO_BiToBr = sKnots[TO_BiToBr_OFFSET][0];

        auto timeByBr = _frame->GetTimestamp() + TO_BiToBr;

        // calculate the so3 offset
        std::pair<std::size_t, T> iuCur;
        _so3Meta.ComputeSplineIndex(timeByBr, iuCur.first, iuCur.second);
        SO3_OFFSET = iuCur.first;

        Sophus::SO3<T> SO3_BrToBr0;
        Sophus::SO3Tangent<T> SO3_VEL_BrToBr0InBr;
        ns_ctraj::CeresSplineHelperJet<T, Order>::EvaluateLie(
            sKnots + SO3_OFFSET, iuCur.second, _so3DtInv, &SO3_BrToBr0, &SO3_VEL_BrToBr0InBr);

        Eigen::Map<const Eigen::Vector3<T>> gyroBias(sKnots[GYRO_BIAS_OFFSET]);
        auto gyroCoeff = sKnots[GYRO_MAP_COEFF_OFFSET];
        Eigen::Matrix33<T> gyroMapMat = Eigen::Matrix33<T>::Zero();
        gyroMapMat.diagonal() = Eigen::Map<const Eigen::Vector3<T>>(gyroCoeff, 3);
        gyroMapMat(0, 1) = *(gyroCoeff + 3);
        gyroMapMat(0, 2) = *(gyroCoeff + 4);
        gyroMapMat(1, 2) = *(gyroCoeff + 5);

        Eigen::Map<Sophus::SO3<T> const> const SO3_AtoG(sKnots[SO3_AtoG_OFFSET]);
        Eigen::Map<Sophus::SO3<T> const> const SO3_BiToBr(sKnots[SO3_BiToBr_OFFSET]);
        Sophus::SO3<T> SO3_BiToBr0 = SO3_BrToBr0 * SO3_BiToBr;
        Sophus::SO3Tangent<T> SO3_VEL_BrToBr0InBr0 = SO3_BrToBr0 * SO3_VEL_BrToBr0InBr;

        Eigen::Vector3<T> pred =
            (gyroMapMat * (SO3_AtoG * SO3_BiToBr0.inverse() * SO3_VEL_BrToBr0InBr0)).eval() +
            gyroBias;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = pred - _frame->GetGyro().template cast<T>();
        residuals = T(_weight) * residuals;

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

extern template struct IMUGyroFactor<Configor ::Prior::SplineOrder>;
}  // namespace ns_ikalibr

#endif  // IKALIBR_IMU_GYRO_FACTOR_HPP
