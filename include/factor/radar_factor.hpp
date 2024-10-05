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

#ifndef IKALIBR_RADAR_FACTOR_HPP
#define IKALIBR_RADAR_FACTOR_HPP

#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "sensor/radar.h"
#include "util/utils.h"
#include "config/configor.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
template <int Order, int TimeDeriv>
struct RadarFactor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta, _scaleMeta;

    RadarTarget::Ptr _frame;

    double _so3DtInv, _scaleDtInv;
    double _weight;

public:
    explicit RadarFactor(const ns_ctraj::SplineMeta<Order> &so3Meta,
                         const ns_ctraj::SplineMeta<Order> &scaleMeta,
                         RadarTarget::Ptr frame,
                         double weight)
        : _so3Meta(so3Meta),
          _scaleMeta(scaleMeta),
          _frame(std::move(frame)),
          _so3DtInv(1.0 / _so3Meta.segments.front().dt),
          _scaleDtInv(1.0 / _scaleMeta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &so3Meta,
                       const ns_ctraj::SplineMeta<Order> &scaleMeta,
                       const RadarTarget::Ptr &frame,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<RadarFactor>(
            new RadarFactor(so3Meta, scaleMeta, frame, weight));
    }

    static std::size_t TypeHashCode() { return typeid(RadarFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_RjToBr | POS_RjInBr | TO_RjToBr ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        std::size_t SO3_OFFSET, LIN_SCALE_OFFSET;
        std::size_t SO3_RjToBr_OFFSET = _so3Meta.NumParameters() + _scaleMeta.NumParameters();
        std::size_t POS_RjInBr_OFFSET = SO3_RjToBr_OFFSET + 1;
        std::size_t TO_RjToBr_OFFSET = POS_RjInBr_OFFSET + 1;

        // get value
        Eigen::Map<const Sophus::SO3<T>> SO3_RjToBr(sKnots[SO3_RjToBr_OFFSET]);
        Eigen::Map<const Eigen::Vector3<T>> POS_RjInBr(sKnots[POS_RjInBr_OFFSET]);
        T TO_RjToBr = sKnots[TO_RjToBr_OFFSET][0];

        auto timeByBr = _frame->GetTimestamp() + TO_RjToBr;

        // calculate the so3 and pos offset
        std::pair<std::size_t, T> iuSo3, iuScale;
        _so3Meta.ComputeSplineIndex(timeByBr, iuSo3.first, iuSo3.second);
        _scaleMeta.ComputeSplineIndex(timeByBr, iuScale.first, iuScale.second);

        SO3_OFFSET = iuSo3.first;
        LIN_SCALE_OFFSET = iuScale.first + _so3Meta.NumParameters();

        // query
        Sophus::SO3<T> SO3_BrToBr0;
        Eigen::Vector3<T> ANG_VEL_BrToBr0InBr;
        ns_ctraj::CeresSplineHelperJet<T, Order>::EvaluateLie(
            sKnots + SO3_OFFSET, iuSo3.second, _so3DtInv, &SO3_BrToBr0, &ANG_VEL_BrToBr0InBr);
        Eigen::Vector3<T> ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * ANG_VEL_BrToBr0InBr;

        Eigen::Vector3<T> LIN_VEL_BrInBr0;
        ns_ctraj::CeresSplineHelperJet<T, Order>::template Evaluate<3, TimeDeriv>(
            sKnots + LIN_SCALE_OFFSET, iuScale.second, _scaleDtInv, &LIN_VEL_BrInBr0);

        Eigen::Vector3<T> tarInRj = _frame->GetTargetXYZ().cast<T>();
        Eigen::Vector1<T> v1 =
            -tarInRj.transpose() * SO3_RjToBr.matrix().transpose() *
            SO3_BrToBr0.matrix().transpose() *
            (-Sophus::SO3<T>::hat(SO3_BrToBr0 * POS_RjInBr) * ANG_VEL_BrToBr0InBr0 +
             LIN_VEL_BrInBr0);

        T v2 = static_cast<T>(_frame->GetRadialVelocity());

        Eigen::Map<Eigen::Vector1<T>> residuals(sResiduals);
        residuals.template block<1, 1>(0, 0) =
            T(_weight) * Eigen::Vector1<T>(_frame->GetInvRange() * v1(0) - v2);

        return true;
    }
};

extern template struct RadarFactor<Configor::Prior::SplineOrder, 2>;
extern template struct RadarFactor<Configor::Prior::SplineOrder, 1>;
extern template struct RadarFactor<Configor::Prior::SplineOrder, 0>;
}  // namespace ns_ikalibr

#endif  // IKALIBR_RADAR_FACTOR_HPP
