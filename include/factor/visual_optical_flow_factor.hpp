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

#ifndef VISUAL_OPTICAL_FLOW_FACTOR_HPP
#define VISUAL_OPTICAL_FLOW_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"
#include "config/configor.h"
#include "factor/data_correspondence.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

template <int Order, int TimeDeriv, bool IsInvDepth = true>
struct VisualOpticalFlowFactor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta, _scaleMeta;

    OpticalFlowCorr::Ptr _corr;

    double _so3DtInv, _scaleDtInv;
    double _weight;

public:
    explicit VisualOpticalFlowFactor(const ns_ctraj::SplineMeta<Order> &so3Meta,
                                     const ns_ctraj::SplineMeta<Order> &scaleMeta,
                                     OpticalFlowCorr::Ptr corr,
                                     double weight)
        : _so3Meta(so3Meta),
          _scaleMeta(scaleMeta),
          _corr(std::move(corr)),
          _so3DtInv(1.0 / _so3Meta.segments.front().dt),
          _scaleDtInv(1.0 / _scaleMeta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &so3Meta,
                       const ns_ctraj::SplineMeta<Order> &scaleMeta,
                       const OpticalFlowCorr::Ptr &corr,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<VisualOpticalFlowFactor>(
            new VisualOpticalFlowFactor(so3Meta, scaleMeta, corr, weight));
    }

    static std::size_t TypeHashCode() { return typeid(VisualOpticalFlowFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
     *   READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        std::size_t SO3_OFFSET, LIN_SCALE_OFFSET;
        std::size_t SO3_CmToBr_OFFSET = _so3Meta.NumParameters() + _scaleMeta.NumParameters();
        std::size_t POS_CmInBr_OFFSET = SO3_CmToBr_OFFSET + 1;
        std::size_t TO_CmToBr_OFFSET = POS_CmInBr_OFFSET + 1;
        std::size_t READOUT_TIME_OFFSET = TO_CmToBr_OFFSET + 1;
        std::size_t FX_OFFSET = READOUT_TIME_OFFSET + 1;
        std::size_t FY_OFFSET = FX_OFFSET + 1;
        std::size_t CX_OFFSET = FY_OFFSET + 1;
        std::size_t CY_OFFSET = CX_OFFSET + 1;
        std::size_t DEPTH_INFO_OFFSET = CY_OFFSET + 1;

        // get value
        Eigen::Map<const Sophus::SO3<T>> SO3_CmToBr(sKnots[SO3_CmToBr_OFFSET]);
        Eigen::Map<const Eigen::Vector3<T>> POS_CmInBr(sKnots[POS_CmInBr_OFFSET]);
        Sophus::SO3<T> SO3_BrToCm = SO3_CmToBr.inverse();

        T TO_CmToBr = sKnots[TO_CmToBr_OFFSET][0];
        T READOUT_TIME = sKnots[READOUT_TIME_OFFSET][0];

        T FX = sKnots[FX_OFFSET][0];
        T FY = sKnots[FY_OFFSET][0];
        T CX = sKnots[CX_OFFSET][0];
        T CY = sKnots[CY_OFFSET][0];

        T DEPTH_INFO = sKnots[DEPTH_INFO_OFFSET][0];

        auto timeByBr = _corr->MidPointTime(READOUT_TIME) + TO_CmToBr;

        // calculate the so3 and pos offset
        std::pair<std::size_t, T> iuSo3, iuScale;
        _so3Meta.template ComputeSplineIndex(timeByBr, iuSo3.first, iuSo3.second);
        _scaleMeta.template ComputeSplineIndex(timeByBr, iuScale.first, iuScale.second);

        SO3_OFFSET = iuSo3.first;
        LIN_SCALE_OFFSET = iuScale.first + _so3Meta.NumParameters();

        // query
        Sophus::SO3<T> SO3_BrToBr0;
        Eigen::Vector3<T> ANG_VEL_BrToBr0InBr;
        ns_ctraj::CeresSplineHelperJet<T, Order>::template EvaluateLie(
            sKnots + SO3_OFFSET, iuSo3.second, _so3DtInv, &SO3_BrToBr0, &ANG_VEL_BrToBr0InBr);

        Eigen::Vector3<T> ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * ANG_VEL_BrToBr0InBr;
        Eigen::Vector3<T> ANG_VEL_CmToBr0InCm = SO3_BrToCm * ANG_VEL_BrToBr0InBr;

        Eigen::Vector3<T> LIN_VEL_BrToBr0InBr0;
        ns_ctraj::CeresSplineHelperJet<T, Order>::template Evaluate<3, TimeDeriv>(
            sKnots + LIN_SCALE_OFFSET, iuScale.second, _scaleDtInv, &LIN_VEL_BrToBr0InBr0);

        Eigen::Vector3<T> LIN_VEL_CmToBr0InBr0 =
            -Sophus::SO3<T>::hat(SO3_BrToBr0 * POS_CmInBr) * ANG_VEL_BrToBr0InBr0 +
            LIN_VEL_BrToBr0InBr0;

        Eigen::Vector3<T> LIN_VEL_CmToBr0InCm =
            SO3_BrToCm * SO3_BrToBr0.inverse() * LIN_VEL_CmToBr0InBr0;

        Eigen::Matrix<T, 2, 3> subAMat, subBMat;
        OpticalFlowCorr::SubMats<T>(&FX, &FY, &CX, &CY, _corr->MidPoint().cast<T>(), &subAMat,
                                    &subBMat);

        Eigen::Vector2<T> pred;
        if constexpr (IsInvDepth) {
            // inverse depth
            pred = DEPTH_INFO * subAMat * LIN_VEL_CmToBr0InCm + subBMat * ANG_VEL_CmToBr0InCm;
        } else {
            // depth
            pred =
                (1.0 / DEPTH_INFO) * subAMat * LIN_VEL_CmToBr0InCm + subBMat * ANG_VEL_CmToBr0InCm;
        }

        Eigen::Map<Eigen::Vector2<T>> residuals(sResiduals);
        residuals = T(_weight) * (pred - _corr->template MidPointVel(READOUT_TIME));

        return true;
    }
};

extern template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 2, true>;
extern template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 2, false>;
extern template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 1, true>;
extern template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 1, false>;
extern template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 0, true>;
extern template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 0, false>;
}  // namespace ns_ikalibr

#endif  // VISUAL_OPTICAL_FLOW_FACTOR_HPP
