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

#ifndef IKALIBR_HAND_EYE_ROT_ALIGN_FACTOR_HPP
#define IKALIBR_HAND_EYE_ROT_ALIGN_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"
#include "config/configor.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
template <int Order>
struct HandEyeRotationAlignFactor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta;

    Sophus::SO3d SO3_CurLkToLastLk;
    double _tLastByLk, _tCurByLk;

    double _so3DtInv;
    double _weight;

public:
    explicit HandEyeRotationAlignFactor(ns_ctraj::SplineMeta<Order> so3Meta,
                                        double tLastByLk,
                                        double tCurByLk,
                                        const Sophus::SO3d &so3CurLkToLastLk,
                                        double weight)
        : _so3Meta(std::move(so3Meta)),
          SO3_CurLkToLastLk(so3CurLkToLastLk),
          _tLastByLk(tLastByLk),
          _tCurByLk(tCurByLk),
          _so3DtInv(1.0 / _so3Meta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &so3Meta,
                       double tLastByLk,
                       double tCurByLk,
                       const Sophus::SO3d &so3CurLkToLastLk,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<HandEyeRotationAlignFactor>(
            new HandEyeRotationAlignFactor(so3Meta, tLastByLk, tCurByLk, so3CurLkToLastLk, weight));
    }

    static std::size_t TypeHashCode() { return typeid(HandEyeRotationAlignFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 | SO3_LkToBr | TO_LkToBr ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        // array offset
        std::size_t CUR_SO3_OFFSET, LAST_SO3_OFFSET;
        std::size_t SO3_LkToBr_OFFSET = _so3Meta.NumParameters();
        std::size_t TO_LkToBr_OFFSET = SO3_LkToBr_OFFSET + 1;

        Eigen::Map<Sophus::SO3<T> const> const SO3_LkToBr(sKnots[SO3_LkToBr_OFFSET]);
        T TO_LkToBr = sKnots[TO_LkToBr_OFFSET][0];

        auto tLastByBr = _tLastByLk + TO_LkToBr;
        auto tCurByBr = _tCurByLk + TO_LkToBr;

        // calculate the so3 offset
        std::pair<std::size_t, T> iuLast;
        _so3Meta.ComputeSplineIndex(tLastByBr, iuLast.first, iuLast.second);
        LAST_SO3_OFFSET = iuLast.first;

        Sophus::SO3<T> SO3_LastBrToBr0;
        ns_ctraj::CeresSplineHelperJet<T, Order>::EvaluateLie(
            sKnots + LAST_SO3_OFFSET, iuLast.second, _so3DtInv, &SO3_LastBrToBr0);

        std::pair<std::size_t, T> iuCur;
        _so3Meta.ComputeSplineIndex(tCurByBr, iuCur.first, iuCur.second);
        CUR_SO3_OFFSET = iuCur.first;

        Sophus::SO3<T> SO3_CurBrToBr0;
        ns_ctraj::CeresSplineHelperJet<T, Order>::EvaluateLie(sKnots + CUR_SO3_OFFSET, iuCur.second,
                                                              _so3DtInv, &SO3_CurBrToBr0);

        Sophus::SO3<T> left = SO3_LkToBr * SO3_CurLkToLastLk;
        Sophus::SO3<T> right = (SO3_LastBrToBr0.inverse() * SO3_CurBrToBr0) * SO3_LkToBr;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = T(_weight) * (right.inverse() * left).log();

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

extern template struct HandEyeRotationAlignFactor<Configor::Prior::SplineOrder>;
}  // namespace ns_ikalibr

#endif  // IKALIBR_HAND_EYE_ROT_ALIGN_FACTOR_HPP
