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

#ifndef IKALIBR_LIN_SCALE_FACTOR_HPP
#define IKALIBR_LIN_SCALE_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"
#include "config/configor.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
template <int Order, int TimeDeriv>
struct LinearScaleDerivFactor {
private:
    ns_ctraj::SplineMeta<Order> _scaleMeta;

    double _time;
    Eigen::Vector3d _linScaleOfDeriv;
    double _scaleDtInv;
    double _weight;

public:
    explicit LinearScaleDerivFactor(ns_ctraj::SplineMeta<Order> linAcceMeta,
                                    double time,
                                    Eigen::Vector3d linScaleOfDeriv,
                                    double weight)
        : _scaleMeta(std::move(linAcceMeta)),
          _time(time),
          _linScaleOfDeriv(std::move(linScaleOfDeriv)),
          _scaleDtInv(1.0 / _scaleMeta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &splineMeta,
                       double time,
                       const Eigen::Vector3d &linScaleOfDeriv,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<LinearScaleDerivFactor>(
            new LinearScaleDerivFactor(splineMeta, time, linScaleOfDeriv, weight));
    }

    static std::size_t TypeHashCode() { return typeid(LinearScaleDerivFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ LIN_SCALE | ... | LIN_SCALE ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        std::size_t LIN_SCALE_OFFSET;

        // calculate the so3 and pos offset
        std::pair<std::size_t, double> iuScale;
        _scaleMeta.ComputeSplineIndex(_time, iuScale.first, iuScale.second);
        LIN_SCALE_OFFSET = iuScale.first;

        Eigen::Vector3<T> linScaleOfDeriv;
        ns_ctraj::CeresSplineHelper<Order>::template Evaluate<T, 3, TimeDeriv>(
            sKnots + LIN_SCALE_OFFSET, iuScale.second, _scaleDtInv, &linScaleOfDeriv);

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = T(_weight) * (linScaleOfDeriv - _linScaleOfDeriv);

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

extern template struct LinearScaleDerivFactor<Configor::Prior::SplineOrder, 2>;
extern template struct LinearScaleDerivFactor<Configor::Prior::SplineOrder, 1>;
extern template struct LinearScaleDerivFactor<Configor::Prior::SplineOrder, 0>;
}  // namespace ns_ikalibr

#endif  // IKALIBR_LIN_SCALE_FACTOR_HPP
