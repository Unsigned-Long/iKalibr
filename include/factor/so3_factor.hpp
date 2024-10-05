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

#ifndef IKALIBR_SO3_FACTOR_HPP
#define IKALIBR_SO3_FACTOR_HPP

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
template <int Order>
struct SO3Factor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta;

    double _time;
    Sophus::SO3d _so3;
    double _so3DtInv;
    double _weight;

public:
    explicit SO3Factor(ns_ctraj::SplineMeta<Order> so3Meta,
                       double time,
                       const Sophus::SO3d &so3,
                       double weight)
        : _so3Meta(std::move(so3Meta)),
          _time(time),
          _so3(so3),
          _so3DtInv(1.0 / _so3Meta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &so3Meta,
                       double time,
                       const Sophus::SO3d &so3,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<SO3Factor>(
            new SO3Factor(so3Meta, time, so3, weight));
    }

    static std::size_t TypeHashCode() { return typeid(SO3Factor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        // array offset
        std::size_t SO3_OFFSET;

        // calculate the so3 offset
        std::pair<std::size_t, double> iuCur;
        _so3Meta.ComputeSplineIndex(_time, iuCur.first, iuCur.second);
        SO3_OFFSET = iuCur.first;

        Sophus::SO3<T> SO3_BrToBr0;
        ns_ctraj::CeresSplineHelper<Order>::EvaluateLie(sKnots + SO3_OFFSET, iuCur.second,
                                                        _so3DtInv, &SO3_BrToBr0);

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = T(_weight) * (_so3.inverse().cast<T>() * SO3_BrToBr0).log();

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

extern template struct SO3Factor<Configor::Prior::SplineOrder>;
}  // namespace ns_ikalibr

#endif  // IKALIBR_SO3_FACTOR_HPP
