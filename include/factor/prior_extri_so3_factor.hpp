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

#ifndef IKALIBR_PRIOR_EXTRI_SO3_FACTOR_HPP
#define IKALIBR_PRIOR_EXTRI_SO3_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct PriorExtriSO3Factor {
private:
    const Sophus::SO3d SO3_Sen1ToSen2;
    double _weight;

public:
    explicit PriorExtriSO3Factor(const Sophus::SO3d &Sen1ToSen2, double weight)
        : SO3_Sen1ToSen2(Sen1ToSen2),
          _weight(weight) {}

    static auto Create(const Sophus::SO3d &Sen1ToSen2, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<PriorExtriSO3Factor>(
            new PriorExtriSO3Factor(Sen1ToSen2, weight));
    }

    static std::size_t TypeHashCode() { return typeid(PriorExtriSO3Factor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3_Sen1ToRef | SO3_Sen2ToRef ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        Eigen::Map<Sophus::SO3<T> const> const SO3_Sen1ToRef(sKnots[0]);
        Eigen::Map<Sophus::SO3<T> const> const SO3_Sen2ToRef(sKnots[1]);

        Sophus::SO3<T> SO3_Sen1ToSen2_Pred = SO3_Sen2ToRef.inverse() * SO3_Sen1ToRef;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = T(_weight) * (SO3_Sen1ToSen2.inverse() * SO3_Sen1ToSen2_Pred).log();

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_PRIOR_EXTRI_SO3_FACTOR_HPP
