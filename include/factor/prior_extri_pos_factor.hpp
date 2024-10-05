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

#ifndef IKALIBR_PRIOR_EXTRI_POS_FACTOR_HPP
#define IKALIBR_PRIOR_EXTRI_POS_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct PriorExtriPOSFactor {
private:
    const Eigen::Vector3d POS_Sen1InSen2;
    double _weight;

public:
    explicit PriorExtriPOSFactor(Eigen::Vector3d Sen1InSen2, double weight)
        : POS_Sen1InSen2(std::move(Sen1InSen2)),
          _weight(weight) {}

    static auto Create(const Eigen::Vector3d &Sen1InSen2, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<PriorExtriPOSFactor>(
            new PriorExtriPOSFactor(Sen1InSen2, weight));
    }

    static std::size_t TypeHashCode() { return typeid(PriorExtriPOSFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ POS_Sen1InRef | SO3_Sen2ToRef | POS_Sen2InRef ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        // t_ref_1 = R_ref_2 * t_2_1 + t_ref_2
        // t_2_1 = (R_ref_2)^(-1) * t_ref_1 - (R_ref_2)^(-1) * t_ref_2

        Eigen::Map<Eigen::Vector3<T> const> const POS_Sen1InRef(sKnots[0]);

        Eigen::Map<Sophus::SO3<T> const> const SO3_Sen2ToRef(sKnots[1]);
        Eigen::Map<Eigen::Vector3<T> const> const POS_Sen2InRef(sKnots[2]);

        Sophus::SO3<T> SO3_RefToSen2 = SO3_Sen2ToRef.inverse();
        Eigen::Vector3<T> POS_Sen1InSen2_Pred =
            SO3_RefToSen2 * POS_Sen1InRef - SO3_RefToSen2 * POS_Sen2InRef;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = T(_weight) * (POS_Sen1InSen2_Pred - POS_Sen1InSen2);

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_PRIOR_EXTRI_POS_FACTOR_HPP
