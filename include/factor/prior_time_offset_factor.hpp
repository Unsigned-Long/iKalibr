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

#ifndef IKALIBR_PRIOR_TIME_OFFSET_FACTOR_HPP
#define IKALIBR_PRIOR_TIME_OFFSET_FACTOR_HPP

#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct PriorTimeOffsetFactor {
private:
    const double TO_Sen1ToSen2;
    double _weight;

public:
    explicit PriorTimeOffsetFactor(const double &Sen1ToSen2, double weight)
        : TO_Sen1ToSen2(Sen1ToSen2),
          _weight(weight) {}

    static auto Create(const double &Sen1ToSen2, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<PriorTimeOffsetFactor>(
            new PriorTimeOffsetFactor(Sen1ToSen2, weight));
    }

    static std::size_t TypeHashCode() { return typeid(PriorTimeOffsetFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ TO_Sen1ToRef | TO_Sen2ToRef ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        T const TO_Sen1ToRef = sKnots[0][0];
        T const TO_Sen2ToRef = sKnots[1][0];

        T TO_Sen1ToSen2_Pred = TO_Sen1ToRef - TO_Sen2ToRef;

        Eigen::Map<Eigen::Vector1<T>> residuals(sResiduals);
        residuals(0) = T(_weight) * (TO_Sen1ToSen2 - TO_Sen1ToSen2_Pred);

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_PRIOR_TIME_OFFSET_FACTOR_HPP
