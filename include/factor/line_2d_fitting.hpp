// iKalibr: Unified Targetless Spatiotemporal Calibration Framework
// Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China
// https://github.com/Unsigned-Long/iKalibr.git
// Author: Shuolong Chen (shlchen@whu.edu.cn)
// GitHub: https://github.com/Unsigned-Long
//  ORCID: 0000-0002-5283-9057
// Purpose: See .h/.hpp file.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * The names of its contributors can not be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
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

#ifndef LINE_2D_FITTING_HPP
#define LINE_2D_FITTING_HPP

#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct Line2DFittingFactor {
private:
    Eigen::Vector2d p;
    Eigen::Vector2d dir;
    double lambda;
    double weight;

public:
    explicit Line2DFittingFactor(const Eigen::Vector2d &p,
                                 const Eigen::Vector2d &dir,
                                 double lambda,
                                 double weight)
        : p(p),
          dir(dir),
          lambda(lambda),
          weight(weight) {}

    static auto Create(const Eigen::Vector2d &p,
                       const Eigen::Vector2d &dir,
                       double lambda,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<Line2DFittingFactor>(
            new Line2DFittingFactor(p, dir, lambda, weight));
    }

    static std::size_t TypeHashCode() { return typeid(Line2DFittingFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ COS, SIN | RHO ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        T cos = sKnots[0][0];
        T sin = sKnots[0][1];
        T rho = sKnots[1][0];

        T e1 = p(0) * cos + p(1) * sin - rho;
        T e2 = lambda * (1.0 - dir(0) * cos - dir(1) * sin);

        Eigen::Map<Eigen::Vector2<T>> residuals(sResiduals);
        residuals.template block<2, 1>(0, 0) = T(weight) * Eigen::Vector2<T>(e1, e2);

        return true;
    }
};
}  // namespace ns_ikalibr

#endif  // LINE_2D_FITTING_HPP
