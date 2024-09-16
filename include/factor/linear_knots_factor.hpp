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

#ifndef IKALIBR_LINEAR_KNOTS_FACTOR_HPP
#define IKALIBR_LINEAR_KNOTS_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct RdLinearKnotsFactor {
private:
    double weight;

public:
    explicit RdLinearKnotsFactor(double weight)
        : weight(weight) {}

    static auto Create(double weight) {
        return new ceres::DynamicAutoDiffCostFunction<RdLinearKnotsFactor>(
            new RdLinearKnotsFactor(weight));
    }

    static std::size_t TypeHashCode() { return typeid(RdLinearKnotsFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ VEL | VEL | VEL ]
     */
    template <class T>
    bool operator()(T const *const *parBlocks, T *sResiduals) const {
        Eigen::Map<const Eigen::Vector3<T>> firKnot(parBlocks[0]);
        Eigen::Map<const Eigen::Vector3<T>> sedKnot(parBlocks[1]);
        Eigen::Map<const Eigen::Vector3<T>> thdKnot(parBlocks[2]);

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = weight * (firKnot - 2.0 * sedKnot + thdKnot);

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct So3LinearKnotsFactor {
private:
    double weight;

public:
    explicit So3LinearKnotsFactor(double weight)
        : weight(weight) {}

    static auto Create(double weight) {
        return new ceres::DynamicAutoDiffCostFunction<So3LinearKnotsFactor>(
            new So3LinearKnotsFactor(weight));
    }

    static std::size_t TypeHashCode() { return typeid(So3LinearKnotsFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | SO3 | SO3 ]
     */
    template <class T>
    bool operator()(T const *const *parBlocks, T *sResiduals) const {
        Eigen::Map<Sophus::SO3<T> const> const SO3_firToW(parBlocks[0]);
        Eigen::Map<Sophus::SO3<T> const> const SO3_sedToW(parBlocks[1]);
        Eigen::Map<Sophus::SO3<T> const> const SO3_thdToW(parBlocks[2]);

        Sophus::SO3<T> SO3_firToSed = SO3_sedToW.inverse() * SO3_firToW;
        Sophus::SO3<T> SO3_SedToThd = SO3_thdToW.inverse() * SO3_sedToW;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = weight * (SO3_SedToThd.inverse() * SO3_firToSed).log();

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_LINEAR_KNOTS_FACTOR_HPP
