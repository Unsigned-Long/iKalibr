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

#ifndef IKALIBR_INERTIAL_ALIGN_FACTOR_HPP
#define IKALIBR_INERTIAL_ALIGN_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct InertialAlignHelper {
public:
    using Ptr = std::shared_ptr<InertialAlignHelper>;

public:
    double dt;

    Eigen::Vector3d velVec;
    Eigen::Matrix3d velMat;

public:
    InertialAlignHelper(double timeDist,
                        const std::pair<Eigen::Vector3d, Eigen::Matrix3d> &velVecMat) {
        dt = timeDist;

        velVec = velVecMat.first;
        velMat = velVecMat.second;
    }
};

struct InertialAlignFactor {
private:
    InertialAlignHelper helper;
    double weight;

public:
    InertialAlignFactor(InertialAlignHelper helper, double weight)
        : helper(std::move(helper)),
          weight(weight) {}

    static auto Create(const InertialAlignHelper &helper, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<InertialAlignFactor>(
            new InertialAlignFactor(helper, weight));
    }

    static std::size_t TypeHashCode() { return typeid(InertialAlignFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ POS_BiInBr | START_VEL | END_VEL | GRAVITY ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        Eigen::Map<const Eigen::Vector3<T>> POS_BiInBr(sKnots[0]);
        Eigen::Map<const Eigen::Vector3<T>> START_VEL(sKnots[1]);
        Eigen::Map<const Eigen::Vector3<T>> END_VEL(sKnots[2]);
        Eigen::Map<const Eigen::Vector3<T>> GRAVITY(sKnots[3]);

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = (helper.velVec.cast<T>() - helper.velMat.cast<T>() * POS_BiInBr) -
                    (END_VEL - START_VEL - GRAVITY * helper.dt);
        residuals = T(weight) * residuals;

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_INERTIAL_ALIGN_FACTOR_HPP
