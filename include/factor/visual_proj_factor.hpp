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

#ifndef IKALIBR_VISUAL_PROJ_FACTOR_HPP
#define IKALIBR_VISUAL_PROJ_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "veta/camera/pinhole_brown.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct VisualProjFactor {
private:
    Eigen::Vector2d _feat;
    ns_veta::PinholeIntrinsic::Ptr _intri;
    double fx, fy, cx, cy;
    double _weight;

public:
    explicit VisualProjFactor(Eigen::Vector2d feat,
                              const ns_veta::PinholeIntrinsic::Ptr &intri,
                              double weight)
        : _feat(std::move(feat)),
          fx(intri->FocalX()),
          fy(intri->FocalY()),
          cx(intri->PrincipalPoint()(0)),
          cy(intri->PrincipalPoint()(1)),
          _weight(weight) {}

    static auto Create(const Eigen::Vector2d &feat,
                       const ns_veta::PinholeIntrinsic::Ptr &intri,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<VisualProjFactor>(
            new VisualProjFactor(feat, intri, weight));
    }

    static std::size_t TypeHashCode() { return typeid(VisualProjFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3_CurCToW | POS_CurCInW | POS_LMInW ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        std::size_t SO3_CurCToW_OFFSET = 0;
        std::size_t POS_CurCInW_OFFSET = 1;
        std::size_t POS_LMInW_OFFSET = 2;

        // get value
        Eigen::Map<const Sophus::SO3<T>> SO3_CurCToW(sKnots[SO3_CurCToW_OFFSET]);
        Eigen::Map<const Eigen::Vector3<T>> POS_CurCInW(sKnots[POS_CurCInW_OFFSET]);
        Eigen::Map<const Eigen::Vector3<T>> POS_LMInW(sKnots[POS_LMInW_OFFSET]);

        Eigen::Vector3<T> POS_LMInCurC = SO3_CurCToW.inverse() * (POS_LMInW - POS_CurCInW);

        T invZ = 1.0 / POS_LMInCurC(2);

        Eigen::Map<Eigen::Vector2<T>> residuals(sResiduals);
        residuals(0) = fx * POS_LMInCurC(0) * invZ + cx - _feat(0);
        residuals(1) = fy * POS_LMInCurC(1) * invZ + cy - _feat(1);
        residuals = T(_weight) * residuals;

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_VISUAL_PROJ_FACTOR_HPP
