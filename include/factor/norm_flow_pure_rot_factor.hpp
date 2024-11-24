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

#ifndef NORM_FLOW_ROTATION_FACTOR_HPP
#define NORM_FLOW_ROTATION_FACTOR_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"
#include "factor/data_correspondence.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

template <int Order>
struct NormFlowPureRotFactor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta;
    NormFlow::Ptr _nf{};

    double _so3DtInv;
    double _weight;

public:
    explicit NormFlowPureRotFactor(ns_ctraj::SplineMeta<Order> so3Meta,
                                   NormFlow::Ptr nf,
                                   double weight)
        : _so3Meta(std::move(so3Meta)),
          _nf(std::move(nf)),
          _so3DtInv(1.0 / _so3Meta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &so3Meta,
                       const NormFlow::Ptr &nf,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<NormFlowPureRotFactor>(
            new NormFlowPureRotFactor(so3Meta, nf, weight));
    }

    static std::size_t TypeHashCode() { return typeid(NormFlowPureRotFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 | SO3_EsToBr | TO_EsToBr | FX | FY | CX | CY ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        // array offset
        std::size_t SO3_EsToBr_OFFSET = _so3Meta.NumParameters();
        std::size_t TO_EsToBr_OFFSET = SO3_EsToBr_OFFSET + 1;
        std::size_t FX_OFFSET = TO_EsToBr_OFFSET + 1;
        std::size_t FY_OFFSET = FX_OFFSET + 1;
        std::size_t CX_OFFSET = FY_OFFSET + 1;
        std::size_t CY_OFFSET = CX_OFFSET + 1;

        Eigen::Map<Sophus::SO3<T> const> const SO3_EsToBr(sKnots[SO3_EsToBr_OFFSET]);
        Sophus::SO3<T> SO3_BrToEs = SO3_EsToBr.inverse();

        T TO_EsToBr = sKnots[TO_EsToBr_OFFSET][0];
        T FX = sKnots[FX_OFFSET][0];
        T FY = sKnots[FY_OFFSET][0];
        T CX = sKnots[CX_OFFSET][0];
        T CY = sKnots[CY_OFFSET][0];

        T timeByBr = _nf->timestamp + TO_EsToBr;

        // calculate the so3 offset
        std::pair<std::size_t, T> iuCur;
        _so3Meta.ComputeSplineIndex(timeByBr, iuCur.first, iuCur.second);
        std::size_t SO3_OFFSET = iuCur.first;

        Sophus::SO3<T> SO3_BrToBr0;
        Sophus::SO3Tangent<T> ANG_VEL_BrToBr0InBr;
        ns_ctraj::CeresSplineHelperJet<T, Order>::EvaluateLie(
            sKnots + SO3_OFFSET, iuCur.second, _so3DtInv, &SO3_BrToBr0, &ANG_VEL_BrToBr0InBr);

        Eigen::Vector3<T> ANG_VEL_EsToBr0InEs = SO3_BrToEs * ANG_VEL_BrToBr0InBr;

        Eigen::Matrix<T, 2, 3> subBMat;
        OpticalFlowCorr::SubBMat<T>(&FX, &FY, &CX, &CY, _nf->p.cast<T>(), &subBMat);

        Eigen::Vector1<T> pred = _nf->nfDir.cast<T>().transpose() * (subBMat * ANG_VEL_EsToBr0InEs);

        Eigen::Map<Eigen::Vector1<T>> residuals(sResiduals);
        residuals.template block<1, 1>(0, 0) =
            T(_weight) * Eigen::Vector1<T>(_nf->nfNorm - pred(0));

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

extern template struct NormFlowPureRotFactor<Configor ::Prior::SplineOrder>;
}  // namespace ns_ikalibr

#endif  // NORM_FLOW_ROTATION_FACTOR_HPP
