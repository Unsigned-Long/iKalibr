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

#ifndef IKALIBR_POINT_TO_SURFEL_FACTOR_HPP
#define IKALIBR_POINT_TO_SURFEL_FACTOR_HPP

#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "util/utils.h"
#include "config/configor.h"
#include "factor/data_correspondence.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

template <int Order, int TimeDeriv>
struct PointToSurfelFactor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta, _scaleMeta;
    PointToSurfelCorr::Ptr _ptsCorr;

    double _so3DtInv, _scaleDtInv;
    double _weight;

public:
    explicit PointToSurfelFactor(const ns_ctraj::SplineMeta<Order> &so3Meta,
                                 const ns_ctraj::SplineMeta<Order> &scaleMeta,
                                 PointToSurfelCorr::Ptr ptsCorr,
                                 double weight)
        : _so3Meta(so3Meta),
          _scaleMeta(scaleMeta),
          _ptsCorr(std::move(ptsCorr)),
          _so3DtInv(1.0 / _so3Meta.segments.front().dt),
          _scaleDtInv(1.0 / _scaleMeta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &so3Meta,
                       const ns_ctraj::SplineMeta<Order> &scaleMeta,
                       const PointToSurfelCorr::Ptr &ptsCorr,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<PointToSurfelFactor>(
            new PointToSurfelFactor(so3Meta, scaleMeta, ptsCorr, weight));
    }

    static std::size_t TypeHashCode() { return typeid(PointToSurfelFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_LkToBr | POS_LkInBr | TO_LkToBr ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        std::size_t SO3_OFFSET;
        std::size_t LIN_SCALE_OFFSET;

        std::size_t SO3_LkToBr_OFFSET = _so3Meta.NumParameters() + _scaleMeta.NumParameters();
        std::size_t POS_LkInBr_OFFSET = SO3_LkToBr_OFFSET + 1;
        std::size_t TO_LkToBr_OFFSET = POS_LkInBr_OFFSET + 1;

        // get value
        Eigen::Map<const Sophus::SO3<T>> SO3_LkToBr(sKnots[SO3_LkToBr_OFFSET]);
        Eigen::Map<const Eigen::Vector3<T>> POS_LkInBr(sKnots[POS_LkInBr_OFFSET]);
        T TO_LkToBr = sKnots[TO_LkToBr_OFFSET][0];

        auto timeByBr = _ptsCorr->timestamp + TO_LkToBr;

        // calculate the so3 and lin scale offset
        std::pair<std::size_t, T> iuSo3, iuScale;
        _so3Meta.ComputeSplineIndex(timeByBr, iuSo3.first, iuSo3.second);
        _scaleMeta.ComputeSplineIndex(timeByBr, iuScale.first, iuScale.second);

        SO3_OFFSET = iuSo3.first;
        LIN_SCALE_OFFSET = iuScale.first + _so3Meta.NumParameters();

        Sophus::SO3<T> SO3_BrToBr0;
        ns_ctraj::CeresSplineHelperJet<T, Order>::EvaluateLie(sKnots + SO3_OFFSET, iuSo3.second,
                                                              _so3DtInv, &SO3_BrToBr0);

        Eigen::Vector3<T> POS_BrInBr0;
        ns_ctraj::CeresSplineHelperJet<T, Order>::template Evaluate<3, TimeDeriv>(
            sKnots + LIN_SCALE_OFFSET, iuScale.second, _scaleDtInv, &POS_BrInBr0);

        // construct the residuals
        Eigen::Vector3<T> pointInBr =
            SO3_LkToBr * _ptsCorr->pInScan.template cast<T>() + POS_LkInBr;
        Eigen::Vector3<T> pointInBr0 = SO3_BrToBr0 * pointInBr + POS_BrInBr0;

        Eigen::Vector3<T> planeNorm = _ptsCorr->surfelInW.head(3).template cast<T>();
        T distance = pointInBr0.dot(planeNorm) + T(_ptsCorr->surfelInW(3));

        Eigen::Map<Eigen::Matrix11<T>> residuals(sResiduals);
        residuals.template block<1, 1>(0, 0) = T(_weight) * Eigen::Matrix11<T>(distance);

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

extern template struct PointToSurfelFactor<Configor::Prior::SplineOrder, 2>;
extern template struct PointToSurfelFactor<Configor::Prior::SplineOrder, 1>;
extern template struct PointToSurfelFactor<Configor::Prior::SplineOrder, 0>;
}  // namespace ns_ikalibr
#endif  // IKALIBR_POINT_TO_SURFEL_FACTOR_HPP
