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

#ifndef IKALIBR_VISUAL_REPROJ_FACTOR_HPP
#define IKALIBR_VISUAL_REPROJ_FACTOR_HPP

#include <utility>
#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/ceres.h"
#include "veta/veta.h"
#include "util/utils.hpp"

_3_

namespace ns_ikalibr {
    struct VisualReProjCorr {
    public:
        using Ptr = std::shared_ptr<VisualReProjCorr>;

    public:
        double ti, tj;

        Eigen::Vector2d fi, fj;

        // row / image height
        double li, lj;

        double weight;

    public:
        VisualReProjCorr(double ti, double tj, Eigen::Vector2d fi, Eigen::Vector2d fj,
                         double li, double lj, double weight)
                : ti(ti), tj(tj), fi(std::move(fi)), fj(std::move(fj)), li(li), lj(lj), weight(weight) {}

        VisualReProjCorr() = default;
    };

    struct VisualReProjCorrSeq {
    public:
        using Ptr = std::shared_ptr<VisualReProjCorrSeq>;

    public:
        std::vector<VisualReProjCorr> corrs;

        std::unique_ptr<double> invDepthFir;

        ns_veta::IndexT lmId{};

        ns_veta::IndexT firObvViewId{};
        ns_veta::Observation firObv;

    public:
        VisualReProjCorrSeq() = default;
    };

    template<int Order, int TimeDeriv>
    struct VisualReProjFactor {
    private:
        ns_ctraj::SplineMeta<Order> _so3Meta, _scaleMeta;
        VisualReProjCorr _corr;

        double _so3DtInv, _scaleDtInv;
        double _weight;
    public:
        explicit VisualReProjFactor(ns_ctraj::SplineMeta<Order> rotMeta, ns_ctraj::SplineMeta<Order> linScaleMeta,
                                    VisualReProjCorr visualCorr, double weight)
                : _so3Meta(rotMeta), _scaleMeta(std::move(linScaleMeta)), _corr(std::move(visualCorr)),
                  _so3DtInv(1.0 / rotMeta.segments.front().dt),
                  _scaleDtInv(1.0 / _scaleMeta.segments.front().dt), _weight(weight) {}

        static auto
        Create(const ns_ctraj::SplineMeta<Order> &rotMeta, const ns_ctraj::SplineMeta<Order> &linScaleMeta,
               const VisualReProjCorr &visualCorr, double weight) {
            return new ceres::DynamicAutoDiffCostFunction<VisualReProjFactor>(
                    new VisualReProjFactor(rotMeta, linScaleMeta, visualCorr, weight)
            );
        }

        static std::size_t TypeHashCode() {
            return typeid(VisualReProjFactor).hash_code();
        }

        template<class T>
        void ComputeSE3BrToBr0(T const *const *sKnots, T *timeByBr,
                               Sophus::SO3<T> *SO3_BrToBr0, Eigen::Vector3<T> *POS_BrInBr0) const {
            std::pair<std::size_t, T> iuSo3, iuScale;
            _so3Meta.template ComputeSplineIndex(*timeByBr, iuSo3.first, iuSo3.second);
            _scaleMeta.template ComputeSplineIndex(*timeByBr, iuScale.first, iuScale.second);

            std::size_t SO3_OFFSET = iuSo3.first;
            std::size_t LIN_SCALE_OFFSET = iuScale.first + _so3Meta.NumParameters();

            ns_ctraj::CeresSplineHelperJet<T, Order>::template EvaluateLie(
                    sKnots + SO3_OFFSET, iuSo3.second, _so3DtInv, SO3_BrToBr0
            );

            ns_ctraj::CeresSplineHelperJet<T, Order>::template Evaluate<3, TimeDeriv>(
                    sKnots + LIN_SCALE_OFFSET, iuScale.second, _scaleDtInv, POS_BrInBr0
            );
        }

        template<class T>
        static void TransformImgToCam(const T *FX_INV, const T *FY_INV, const T *CX, const T *CY,
                                      const Eigen::Vector2<T> &feat, Eigen::Vector3<T> *P) {
            P->operator()(0) = (feat(0) - *CX) * *FX_INV;
            P->operator()(1) = (feat(1) - *CY) * *FY_INV;
            P->operator()(2) = (T) 1.0;
        }

        template<class T>
        static void TransformCamToImg(const T *FX, const T *FY, const T *CX, const T *CY,
                                      const Eigen::Vector3<T> &P, Eigen::Vector2<T> *feat) {
            feat->operator()(0) = *FX * P(0) + *CX;
            feat->operator()(1) = *FY * P(1) + *CY;
        }

    public:
        /**
         * param blocks:
         * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr | READOUT_TIME |
         *   FX | FY | CX | CY | GLOBAL_SCALE | INV_DEPTH ]
         */
        template<class T>
        bool operator()(T const *const *sKnots, T *sResiduals) const {

            std::size_t SO3_CmToBr_OFFSET = _so3Meta.NumParameters() + _scaleMeta.NumParameters();
            std::size_t POS_CmInBr_OFFSET = SO3_CmToBr_OFFSET + 1;
            std::size_t TO_CmToBr_OFFSET = POS_CmInBr_OFFSET + 1;
            std::size_t READOUT_TIME_OFFSET = TO_CmToBr_OFFSET + 1;
            std::size_t FX_OFFSET = READOUT_TIME_OFFSET + 1;
            std::size_t FY_OFFSET = FX_OFFSET + 1;
            std::size_t CX_OFFSET = FY_OFFSET + 1;
            std::size_t CY_OFFSET = CX_OFFSET + 1;
            std::size_t GLOBAL_SCALE_OFFSET = CY_OFFSET + 1;
            std::size_t INV_DEPTH_OFFSET = GLOBAL_SCALE_OFFSET + 1;

            // get value
            Eigen::Map<const Sophus::SO3<T>> SO3_CmToBr(sKnots[SO3_CmToBr_OFFSET]);
            Eigen::Map<const Eigen::Vector3<T>> POS_CmInBr(sKnots[POS_CmInBr_OFFSET]);
            Sophus::SE3<T> SE3_CmToBr(SO3_CmToBr, POS_CmInBr);

            T TO_CmToBr = sKnots[TO_CmToBr_OFFSET][0];
            T READOUT_TIME = sKnots[READOUT_TIME_OFFSET][0];

            T FX = sKnots[FX_OFFSET][0];
            T FX_INV = (T) 1.0 / FX;
            T FY = sKnots[FY_OFFSET][0];
            T FY_INV = (T) 1.0 / FY;
            T CX = sKnots[CX_OFFSET][0];
            T CY = sKnots[CY_OFFSET][0];

            T GLOBAL_SCALE = sKnots[GLOBAL_SCALE_OFFSET][0];
            T INV_DEPTH = sKnots[INV_DEPTH_OFFSET][0];
            T DEPTH = (T) 1.0 / INV_DEPTH;

            // calculate the so3 and lin scale offset for i-feat
            T timeIByBr = _corr.ti + TO_CmToBr + _corr.li * READOUT_TIME;
            Sophus::SE3<T> SE3_BrToBr0_I;
            ComputeSE3BrToBr0<T>(sKnots, &timeIByBr, &SE3_BrToBr0_I.so3(), &SE3_BrToBr0_I.translation());

            // calculate the so3 and lin scale offset for j-feat
            auto timeJByBr = _corr.tj + TO_CmToBr + _corr.lj * READOUT_TIME;
            Sophus::SE3<T> SE3_BrToBr0_J;
            ComputeSE3BrToBr0<T>(sKnots, &timeJByBr, &SE3_BrToBr0_J.so3(), &SE3_BrToBr0_J.translation());

            Sophus::SE3<T> SE3_BrIToBrJ = SE3_BrToBr0_J.inverse() * SE3_BrToBr0_I;
            Sophus::SE3<T> SE3_CmIToCmJ = SE3_CmToBr.inverse() * SE3_BrIToBrJ * SE3_CmToBr;

            Eigen::Vector3<T> PI;
            TransformImgToCam<T>(&FX_INV, &FY_INV, &CX, &CY, _corr.fi.cast<T>(), &PI);
            PI *= DEPTH * GLOBAL_SCALE;

            Eigen::Vector3<T> PJ = SE3_CmIToCmJ * PI;
            PJ /= PJ(2);
            Eigen::Vector2<T> fjPred;
            TransformCamToImg<T>(&FX, &FY, &CX, &CY, PJ, &fjPred);

            Eigen::Map<Eigen::Vector2<T>> residuals(sResiduals);
            residuals = fjPred - _corr.fj.cast<T>();
            residuals = T(_weight) * T(_corr.weight) * residuals;

            return true;
        }

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

}

#endif //IKALIBR_VISUAL_REPROJ_FACTOR_HPP
