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

#ifndef PPP_TRIFOCAL_TENSOR_VEL_FACTOR_HPP
#define PPP_TRIFOCAL_TENSOR_VEL_FACTOR_HPP

#include "util/utils.h"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/core/spline_bundle.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "factor/data_correspondence.h"
#include "veta/camera/pinhole.h"
#include "factor/visual_optical_flow_factor.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
/**
 * attention: this trifocal tensor velocity factor use the second view as the reference
 */
template <int Order>
struct PPPTrifocalTensorVelFactorHelper {
public:
    using Ptr = std::shared_ptr<PPPTrifocalTensorVelFactorHelper>;
    using SplineBundleType = ns_ctraj::SplineBundle<Order>;
    using So3SplineType = typename SplineBundleType::So3SplineType;

    constexpr static int ONE = 0;
    constexpr static int TWO = 1;
    constexpr static int THREE = 2;

public:
    // bearing vectors
    std::array<Eigen::Vector3d, 3> bvs;

    // t1 - t2
    double dt12;
    // t3 - t2
    double dt32;

    Eigen::Matrix3d ROT_Cam2ToCam1;
    Eigen::Matrix3d ROT_Cam2ToCam3;

    std::array<Eigen::Vector3d, 3> ROT_VEC_Cam2ToCam1;
    std::array<Eigen::Vector3d, 3> ROT_VEC_Cam2ToCam3;

public:
    PPPTrifocalTensorVelFactorHelper(const So3SplineType &so3Spline,
                                     double TO_CamToBr,
                                     double RS_READOUT,
                                     Sophus::SO3d SO3_CamToBr,
                                     const OpticalFlowCorr::Ptr &corr,
                                     const ns_veta::PinholeIntrinsic::Ptr &intri) {
        std::array<Sophus::SO3d, 3> SO3_CamToW;

        for (int i = 0; i < 3; ++i) {
            const double x = corr->xTraceAry[i], y = corr->yTraceAry[i];
            bvs[i] = intri->ImgToCam({x, y}).homogeneous().normalized();

            const double t = corr->timeAry[i];
            double timeByBr = t + TO_CamToBr + RS_READOUT * corr->rdFactorAry[i];
            SO3_CamToW[i] = so3Spline.Evaluate(timeByBr) * SO3_CamToBr;
        }

        ROT_Cam2ToCam1 = (SO3_CamToW[ONE].inverse() * SO3_CamToW[TWO]).matrix();
        ROT_Cam2ToCam3 = (SO3_CamToW[THREE].inverse() * SO3_CamToW[TWO]).matrix();

        for (int i = 0; i < 3; ++i) {
            ROT_VEC_Cam2ToCam1[i] = ROT_Cam2ToCam1.col(i);
            ROT_VEC_Cam2ToCam3[i] = ROT_Cam2ToCam3.col(i);
        }

        dt12 = corr->timeAry[ONE] - corr->timeAry[TWO];
        dt32 = corr->timeAry[THREE] - corr->timeAry[TWO];
    }
};

template <int Order>
struct PPPTrifocalTensorVelFactor {
public:
    constexpr static int ONE = 0;
    constexpr static int TWO = 1;
    constexpr static int THREE = 2;

private:
    PPPTrifocalTensorVelFactorHelper<Order> _helper;

    double _weight;

public:
    explicit PPPTrifocalTensorVelFactor(const PPPTrifocalTensorVelFactorHelper<Order> &helper,
                                        double weight)
        : _helper(helper),
          _weight(weight) {}

    static auto Create(const PPPTrifocalTensorVelFactorHelper<Order> &helper, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<PPPTrifocalTensorVelFactor>(
            new PPPTrifocalTensorVelFactor(helper, weight));
    }

    static std::size_t TypeHashCode() { return typeid(PPPTrifocalTensorVelFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ LIN_VEL_DIR_CmToWInCm ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        Eigen::Map<const Eigen::Vector3<T>> LIN_VEL_DIR_CmToWInCm(sKnots[0]);

        Eigen::Vector3<T> bv1 = _helper.bvs[ONE].template cast<T>();
        Eigen::Vector3<T> bv2 = _helper.bvs[TWO].template cast<T>();
        Eigen::Vector3<T> bv3 = _helper.bvs[THREE].template cast<T>();

        Eigen::Matrix33<T> m1 = Eigen::Matrix33<T>::Zero();

        for (int i = 0; i < 3; ++i) {
            Eigen::Matrix33<T> v1 =
                _helper.ROT_VEC_Cam2ToCam1[i] *
                (_helper.ROT_Cam2ToCam3 * LIN_VEL_DIR_CmToWInCm * (-_helper.dt32)).transpose();
            Eigen::Matrix33<T> v2 =
                (_helper.ROT_Cam2ToCam1 * LIN_VEL_DIR_CmToWInCm * (-_helper.dt12)) *
                _helper.ROT_VEC_Cam2ToCam3[i].transpose();
            Eigen::Matrix33<T> Ti = v1 - v2;
            m1 += bv2(i) * Ti;
        }
        Eigen::Matrix33<T> m2 = Sophus::SO3<T>::hat(bv1) * m1 * Sophus::SO3<T>::hat(bv3);
        Eigen::Vector4<T> resVec;
        resVec(0) = m2(0, 0);
        resVec(1) = m2(0, 2);
        resVec(2) = m2(2, 0);
        resVec(3) = m2(2, 2);

        Eigen::Map<Eigen::Vector4<T>> residuals(sResiduals);
        residuals = T(_weight) * resVec;
        return true;
    }
};

template <int Order, int TimeDeriv, bool TruePosFalseVel>
struct PPPTrifocalTensorFactor {
private:
    ns_ctraj::SplineMeta<Order> _so3Meta, _scaleMeta;
    OpticalFlowCorr::Ptr _corr;

    double _so3DtInv, _scaleDtInv;
    double _weight;

public:
    explicit PPPTrifocalTensorFactor(ns_ctraj::SplineMeta<Order> rotMeta,
                                     ns_ctraj::SplineMeta<Order> linScaleMeta,
                                     OpticalFlowCorr::Ptr ofCorr,
                                     double weight)
        : _so3Meta(rotMeta),
          _scaleMeta(std::move(linScaleMeta)),
          _corr(std::move(ofCorr)),
          _so3DtInv(1.0 / rotMeta.segments.front().dt),
          _scaleDtInv(1.0 / _scaleMeta.segments.front().dt),
          _weight(weight) {}

    static auto Create(const ns_ctraj::SplineMeta<Order> &rotMeta,
                       const ns_ctraj::SplineMeta<Order> &linScaleMeta,
                       const OpticalFlowCorr::Ptr &ofCorr,
                       double weight) {
        return new ceres::DynamicAutoDiffCostFunction<PPPTrifocalTensorFactor>(
            new PPPTrifocalTensorFactor(rotMeta, linScaleMeta, ofCorr, weight));
    }

    static std::size_t TypeHashCode() { return typeid(PPPTrifocalTensorFactor).hash_code(); }

    template <class T>
    void ComputeSE3BrToBr0ByPosSpline(T const *const *sKnots,
                                      T *timeByBr,
                                      Sophus::SO3<T> *SO3_BrToBr0,
                                      Eigen::Vector3<T> *POS_BrInBr0) const {
        using ReprojHelper = VisualOpticalFlowReProjFactor<Order, TimeDeriv, true, TruePosFalseVel>;
        ReprojHelper::template ComputeSE3BrToBr0ByPosSpline(sKnots, timeByBr, SO3_BrToBr0,
                                                            POS_BrInBr0, _so3Meta, _so3DtInv,
                                                            _scaleMeta, _scaleDtInv);
    }

    template <class T>
    void ComputeSE3Br1ToBr2ByVelSpline(T const *const *sKnots,
                                       T *t1ByBr,
                                       T *t2ByBr,
                                       Sophus::SO3<T> *SO3_Br1ToBr2,
                                       Eigen::Vector3<T> *POS_Br1InBr2) const {
        using ReprojHelper = VisualOpticalFlowReProjFactor<Order, TimeDeriv, true, TruePosFalseVel>;
        ReprojHelper::template ComputeSE3Br1ToBr2ByVelSpline(sKnots, t1ByBr, t2ByBr, SO3_Br1ToBr2,
                                                             POS_Br1InBr2, _so3Meta, _so3DtInv,
                                                             _scaleMeta, _scaleDtInv);
    }

public:
    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
     * READOUT_TIME | FX | FY | CX | CY ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        std::size_t SO3_CmToBr_OFFSET = _so3Meta.NumParameters() + _scaleMeta.NumParameters();
        std::size_t POS_CmInBr_OFFSET = SO3_CmToBr_OFFSET + 1;
        std::size_t TO_CmToBr_OFFSET = POS_CmInBr_OFFSET + 1;
        std::size_t READOUT_TIME_OFFSET = TO_CmToBr_OFFSET + 1;
        std::size_t FX_OFFSET = READOUT_TIME_OFFSET + 1;
        std::size_t FY_OFFSET = FX_OFFSET + 1;
        std::size_t CX_OFFSET = FY_OFFSET + 1;
        std::size_t CY_OFFSET = CX_OFFSET + 1;

        // get value
        Eigen::Map<const Sophus::SO3<T>> SO3_CmToBr(sKnots[SO3_CmToBr_OFFSET]);
        Eigen::Map<const Eigen::Vector3<T>> POS_CmInBr(sKnots[POS_CmInBr_OFFSET]);
        Sophus::SE3<T> SE3_CmToBr(SO3_CmToBr, POS_CmInBr);

        T TO_CmToBr = sKnots[TO_CmToBr_OFFSET][0];
        T READOUT_TIME = sKnots[READOUT_TIME_OFFSET][0];

        T FX = sKnots[FX_OFFSET][0];
        T FX_INV = (T)1.0 / FX;
        T FY = sKnots[FY_OFFSET][0];
        T FY_INV = (T)1.0 / FY;
        T CX = sKnots[CX_OFFSET][0];
        T CY = sKnots[CY_OFFSET][0];

        Sophus::SE3<T> SE3_BrMidToBrFir, SE3_BrMidToBrLast;

        T timeByBrFir = _corr->timeAry.at(0) + TO_CmToBr + _corr->rdFactorAry.at(0) * READOUT_TIME;
        T timeByBrMid = _corr->timeAry.at(1) + TO_CmToBr + _corr->rdFactorAry.at(1) * READOUT_TIME;
        T timeByBrLast = _corr->timeAry.at(2) + TO_CmToBr + _corr->rdFactorAry.at(2) * READOUT_TIME;

        if constexpr (TruePosFalseVel) {
            // calculate the so3 and lin scale offset for first feat
            Sophus::SE3<T> SE3_BrToBr0_Fir;
            ComputeSE3BrToBr0ByPosSpline(sKnots, &timeByBrFir, &SE3_BrToBr0_Fir.so3(),
                                         &SE3_BrToBr0_Fir.translation());

            // calculate the so3 and lin scale offset for middle feat
            Sophus::SE3<T> SE3_BrToBr0_Mid;
            ComputeSE3BrToBr0ByPosSpline<T>(sKnots, &timeByBrMid, &SE3_BrToBr0_Mid.so3(),
                                            &SE3_BrToBr0_Mid.translation());

            // calculate the so3 and lin scale offset for last feat
            Sophus::SE3<T> SE3_BrToBr0_Last;
            ComputeSE3BrToBr0ByPosSpline<T>(sKnots, &timeByBrLast, &SE3_BrToBr0_Last.so3(),
                                            &SE3_BrToBr0_Last.translation());

            SE3_BrMidToBrFir = SE3_BrToBr0_Fir.inverse() * SE3_BrToBr0_Mid;
            SE3_BrMidToBrLast = SE3_BrToBr0_Last.inverse() * SE3_BrToBr0_Mid;
        } else {
            ComputeSE3Br1ToBr2ByVelSpline<T>(sKnots, &timeByBrMid, &timeByBrFir,
                                             &SE3_BrMidToBrFir.so3(),
                                             &SE3_BrMidToBrFir.translation());

            ComputeSE3Br1ToBr2ByVelSpline<T>(sKnots, &timeByBrMid, &timeByBrLast,
                                             &SE3_BrMidToBrLast.so3(),
                                             &SE3_BrMidToBrLast.translation());
        }

        Sophus::SE3<T> SE3_CmMidToCmFir = SE3_CmToBr.inverse() * SE3_BrMidToBrFir * SE3_CmToBr;
        Sophus::SE3<T> SE3_CmMidToCmLast = SE3_CmToBr.inverse() * SE3_BrMidToBrLast * SE3_CmToBr;

        Eigen::Matrix33<T> SO3_CmMidToCmFir = SE3_CmMidToCmFir.so3().matrix();
        Eigen::Matrix33<T> SO3_CmMidToCmLast = SE3_CmMidToCmLast.so3().matrix();

        Eigen::Vector3<T> POS_CmMidInCmFir = SE3_CmMidToCmFir.translation();
        Eigen::Vector3<T> POS_CmMidInCmLast = SE3_CmMidToCmLast.translation();

        std::array<Eigen::Vector3<T>, 3> bvs;

        for (int i = 0; i < 3; ++i) {
            const double x = _corr->xTraceAry[i], y = _corr->yTraceAry[i];
            Eigen::Vector2<T> f = Eigen::Vector2d(x, y).cast<T>();
            Eigen::Vector3<T> p;
            VisualReProjCorr::TransformImgToCam<T>(&FX_INV, &FY_INV, &CX, &CY, f, &p);
            bvs[i] = p.normalized();
        }

        Eigen::Matrix33<T> m1 = Eigen::Matrix33<T>::Zero();

        for (int i = 0; i < 3; ++i) {
            Eigen::Matrix33<T> v1 = SO3_CmMidToCmFir.col(i) * POS_CmMidInCmLast.transpose();
            Eigen::Matrix33<T> v2 = POS_CmMidInCmFir * SO3_CmMidToCmLast.col(i).transpose();
            Eigen::Matrix33<T> Ti = v1 - v2;
            m1 += bvs.at(1)(i) * Ti;
        }

        Eigen::Matrix33<T> m2 =
            Sophus::SO3<T>::hat(bvs.at(0)) * m1 * Sophus::SO3<T>::hat(bvs.at(2));

        Eigen::Vector4<T> resVec;
        resVec(0) = m2(0, 0);
        resVec(1) = m2(0, 2);
        resVec(2) = m2(2, 0);
        resVec(3) = m2(2, 2);

        Eigen::Map<Eigen::Vector4<T>> residuals(sResiduals);
        residuals = T(_weight) * resVec;

        return true;
    }
};
}  // namespace ns_ikalibr

#endif  // PPP_TRIFOCAL_TENSOR_VEL_FACTOR_HPP
