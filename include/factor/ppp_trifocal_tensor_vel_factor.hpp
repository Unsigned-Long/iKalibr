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

        Eigen::Matrix3<T> m1 = Eigen::Matrix3<T>::Zero();

        for (int i = 0; i < 3; ++i) {
            Eigen::Matrix3<T> v1 =
                _helper.ROT_VEC_Cam2ToCam1[i] *
                (_helper.ROT_Cam2ToCam3 * LIN_VEL_DIR_CmToWInCm * (-_helper.dt32)).transpose();
            Eigen::Matrix3<T> v2 =
                (_helper.ROT_Cam2ToCam1 * LIN_VEL_DIR_CmToWInCm * (-_helper.dt12)) *
                _helper.ROT_VEC_Cam2ToCam3[i].transpose();
            Eigen::Matrix3<T> Ti = v1 - v2;
            m1 += bv2(i) * Ti;
        }
        Eigen::Matrix3<T> m2 = Sophus::SO3<T>::hat(bv1) * m1 * Sophus::SO3<T>::hat(bv3);

        // Eigen::Map<Eigen::Vector9<T>> residuals(sResiduals);
        // residuals = Eigen::Map<Eigen::Vector9<T>>(m2.data(), 9);

        Eigen::Map<Eigen::Vector1<T>> residuals(sResiduals);
        residuals.template block<1, 1>(0, 0) =
            T(_weight) * Eigen::Vector1<T>(ceres::sqrt((m2.array() * m2.array()).sum()));

        // Eigen::Map<Eigen::Vector1<T>> residuals(sResiduals);
        // T m2 = bv3.transpose() * m1 * bv1;
        // residuals.template block<1, 1>(0, 0) = T(_weight) * Eigen::Vector1<T>(m2);
        return true;
    }
};
}  // namespace ns_ikalibr

#endif  // PPP_TRIFOCAL_TENSOR_VEL_FACTOR_HPP
