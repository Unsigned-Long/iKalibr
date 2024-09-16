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

#ifndef IKALIBR_RADAR_INERTIAL_ALIGN_FACTOR_HPP
#define IKALIBR_RADAR_INERTIAL_ALIGN_FACTOR_HPP

#include "util/utils.h"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/core/spline_bundle.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "sensor/radar.h"
#include "core/radar_velocity_sac.h"
#include "opengv/sac/Ransac.hpp"
#include "spdlog/spdlog.h"
#include "config/configor.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

template <int Order>
struct RadarInertialAlignHelper {
public:
    using Ptr = std::shared_ptr<RadarInertialAlignHelper>;
    using SplineBundleType = ns_ctraj::SplineBundle<Order>;
    using So3SplineType = typename SplineBundleType::So3SplineType;

public:
    double dt;
    Eigen::Vector3d velVec;
    Eigen::Matrix3d velMat;

    Eigen::Vector3d sDVec;
    Eigen::Vector3d eDVec;

    Eigen::Vector3d sANG_VEL_BcToBc0;
    Eigen::Vector3d eANG_VEL_BcToBc0;

    Sophus::SO3d sSO3_BrToBr0;
    Sophus::SO3d eSO3_BrToBr0;

public:
    RadarInertialAlignHelper(const So3SplineType &so3Spline,
                             const RadarTargetArray::Ptr &sRadarAry,
                             const RadarTargetArray::Ptr &eRadarAry,
                             double TO_RjToBr,
                             const std::pair<Eigen::Vector3d, Eigen::Matrix3d> &velVecMat,
                             double ransacThd) {
        double st = sRadarAry->GetTimestamp(), et = eRadarAry->GetTimestamp();
        dt = et - st;

        velVec = velVecMat.first;
        velMat = velVecMat.second;

        sDVec = RadarVelocityFromStaticTargetArray(sRadarAry, ransacThd);
        eDVec = RadarVelocityFromStaticTargetArray(eRadarAry, ransacThd);

        sANG_VEL_BcToBc0 = AngularVelBrToBr0(so3Spline, st + TO_RjToBr);
        eANG_VEL_BcToBc0 = AngularVelBrToBr0(so3Spline, et + TO_RjToBr);

        sSO3_BrToBr0 = so3Spline.Evaluate(st + TO_RjToBr);
        eSO3_BrToBr0 = so3Spline.Evaluate(et + TO_RjToBr);
    }

protected:
    static Eigen::Vector3d AngularVelBrToBr0(const So3SplineType &so3Spline, double t) {
        return so3Spline.Evaluate(t) * so3Spline.VelocityBody(t);
    }

    static Eigen::Vector3d RadarVelocityFromStaticTargetArray(const RadarTargetArray::Ptr &ary,
                                                              double ransacThd) {
        opengv::sac::Ransac<RadarVelocitySacProblem> ransac;
        std::shared_ptr<RadarVelocitySacProblem> probPtr(new RadarVelocitySacProblem(ary));
        ransac.sac_model_ = probPtr;
        ransac.threshold_ = ransacThd;
        ransac.max_iterations_ = 20;
        bool res = ransac.computeModel();
        if (res) {
            // spdlog::info("inlier rate: {}/{}", ransac.inliers_.size(), ary->GetTargets().size());
            Eigen::Vector3d radarVel;
            probPtr->optimizeModelCoefficients(ransac.inliers_, ransac.model_coefficients_,
                                               radarVel);

            // std::cout << " RANSAC: " << radarVel.transpose() << std::endl;
            // std::cout << "LS PROB: " << ary->RadarVelocityFromStaticTargetArray().transpose() <<
            // std::endl;

            return radarVel;
        } else {
            spdlog::warn("compute velocity using RANSAC failed, try to use all targets to fit...");
            return ary->RadarVelocityFromStaticTargetArray();
        }
    }
};

extern template struct RadarInertialAlignHelper<Configor::Prior::SplineOrder>;

template <int Order>
struct RadarInertialAlignFactor {
private:
    RadarInertialAlignHelper<Order> _helper;

    double _weight;

public:
    explicit RadarInertialAlignFactor(const RadarInertialAlignHelper<Order> &helper, double weight)
        : _helper(helper),
          _weight(weight) {}

    static auto Create(const RadarInertialAlignHelper<Order> &helper, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<RadarInertialAlignFactor>(
            new RadarInertialAlignFactor(helper, weight));
    }

    static std::size_t TypeHashCode() { return typeid(RadarInertialAlignFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ POS_BiInBr | SO3_RjToBr | POS_RjInBr | GRAVITY ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        Eigen::Map<const Eigen::Vector3<T>> POS_BiInBr(sKnots[0]);
        Eigen::Map<Sophus::SO3<T> const> const SO3_RjToBr(sKnots[1]);
        Eigen::Map<const Eigen::Vector3<T>> POS_RjInBr(sKnots[2]);
        Eigen::Map<const Eigen::Vector3<T>> GRAVITY(sKnots[3]);

        Eigen::Vector3<T> sLIN_VEL_BrToBr0InBr0 =
            _helper.sSO3_BrToBr0 * SO3_RjToBr * _helper.sDVec +
            Sophus::SO3<T>::hat(_helper.sSO3_BrToBr0 * POS_RjInBr) * _helper.sANG_VEL_BcToBc0;

        Eigen::Vector3<T> eLIN_VEL_BrToBr0InBr0 =
            _helper.eSO3_BrToBr0 * SO3_RjToBr * _helper.eDVec +
            Sophus::SO3<T>::hat(_helper.eSO3_BrToBr0 * POS_RjInBr) * _helper.eANG_VEL_BcToBc0;

        Eigen::Vector3<T> pred =
            eLIN_VEL_BrToBr0InBr0 - sLIN_VEL_BrToBr0InBr0 - GRAVITY * _helper.dt;
        Eigen::Vector3<T> mes =
            _helper.velVec.template cast<T>() - _helper.velMat.template cast<T>() * POS_BiInBr;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals.template block<3, 1>(0, 0) = T(_weight) * (pred - mes);

        return true;
    }
};

extern template struct RadarInertialAlignFactor<Configor::Prior::SplineOrder>;

}  // namespace ns_ikalibr
#endif  // IKALIBR_RADAR_INERTIAL_ALIGN_FACTOR_HPP
