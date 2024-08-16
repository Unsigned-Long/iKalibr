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

#ifndef IKALIBR_ESTIMATOR_H
#define IKALIBR_ESTIMATOR_H

#include "config/configor.h"
#include "ctraj/core/spline_bundle.h"
#include "ctraj/core/pose.hpp"
#include "calib/calib_param_manager.h"
#include "calib/calib_data_manager.h"
#include "factor/imu_acce_factor.hpp"
#include "factor/radar_factor.hpp"
#include "factor/lin_scale_factor.hpp"
#include "factor/point_to_surfel_factor.hpp"
#include "factor/visual_reproj_factor.hpp"
#include "factor/rgbd_velocity_factor.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
using namespace magic_enum::bitwise_operators;

struct OptOption {
    // myenumGenor Option OPT_SO3_SPLINE OPT_SCALE_SPLINE OPT_SO3_BiToBr OPT_POS_BiInBr
    // OPT_SO3_RjToBr OPT_POS_RjInBr OPT_SO3_LkToBr OPT_POS_LkInBr OPT_SO3_CmToBr OPT_POS_CmInBr
    // OPT_SO3_DnToBr OPT_POS_DnInBr OPT_TO_BiToBr OPT_TO_RjToBr OPT_TO_LkToBr OPT_TO_CmToBr
    // OPT_TO_DnToBr OPT_GYRO_BIAS OPT_GYRO_MAP_COEFF OPT_ACCE_BIAS OPT_ACCE_MAP_COEFF OPT_SO3_AtoG
    // OPT_GRAVITY OPT_VISUAL_GLOBAL_SCALE OPT_VISUAL_INV_DEPTH OPT_RGBD_DEPTH OPT_RGBD_ALPHA
    // OPT_RGBD_BETA OPT_CAM_FOCAL_LEN OPT_CAM_PRINCIPAL_POINT OPT_RS_CAM_READOUT_TIME
    enum Option : long {
        /**
         * @brief options
         */
        NONE = 1 << 0,

        OPT_SO3_SPLINE = 1 << 1,
        OPT_SCALE_SPLINE = 1 << 2,

        OPT_SO3_BiToBr = 1 << 3,
        OPT_POS_BiInBr = 1 << 4,

        OPT_SO3_RjToBr = 1 << 5,
        OPT_POS_RjInBr = 1 << 6,

        OPT_SO3_LkToBr = 1 << 7,
        OPT_POS_LkInBr = 1 << 8,

        OPT_SO3_CmToBr = 1 << 9,
        OPT_POS_CmInBr = 1 << 10,

        OPT_SO3_DnToBr = 1 << 11,
        OPT_POS_DnInBr = 1 << 12,

        OPT_TO_BiToBr = 1 << 13,
        OPT_TO_RjToBr = 1 << 14,
        OPT_TO_LkToBr = 1 << 15,
        OPT_TO_CmToBr = 1 << 16,
        OPT_TO_DnToBr = 1 << 17,

        OPT_GYRO_BIAS = 1 << 18,
        OPT_GYRO_MAP_COEFF = 1 << 19,
        OPT_ACCE_BIAS = 1 << 20,
        OPT_ACCE_MAP_COEFF = 1 << 21,
        OPT_SO3_AtoG = 1 << 22,

        OPT_GRAVITY = 1 << 23,

        OPT_VISUAL_GLOBAL_SCALE = 1 << 24,
        OPT_VISUAL_INV_DEPTH = 1 << 25,

        OPT_RGBD_DEPTH = 1 << 26,
        OPT_RGBD_ALPHA = 1 << 27,
        OPT_RGBD_BETA = 1 << 28,

        OPT_CAM_FOCAL_LEN = 1 << 29,
        OPT_CAM_PRINCIPAL_POINT = 1 << 30,
        OPT_RS_CAM_READOUT_TIME = 1 << 31,

        ALL = OPT_SO3_SPLINE | OPT_SCALE_SPLINE | OPT_SO3_BiToBr | OPT_POS_BiInBr | OPT_SO3_RjToBr |
              OPT_POS_RjInBr | OPT_SO3_LkToBr | OPT_POS_LkInBr | OPT_SO3_CmToBr | OPT_POS_CmInBr |
              OPT_SO3_DnToBr | OPT_POS_DnInBr | OPT_TO_BiToBr | OPT_TO_RjToBr | OPT_TO_LkToBr |
              OPT_TO_CmToBr | OPT_TO_DnToBr | OPT_GYRO_BIAS | OPT_GYRO_MAP_COEFF | OPT_ACCE_BIAS |
              OPT_ACCE_MAP_COEFF | OPT_SO3_AtoG | OPT_GRAVITY | OPT_VISUAL_GLOBAL_SCALE |
              OPT_VISUAL_INV_DEPTH | OPT_RGBD_ALPHA | OPT_RGBD_BETA | OPT_CAM_FOCAL_LEN |
              OPT_CAM_PRINCIPAL_POINT | OPT_RS_CAM_READOUT_TIME
    };
};

struct TimeDeriv {
    enum ScaleSplineType : int { LIN_ACCE_SPLINE = 0, LIN_VEL_SPLINE = 1, LIN_POS_SPLINE = 2 };

    enum ScaleType : int {
        LIN_ACCE = 0,
        LIN_VEL = -1,
        LIN_POS = -2,
    };

    template <ScaleSplineType spType, ScaleType sType>
    static constexpr int Deriv() {
        // compute how many times to perform time derivation to obtain "sType" from "spType"
        return spType + sType;
    }
};

struct SpatialTemporalPriori;
using SpatialTemporalPrioriPtr = std::shared_ptr<SpatialTemporalPriori>;

class Estimator : public ceres::Problem {
public:
    using Ptr = std::shared_ptr<Estimator>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;
    using SplineMetaType = ns_ctraj::SplineMeta<Configor::Prior::SplineOrder>;
    using Opt = OptOption::Option;

private:
    SplineBundleType::Ptr splines;
    CalibParamManager::Ptr parMagr;

    // manifolds
    static std::shared_ptr<ceres::EigenQuaternionManifold> QUATER_MANIFOLD;
    static std::shared_ptr<ceres::SphereManifold<3>> GRAVITY_MANIFOLD;

public:
    Estimator(SplineBundleType::Ptr splines, CalibParamManager::Ptr calibParamManager);

    static Ptr Create(const SplineBundleType::Ptr &splines,
                      const CalibParamManager::Ptr &calibParamManager);

    static ceres::Problem::Options DefaultProblemOptions();

    static ceres::Solver::Options DefaultSolverOptions(int threadNum = -1,
                                                       bool toStdout = true,
                                                       bool useCUDA = false);

    ceres::Solver::Summary Solve(
        const ceres::Solver::Options &options = Estimator::DefaultSolverOptions(),
        const SpatialTemporalPrioriPtr &priori = nullptr);

    Eigen::MatrixXd GetHessianMatrix(const std::vector<double *> &consideredParBlocks,
                                     int numThread = 1);

public:
    void AddIMUGyroMeasurement(const IMUFrame::Ptr &imuFrame,
                               const std::string &topic,
                               Opt option,
                               double gyroWeight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | ACCE_BIAS | ACCE_MAP_COEFF | GRAVITY |
     *   SO3_BiToBr | POS_BiInBr | TO_BiToBr ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddIMUAcceMeasurement(const IMUFrame::Ptr &imuFrame,
                               const std::string &topic,
                               Opt option,
                               double acceWeight) {
        // prepare metas for splines
        SplineMetaType so3Meta, scaleMeta;

        // different relative control points finding [single vs. range]
        // for the inertial measurements from the reference IMU, there is no need to consider a time
        // padding, as its time offsets would be fixed as identity
        if (IsOptionWith(Opt::OPT_TO_BiToBr, option) && Configor::DataStream::ReferIMU != topic) {
            double minTime = imuFrame->GetTimestamp() - Configor::Prior::TimeOffsetPadding;
            double maxTime = imuFrame->GetTimestamp() + Configor::Prior::TimeOffsetPadding;
            // invalid time stamp
            if (!splines->TimeInRangeForSo3(minTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(maxTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForRd(minTime, Configor::Preference::SCALE_SPLINE) ||
                !splines->TimeInRangeForRd(maxTime, Configor::Preference::SCALE_SPLINE)) {
                return;
            }
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{minTime, maxTime}},
                                            so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {{minTime, maxTime}},
                                           scaleMeta);
        } else {
            double curTime = imuFrame->GetTimestamp() + parMagr->TEMPORAL.TO_BiToBr.at(topic);

            // check point time stamp
            if (!splines->TimeInRangeForSo3(curTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForRd(curTime, Configor::Preference::SCALE_SPLINE)) {
                return;
            }
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{curTime, curTime}},
                                            so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {{curTime, curTime}},
                                           scaleMeta);
        }
        // create a cost function
        constexpr int derivIMU = TimeDeriv::Deriv<type, TimeDeriv::LIN_ACCE>();
        auto costFunc = IMUAcceFactor<Configor::Prior::SplineOrder, derivIMU>::Create(
            so3Meta, scaleMeta, imuFrame, acceWeight);

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }
        // pos knots param block [each has three sub params]
        for (int i = 0; i < static_cast<int>(scaleMeta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(3);
        }

        // ACCE_BIAS
        costFunc->AddParameterBlock(3);
        // ACCE_MAP_COEFF
        costFunc->AddParameterBlock(6);
        // GRAVITY
        costFunc->AddParameterBlock(3);
        // SO3_BiToBr
        costFunc->AddParameterBlock(4);
        // POS_BiInBr
        costFunc->AddParameterBlock(3);
        // TO_BiToBr
        costFunc->AddParameterBlock(1);

        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE),
                        so3Meta, !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

        // lin acce knots
        AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                       scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

        // ACCE_BIAS
        auto acceBias = parMagr->INTRI.IMU.at(topic)->ACCE.BIAS.data();
        paramBlockVec.push_back(acceBias);
        // ACCE_MAP_COEFF
        auto aceMapCoeff = parMagr->INTRI.IMU.at(topic)->ACCE.MAP_COEFF.data();
        paramBlockVec.push_back(aceMapCoeff);
        // GRAVITY
        auto gravity = parMagr->GRAVITY.data();
        paramBlockVec.push_back(gravity);
        // SO3_BiToBc
        auto SO3_BiToBc = parMagr->EXTRI.SO3_BiToBr.at(topic).data();
        paramBlockVec.push_back(SO3_BiToBc);
        // POS_BiInBc
        auto POS_BiInBc = parMagr->EXTRI.POS_BiInBr.at(topic).data();
        paramBlockVec.push_back(POS_BiInBc);
        // TIME_OFFSET_BiToBc
        auto TIME_OFFSET_BiToBc = &parMagr->TEMPORAL.TO_BiToBr.at(topic);
        paramBlockVec.push_back(TIME_OFFSET_BiToBc);

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);
        this->SetManifold(gravity, GRAVITY_MANIFOLD.get());
        this->SetManifold(SO3_BiToBc, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_ACCE_BIAS, option)) {
            this->SetParameterBlockConstant(acceBias);
        }

        if (!IsOptionWith(Opt::OPT_ACCE_MAP_COEFF, option)) {
            this->SetParameterBlockConstant(aceMapCoeff);
        }

        if (!IsOptionWith(Opt::OPT_GRAVITY, option)) {
            this->SetParameterBlockConstant(gravity);
        }

        if (!IsOptionWith(Opt::OPT_SO3_BiToBr, option)) {
            this->SetParameterBlockConstant(SO3_BiToBc);
        }

        if (!IsOptionWith(Opt::OPT_POS_BiInBr, option)) {
            this->SetParameterBlockConstant(POS_BiInBc);
        }

        if (!IsOptionWith(Opt::OPT_TO_BiToBr, option)) {
            this->SetParameterBlockConstant(TIME_OFFSET_BiToBc);
        } else {
            // set bound
            this->SetParameterLowerBound(TIME_OFFSET_BiToBc, 0,
                                         -Configor::Prior::TimeOffsetPadding);
            this->SetParameterUpperBound(TIME_OFFSET_BiToBc, 0, Configor::Prior::TimeOffsetPadding);
        }
    }

    void AddInertialAlignment(const std::vector<IMUFrame::Ptr> &data,
                              const std::string &imuTopic,
                              double sTimeByBr,
                              double eTimeByBr,
                              Eigen::Vector3d *sVel,
                              Eigen::Vector3d *eVel,
                              Opt option,
                              double weight);

    void AddLiDARInertialAlignment(const std::vector<IMUFrame::Ptr> &data,
                                   const std::string &lidarTopic,
                                   const std::string &imuTopic,
                                   const ns_ctraj::Posed &sPose,
                                   const ns_ctraj::Posed &ePose,
                                   double mapTime,
                                   Eigen::Vector3d *sVel,
                                   Eigen::Vector3d *eVel,
                                   Estimator::Opt option,
                                   double weight);

    void AddVisualInertialAlignment(const std::vector<IMUFrame::Ptr> &data,
                                    const std::string &camTopic,
                                    const std::string &imuTopic,
                                    const ns_ctraj::Posed &sPose,
                                    const ns_ctraj::Posed &ePose,
                                    double mapTime,
                                    Eigen::Vector3d *sVel,
                                    Eigen::Vector3d *eVel,
                                    double *SCALE,
                                    Estimator::Opt option,
                                    double weight);

    void AddRadarInertialAlignment(const std::vector<IMUFrame::Ptr> &data,
                                   const std::string &imuTopic,
                                   const std::string &radarTopic,
                                   const RadarTargetArray::Ptr &sRadarAry,
                                   const RadarTargetArray::Ptr &eRadarAry,
                                   Estimator::Opt option,
                                   double weight);

    void AddRGBDInertialAlignment(const std::vector<IMUFrame::Ptr> &data,
                                  const std::string &imuTopic,
                                  const std::string &rgbdTopic,
                                  const std::pair<CameraFrame::Ptr, Eigen::Vector3d> &sRGBDAry,
                                  const std::pair<CameraFrame::Ptr, Eigen::Vector3d> &eRGBDAry,
                                  Estimator::Opt option,
                                  double weight);

    void AddRadarInertialRotRoughAlignment(const std::vector<IMUFrame::Ptr> &data,
                                           const std::string &imuTopic,
                                           const std::string &radarTopic,
                                           const RadarTargetArray::Ptr &sRadarAry,
                                           const RadarTargetArray::Ptr &eRadarAry,
                                           Estimator::Opt option,
                                           double weight);

    void AddHandEyeRotationAlignmentForLiDAR(const std::string &lidarTopic,
                                             double tLastByLk,
                                             double tCurByLk,
                                             const Sophus::SO3d &so3LastLkToM,
                                             const Sophus::SO3d &so3CurLkToM,
                                             Estimator::Opt option,
                                             double weight);

    void AddHandEyeRotationAlignmentForCamera(const std::string &camTopic,
                                              double tLastByCm,
                                              double tCurByCm,
                                              const Sophus::SO3d &so3LastCmToW,
                                              const Sophus::SO3d &so3CurCmToW,
                                              Estimator::Opt option,
                                              double weight);

    void AddHandEyeRotationAlignmentForRGBD(const std::string &rgbdTopic,
                                            double tLastByDn,
                                            double tCurByDn,
                                            const Sophus::SO3d &so3LastDnToW,
                                            const Sophus::SO3d &so3CurDnToW,
                                            Estimator::Opt option,
                                            double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_RjToBr | POS_RjInBr | TO_RjToBr ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddRadarMeasurement(const RadarTarget::Ptr &radarFrame,
                             const std::string &topic,
                             Estimator::Opt option,
                             double weight) {
        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        const auto &scaleSpline = splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

        SplineMetaType so3Meta, scaleMeta;

        // different relative control points finding [single vs. range]
        if (IsOptionWith(Opt::OPT_TO_RjToBr, option)) {
            double tMin = radarFrame->GetTimestamp() - Configor::Prior::TimeOffsetPadding;
            double tMax = radarFrame->GetTimestamp() + Configor::Prior::TimeOffsetPadding;
            // invalid time stamp
            if (!splines->TimeInRange(tMin, so3Spline) || !splines->TimeInRange(tMax, so3Spline) ||
                !splines->TimeInRange(tMin, scaleSpline) ||
                !splines->TimeInRange(tMax, scaleSpline)) {
                return;
            }

            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{tMin, tMax}},
                                            so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {{tMin, tMax}},
                                           scaleMeta);
        } else {
            double t = radarFrame->GetTimestamp() + parMagr->TEMPORAL.TO_RjToBr.at(topic);

            // check point time stamp
            if (!splines->TimeInRange(t, so3Spline) || !splines->TimeInRange(t, scaleSpline)) {
                return;
            }

            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{t, t}}, so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {{t, t}}, scaleMeta);
        }

        static constexpr int derivRadar = TimeDeriv::Deriv<type, TimeDeriv::LIN_VEL>();

        // create a cost function
        auto costFunc = RadarFactor<Configor::Prior::SplineOrder, derivRadar>::Create(
            so3Meta, scaleMeta, radarFrame, weight);

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }
        // pos knots param block [each has three sub params]
        for (int i = 0; i < static_cast<int>(scaleMeta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(3);
        }
        costFunc->AddParameterBlock(4);  // SO3_RtoB
        costFunc->AddParameterBlock(3);  // POS_RinB
        costFunc->AddParameterBlock(1);  // TIME_OFFSET_RtoB

        // the Residual
        costFunc->SetNumResiduals(1);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE),
                        so3Meta, !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

        // lin acce knots
        AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                       scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

        auto SO3_RjToBr = parMagr->EXTRI.SO3_RjToBr.at(topic).data();
        paramBlockVec.push_back(SO3_RjToBr);

        auto POS_RjInBr = parMagr->EXTRI.POS_RjInBr.at(topic).data();
        paramBlockVec.push_back(POS_RjInBr);

        auto TO_RjToBr = &parMagr->TEMPORAL.TO_RjToBr.at(topic);
        paramBlockVec.push_back(TO_RjToBr);

        // pass to problem
        // remove dynamic targets (outliers)
        // this->AddResidualBlockToProblem(costFunc, new ceres::HuberLoss(weight * weight * 0.125),
        // paramBlockVec);
        this->AddResidualBlock(
            costFunc, new ceres::CauchyLoss(Configor::Prior::CauchyLossForRadarFactor * weight),
            paramBlockVec);
        this->SetManifold(SO3_RjToBr, QUATER_MANIFOLD.get());

        // lock param or not
        if (!IsOptionWith(Opt::OPT_TO_RjToBr, option)) {
            this->SetParameterBlockConstant(TO_RjToBr);
        } else {
            // set bound
            this->SetParameterLowerBound(TO_RjToBr, 0, -Configor::Prior::TimeOffsetPadding);
            this->SetParameterUpperBound(TO_RjToBr, 0, Configor::Prior::TimeOffsetPadding);
        }
        if (!IsOptionWith(Opt::OPT_SO3_RjToBr, option)) {
            this->SetParameterBlockConstant(SO3_RjToBr);
        }
        if (!IsOptionWith(Opt::OPT_POS_RjInBr, option)) {
            this->SetParameterBlockConstant(POS_RjInBr);
        }
    }

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_RjToBr | POS_RjInBr | TO_RjToBr ]
     */
    template <int TimeDeriv>
    void AddLinearScaleConstraint(double timeByBr,
                                  const Eigen::Vector3d &linScaleOfDeriv,
                                  Estimator::Opt option,
                                  double weight) {
        const auto &scaleSpline = splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
        // check point time stamp
        if (!scaleSpline.TimeStampInRange(timeByBr) || !scaleSpline.TimeStampInRange(timeByBr)) {
            return;
        }

        SplineMetaType scaleMeta;
        splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {{timeByBr, timeByBr}},
                                       scaleMeta);

        // create a cost function
        auto costFunc = LinearScaleDerivFactor<Configor::Prior::SplineOrder, TimeDeriv>::Create(
            scaleMeta, timeByBr, linScaleOfDeriv, weight);

        // pos knots param block [each has three sub params]
        for (int i = 0; i < static_cast<int>(scaleMeta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(3);
        }

        // the Residual
        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // lin acce knots
        AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                       scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);
    }

    void AddSO3Constraint(double timeByBr,
                          const Sophus::SO3d &so3,
                          Estimator::Opt option,
                          double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_LkToBr | POS_LkInBr | TO_LkToBr ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddPointTiSurfelConstraint(const PointToSurfelCorr::Ptr &ptsCorr,
                                    const std::string &topic,
                                    Opt option,
                                    double weight) {
        // prepare metas for splines
        SplineMetaType so3Meta, scaleMeta;

        // different relative control points finding [single vs. range]
        if (IsOptionWith(Opt::OPT_TO_LkToBr, option)) {
            double maxTime = ptsCorr->timestamp + Configor::Prior::TimeOffsetPadding;
            double minTime = ptsCorr->timestamp - Configor::Prior::TimeOffsetPadding;

            // invalid time stamp
            if (!splines->TimeInRangeForSo3(minTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(maxTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForRd(minTime, Configor::Preference::SCALE_SPLINE) ||
                !splines->TimeInRangeForRd(maxTime, Configor::Preference::SCALE_SPLINE)) {
                return;
            }

            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{minTime, maxTime}},
                                            so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {{minTime, maxTime}},
                                           scaleMeta);
        } else {
            double curTime = ptsCorr->timestamp + parMagr->TEMPORAL.TO_LkToBr.at(topic);

            // check point time stamp
            if (!splines->TimeInRangeForSo3(curTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForRd(curTime, Configor::Preference::SCALE_SPLINE)) {
                return;
            }

            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{curTime, curTime}},
                                            so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {{curTime, curTime}},
                                           scaleMeta);
        }
        static constexpr int derivLiDAR = TimeDeriv::Deriv<type, TimeDeriv::LIN_POS>();
        // create a cost function
        auto costFunc = PointToSurfelFactor<Configor::Prior::SplineOrder, derivLiDAR>::Create(
            so3Meta, scaleMeta, ptsCorr, weight);

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }
        // pos knots param block [each has three sub params]
        for (int i = 0; i < static_cast<int>(scaleMeta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(3);
        }

        costFunc->AddParameterBlock(4);
        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(1);

        costFunc->SetNumResiduals(1);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE),
                        so3Meta, !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

        // lin acce knots
        AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                       scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

        // SO3_LkToBr
        auto SO3_LkToBr = parMagr->EXTRI.SO3_LkToBr.at(topic).data();
        paramBlockVec.push_back(SO3_LkToBr);

        // POS_LkInBr
        auto POS_LkInBr = parMagr->EXTRI.POS_LkInBr.at(topic).data();
        paramBlockVec.push_back(POS_LkInBr);

        // TO_LkToBr
        auto TO_LkToBr = &parMagr->TEMPORAL.TO_LkToBr.at(topic);
        paramBlockVec.push_back(TO_LkToBr);

        // pass to problem
        this->AddResidualBlock(costFunc,
                               new ceres::CauchyLoss(Configor::Prior::CauchyLossForLiDARFactor),
                               paramBlockVec);
        this->SetManifold(SO3_LkToBr, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_SO3_LkToBr, option)) {
            this->SetParameterBlockConstant(SO3_LkToBr);
        }

        if (!IsOptionWith(Opt::OPT_POS_LkInBr, option)) {
            this->SetParameterBlockConstant(POS_LkInBr);
        }

        if (!IsOptionWith(Opt::OPT_TO_LkToBr, option)) {
            this->SetParameterBlockConstant(TO_LkToBr);
        } else {
            // set bound
            this->SetParameterLowerBound(TO_LkToBr, 0, -Configor::Prior::TimeOffsetPadding);
            this->SetParameterUpperBound(TO_LkToBr, 0, Configor::Prior::TimeOffsetPadding);
        }
    }

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
     * READOUT_TIME | FX | FY | CX | CY | GLOBAL_SCALE | INV_DEPTH ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddVisualReprojection(const VisualReProjCorr &visualCorr,
                               const std::string &topic,
                               double *globalScale,
                               double *invDepth,
                               Opt option,
                               double weight) {
        // prepare metas for splines
        SplineMetaType so3Meta, scaleMeta;

        const double TO_PADDING = Configor::Prior::TimeOffsetPadding;
        const double RT_PADDING = Configor::Prior::ReadoutTimePadding;
        double *RS_READOUT = &parMagr->TEMPORAL.RS_READOUT.at(topic);
        double *TO_CmToBr = &parMagr->TEMPORAL.TO_CmToBr.at(topic);

        double minTimeI = visualCorr.ti;
        double maxTimeI = visualCorr.ti;

        double minTimeJ = visualCorr.tj;
        double maxTimeJ = visualCorr.tj;

        // different relative control points finding [single vs. range]
        if (IsOptionWith(Opt::OPT_TO_CmToBr, option)) {
            minTimeI -= TO_PADDING;
            maxTimeI += TO_PADDING;

            minTimeJ -= TO_PADDING;
            maxTimeJ += TO_PADDING;
        } else {
            minTimeI += *TO_CmToBr;
            maxTimeI += *TO_CmToBr;

            minTimeJ += *TO_CmToBr;
            maxTimeJ += *TO_CmToBr;
        }

        if (IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option)) {
            minTimeI += std::min(visualCorr.li * 0.0, visualCorr.li * RT_PADDING);
            maxTimeI += std::max(visualCorr.li * 0.0, visualCorr.li * RT_PADDING);

            minTimeJ += std::min(visualCorr.lj * 0.0, visualCorr.lj * RT_PADDING);
            maxTimeJ += std::max(visualCorr.lj * 0.0, visualCorr.lj * RT_PADDING);
        } else {
            minTimeI += visualCorr.li * *RS_READOUT;
            maxTimeI += visualCorr.li * *RS_READOUT;

            minTimeJ += visualCorr.lj * *RS_READOUT;
            maxTimeJ += visualCorr.lj * *RS_READOUT;
        }

        // invalid time stamp
        if (!splines->TimeInRangeForSo3(minTimeI, Configor::Preference::SO3_SPLINE) ||
            !splines->TimeInRangeForSo3(maxTimeI, Configor::Preference::SO3_SPLINE) ||
            !splines->TimeInRangeForRd(minTimeI, Configor::Preference::SCALE_SPLINE) ||
            !splines->TimeInRangeForRd(maxTimeI, Configor::Preference::SCALE_SPLINE)) {
            return;
        }

        // invalid time stamp
        if (!splines->TimeInRangeForSo3(minTimeJ, Configor::Preference::SO3_SPLINE) ||
            !splines->TimeInRangeForSo3(maxTimeJ, Configor::Preference::SO3_SPLINE) ||
            !splines->TimeInRangeForRd(minTimeJ, Configor::Preference::SCALE_SPLINE) ||
            !splines->TimeInRangeForRd(maxTimeJ, Configor::Preference::SCALE_SPLINE)) {
            return;
        }

        std::pair<double, double> timePairI = {minTimeI, maxTimeI};
        std::pair<double, double> timePairJ = {minTimeJ, maxTimeJ};

        if (minTimeI < minTimeJ) {
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                            {timePairI, timePairJ}, so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE,
                                           {timePairI, timePairJ}, scaleMeta);
        } else {
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                            {timePairJ, timePairI}, so3Meta);
            splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE,
                                           {timePairJ, timePairI}, scaleMeta);
        }

        static constexpr int deriv = TimeDeriv::Deriv<type, TimeDeriv::LIN_POS>();
        // create a cost function
        auto costFunc = VisualReProjFactor<Configor::Prior::SplineOrder, deriv>::Create(
            so3Meta, scaleMeta, visualCorr, weight);

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }
        // pos knots param block [each has three sub params]
        for (int i = 0; i < static_cast<int>(scaleMeta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(3);
        }

        costFunc->AddParameterBlock(4);
        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);

        // fx, fy, cx, cy
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);

        // global scale, inv depth
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);

        costFunc->SetNumResiduals(2);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE),
                        so3Meta, !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

        // lin acce knots
        AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                       scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

        auto SO3_CmToBr = parMagr->EXTRI.SO3_CmToBr.at(topic).data();
        paramBlockVec.push_back(SO3_CmToBr);

        auto POS_CmInBr = parMagr->EXTRI.POS_CmInBr.at(topic).data();
        paramBlockVec.push_back(POS_CmInBr);

        paramBlockVec.push_back(TO_CmToBr);
        paramBlockVec.push_back(RS_READOUT);

        auto &intri = parMagr->INTRI.Camera.at(topic);
        paramBlockVec.push_back(intri->FXAddress());
        paramBlockVec.push_back(intri->FYAddress());
        paramBlockVec.push_back(intri->CXAddress());
        paramBlockVec.push_back(intri->CYAddress());

        paramBlockVec.push_back(globalScale);
        paramBlockVec.push_back(invDepth);

        // pass to problem
        this->AddResidualBlock(costFunc,
                               new ceres::CauchyLoss(Configor::Prior::CauchyLossForCameraFactor),
                               paramBlockVec);
        this->SetManifold(SO3_CmToBr, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_SO3_CmToBr, option)) {
            this->SetParameterBlockConstant(SO3_CmToBr);
        }

        if (!IsOptionWith(Opt::OPT_POS_CmInBr, option)) {
            this->SetParameterBlockConstant(POS_CmInBr);
        }

        if (!IsOptionWith(Opt::OPT_TO_CmToBr, option)) {
            this->SetParameterBlockConstant(TO_CmToBr);
        } else {
            // set bound
            this->SetParameterLowerBound(TO_CmToBr, 0, -TO_PADDING);
            this->SetParameterUpperBound(TO_CmToBr, 0, TO_PADDING);
        }

        if (!IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option)) {
            this->SetParameterBlockConstant(RS_READOUT);
        } else {
            // set bound
            this->SetParameterLowerBound(RS_READOUT, 0, 0.0);
            this->SetParameterUpperBound(RS_READOUT, 0, RT_PADDING);
        }

        if (!IsOptionWith(Opt::OPT_CAM_FOCAL_LEN, option)) {
            this->SetParameterBlockConstant(intri->FXAddress());
            this->SetParameterBlockConstant(intri->FYAddress());
        }

        if (!IsOptionWith(Opt::OPT_CAM_PRINCIPAL_POINT, option)) {
            this->SetParameterBlockConstant(intri->CXAddress());
            this->SetParameterBlockConstant(intri->CYAddress());
        }

        if (!IsOptionWith(Opt::OPT_VISUAL_GLOBAL_SCALE, option)) {
            this->SetParameterBlockConstant(globalScale);
        } else {
            // set bound
            this->SetParameterLowerBound(globalScale, 0, 1E-3);
        }

        if (!IsOptionWith(Opt::OPT_VISUAL_INV_DEPTH, option)) {
            this->SetParameterBlockConstant(invDepth);
        } else {
            // set bound
            this->SetParameterLowerBound(invDepth, 0, 1E-3);
        }
    }

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_DnToBr | POS_DnInBr | TO_DnToBr |
     *   READOUT_TIME | FX | FY | CX | CY | ALPHA | BETA | DEPTH ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddRGBDVelocityConstraint(const RGBDVelocityCorr::Ptr &velCorr,
                                   const std::string &topic,
                                   Opt option,
                                   double weight) {
        // invalid depth
        if (velCorr->depth < 1E-3) {
            return;
        }
        // prepare metas for splines
        SplineMetaType so3Meta, scaleMeta;

        const double TO_PADDING = Configor::Prior::TimeOffsetPadding;
        const double RT_PADDING = Configor::Prior::ReadoutTimePadding;
        double *RS_READOUT = &parMagr->TEMPORAL.RS_READOUT.at(topic);
        double *TO_DnToBr = &parMagr->TEMPORAL.TO_DnToBr.at(topic);

        // readout time equals to 0.0 means the raw middle time
        double minTime = velCorr->MidPointTime(0.0);
        double maxTime = velCorr->MidPointTime(0.0);

        // different relative control points finding [single vs. range]
        if (IsOptionWith(Opt::OPT_TO_DnToBr, option)) {
            minTime -= TO_PADDING;
            maxTime += TO_PADDING;
        } else {
            minTime += *TO_DnToBr;
            maxTime += *TO_DnToBr;
        }

        const double midRDFactor = velCorr->MidReadoutFactor();
        if (IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option)) {
            minTime += std::min(midRDFactor * 0.0, midRDFactor * RT_PADDING);
            maxTime += std::max(midRDFactor * 0.0, midRDFactor * RT_PADDING);
        } else {
            minTime += midRDFactor * *RS_READOUT;
            maxTime += midRDFactor * *RS_READOUT;
        }

        // invalid time stamp
        if (!splines->TimeInRangeForSo3(minTime, Configor::Preference::SO3_SPLINE) ||
            !splines->TimeInRangeForSo3(maxTime, Configor::Preference::SO3_SPLINE) ||
            !splines->TimeInRangeForRd(minTime, Configor::Preference::SCALE_SPLINE) ||
            !splines->TimeInRangeForRd(maxTime, Configor::Preference::SCALE_SPLINE)) {
            return;
        }

        std::pair<double, double> timePair = {minTime, maxTime};
        splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {timePair}, so3Meta);
        splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {timePair}, scaleMeta);

        static constexpr int deriv = TimeDeriv::Deriv<type, TimeDeriv::LIN_VEL>();
        // create a cost function
        auto costFunc = RGBDVelocityFactor<Configor::Prior::SplineOrder, deriv>::Create(
            so3Meta, scaleMeta, velCorr, weight);

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }
        // pos knots param block [each has three sub params]
        for (int i = 0; i < static_cast<int>(scaleMeta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(3);
        }

        costFunc->AddParameterBlock(4);
        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);

        // fx, fy, cx, cy
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);

        // alpha, beta, depth
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);
        costFunc->AddParameterBlock(1);

        costFunc->SetNumResiduals(2);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE),
                        so3Meta, !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

        // lin acce knots
        AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                       scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

        auto SO3_DnToBr = parMagr->EXTRI.SO3_DnToBr.at(topic).data();
        paramBlockVec.push_back(SO3_DnToBr);

        auto POS_DnInBr = parMagr->EXTRI.POS_DnInBr.at(topic).data();
        paramBlockVec.push_back(POS_DnInBr);

        paramBlockVec.push_back(TO_DnToBr);
        paramBlockVec.push_back(RS_READOUT);

        auto &intri = parMagr->INTRI.RGBD.at(topic);
        paramBlockVec.push_back(intri->intri->FXAddress());
        paramBlockVec.push_back(intri->intri->FYAddress());
        paramBlockVec.push_back(intri->intri->CXAddress());
        paramBlockVec.push_back(intri->intri->CYAddress());

        paramBlockVec.push_back(&intri->alpha);
        paramBlockVec.push_back(&intri->beta);
        paramBlockVec.push_back(&velCorr->depth);

        // pass to problem
        this->AddResidualBlock(costFunc,
                               new ceres::CauchyLoss(Configor::Prior::CauchyLossForRGBDFactor),
                               paramBlockVec);
        this->SetManifold(SO3_DnToBr, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_SO3_DnToBr, option)) {
            this->SetParameterBlockConstant(SO3_DnToBr);
        }

        if (!IsOptionWith(Opt::OPT_POS_DnInBr, option)) {
            this->SetParameterBlockConstant(POS_DnInBr);
        }

        if (!IsOptionWith(Opt::OPT_TO_DnToBr, option)) {
            this->SetParameterBlockConstant(TO_DnToBr);
        } else {
            // set bound
            this->SetParameterLowerBound(TO_DnToBr, 0, -TO_PADDING);
            this->SetParameterUpperBound(TO_DnToBr, 0, TO_PADDING);
        }

        if (!IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option)) {
            this->SetParameterBlockConstant(RS_READOUT);
        } else {
            // set bound
            this->SetParameterLowerBound(RS_READOUT, 0, 0.0);
            this->SetParameterUpperBound(RS_READOUT, 0, RT_PADDING);
        }

        if (!IsOptionWith(Opt::OPT_CAM_FOCAL_LEN, option)) {
            this->SetParameterBlockConstant(intri->intri->FXAddress());
            this->SetParameterBlockConstant(intri->intri->FYAddress());
        }

        if (!IsOptionWith(Opt::OPT_CAM_PRINCIPAL_POINT, option)) {
            this->SetParameterBlockConstant(intri->intri->CXAddress());
            this->SetParameterBlockConstant(intri->intri->CYAddress());
        }

        if (!IsOptionWith(Opt::OPT_RGBD_ALPHA, option)) {
            this->SetParameterBlockConstant(&intri->alpha);
        }

        if (!IsOptionWith(Opt::OPT_RGBD_BETA, option)) {
            this->SetParameterBlockConstant(&intri->beta);
        }

        if (!IsOptionWith(Opt::OPT_RGBD_DEPTH, option)) {
            this->SetParameterBlockConstant(&velCorr->depth);
        }
    }

    void SetRefIMUParamsConstant();

    void FixFirSO3ControlPoint();

    void AddVisualProjectionFactor(ns_veta::Posed *T_CurCToW,
                                   Eigen::Vector3d *POS_LMInW,
                                   const ns_veta::PinholeIntrinsic::Ptr &intri,
                                   const Eigen::Vector2d &feat,
                                   double weight);

    void AddLinScaleTailConstraint(Opt option,
                                   double weight,
                                   int count = Configor::Prior::SplineOrder);

    void AddSO3TailConstraint(Opt option, double weight, int count = Configor::Prior::SplineOrder);

    void AddLinScaleHeadConstraint(Opt option,
                                   double weight,
                                   int count = Configor::Prior::SplineOrder);

    void AddSO3HeadConstraint(Opt option, double weight, int count = Configor::Prior::SplineOrder);

    void AddPriorExtriSO3Constraint(const Sophus::SO3d &SO3_Sen1ToSen2,
                                    Sophus::SO3d *SO3_Sen1ToRef,
                                    Sophus::SO3d *SO3_Sen2ToRef,
                                    double weight);

    void AddPriorExtriPOSConstraint(const Eigen::Vector3d &POS_Sen1InSen2,
                                    Eigen::Vector3d *POS_Sen1InRef,
                                    Sophus::SO3d *SO3_Sen2ToRef,
                                    Eigen::Vector3d *POS_Sen2InRef,
                                    double weight);

    void AddPriorTimeOffsetConstraint(const double &TO_Sen1ToSen2,
                                      double *TO_Sen1ToRef,
                                      double *TO_Sen2ToRef,
                                      double weight);

protected:
    void AddSo3KnotsData(std::vector<double *> &paramBlockVec,
                         const SplineBundleType::So3SplineType &spline,
                         const SplineMetaType &splineMeta,
                         bool setToConst);

    void AddRdKnotsData(std::vector<double *> &paramBlockVec,
                        const SplineBundleType::RdSplineType &spline,
                        const SplineMetaType &splineMeta,
                        bool setToConst);

    static Eigen::MatrixXd CRSMatrix2EigenMatrix(ceres::CRSMatrix *jacobian_crs_matrix);

    std::pair<Eigen::Vector3d, Eigen::Matrix3d> InertialVelIntegration(
        const std::vector<IMUFrame::Ptr> &data,
        const std::string &imuTopic,
        double sTimeByBi,
        double eTimeByBi);

    std::pair<std::pair<Eigen::Vector3d, Eigen::Matrix3d>,
              std::pair<Eigen::Vector3d, Eigen::Matrix3d>>
    InertialPosIntegration(const std::vector<IMUFrame::Ptr> &data,
                           const std::string &imuTopic,
                           double sTimeByBi,
                           double eTimeByBi);

    std::pair<std::vector<std::pair<double, Eigen::Vector3d>>,
              std::vector<std::pair<double, Eigen::Matrix3d>>>
    InertialIntegrationBase(const std::vector<IMUFrame::Ptr> &data,
                            const std::string &imuTopic,
                            double sTimeByBi,
                            double eTimeByBi);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_ESTIMATOR_H
