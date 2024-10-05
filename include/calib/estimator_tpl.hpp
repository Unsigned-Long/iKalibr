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

#ifndef IKALIBR_ESTIMATOR_TPL_HPP
#define IKALIBR_ESTIMATOR_TPL_HPP

#include "calib/estimator.h"
#include "factor/imu_acce_factor.hpp"
#include "factor/lin_scale_factor.hpp"
#include "factor/point_to_surfel_factor.hpp"
#include "factor/radar_factor.hpp"
#include "factor/visual_optical_flow_factor.hpp"
#include "factor/visual_reproj_factor.hpp"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | ACCE_BIAS | ACCE_MAP_COEFF | GRAVITY |
 *   SO3_BiToBr | POS_BiInBr | TO_BiToBr ]
 */
template <TimeDeriv::ScaleSplineType type>
void Estimator::AddIMUAcceMeasurement(const IMUFrame::Ptr &imuFrame,
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
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

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
        this->SetParameterLowerBound(TIME_OFFSET_BiToBc, 0, -Configor::Prior::TimeOffsetPadding);
        this->SetParameterUpperBound(TIME_OFFSET_BiToBc, 0, Configor::Prior::TimeOffsetPadding);
    }
}

/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_RjToBr | POS_RjInBr | TO_RjToBr ]
 */
template <TimeDeriv::ScaleSplineType type>
void Estimator::AddRadarMeasurement(const RadarTarget::Ptr &radarFrame,
                                    const std::string &topic,
                                    Opt option,
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
            !splines->TimeInRange(tMin, scaleSpline) || !splines->TimeInRange(tMax, scaleSpline)) {
            return;
        }

        splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{tMin, tMax}}, so3Meta);
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
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

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
        costFunc, new ceres::HuberLoss(Configor::Prior::LossForRadarDopplerFactor * weight),
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
void Estimator::AddLinearScaleConstraint(double timeByBr,
                                         const Eigen::Vector3d &linScaleOfDeriv,
                                         Opt option,
                                         double weight) {
    const auto &scaleSpline = splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    // check point time stamp
    if (!scaleSpline.TimeStampInRange(timeByBr)) {
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

/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_LkToBr | POS_LkInBr | TO_LkToBr ]
 */
template <TimeDeriv::ScaleSplineType type>
void Estimator::AddLiDARPointTiSurfelConstraint(const PointToSurfelCorrPtr &ptsCorr,
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
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

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
    this->AddResidualBlock(
        costFunc, new ceres::HuberLoss(Configor::Prior::LossForPointToSurfelFactor * weight),
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
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_DnToBr | POS_DnInBr | TO_DnToBr ]
 */
template <TimeDeriv::ScaleSplineType type>
void Estimator::AddRGBDPointTiSurfelConstraint(const PointToSurfelCorrPtr &ptsCorr,
                                               const std::string &topic,
                                               Opt option,
                                               double weight) {
    // prepare metas for splines
    SplineMetaType so3Meta, scaleMeta;

    // different relative control points finding [single vs. range]
    if (IsOptionWith(Opt::OPT_TO_DnToBr, option)) {
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
        double curTime = ptsCorr->timestamp + parMagr->TEMPORAL.TO_DnToBr.at(topic);

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
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

    // lin acce knots
    AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                   scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

    // SO3_DnToBr
    auto SO3_DnToBr = parMagr->EXTRI.SO3_DnToBr.at(topic).data();
    paramBlockVec.push_back(SO3_DnToBr);

    // POS_DnInBr
    auto POS_DnInBr = parMagr->EXTRI.POS_DnInBr.at(topic).data();
    paramBlockVec.push_back(POS_DnInBr);

    // TO_DnToBr
    auto TO_DnToBr = &parMagr->TEMPORAL.TO_DnToBr.at(topic);
    paramBlockVec.push_back(TO_DnToBr);

    // pass to problem
    this->AddResidualBlock(
        costFunc,
        // we use 'Configor::Prior::LossForLiDARFactor' for rgbds here
        new ceres::HuberLoss(Configor::Prior::LossForPointToSurfelFactor * weight), paramBlockVec);
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
        this->SetParameterLowerBound(TO_DnToBr, 0, -Configor::Prior::TimeOffsetPadding);
        this->SetParameterUpperBound(TO_DnToBr, 0, Configor::Prior::TimeOffsetPadding);
    }
}

/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
 * READOUT_TIME | FX | FY | CX | CY | GLOBAL_SCALE | INV_DEPTH ]
 */
template <TimeDeriv::ScaleSplineType type>
void Estimator::AddVisualReprojection(const VisualReProjCorr::Ptr &visualCorr,
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

    std::pair<double, double> timePairI = ConsideredTimeRangeForCameraStamp(
        visualCorr->ti,                                      // time stamped by the camera
        *RS_READOUT, RT_PADDING, visualCorr->li,             // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_CmToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_CmToBr, option)  // if opt time offset
    );
    std::pair<double, double> timePairJ = ConsideredTimeRangeForCameraStamp(
        visualCorr->tj,                                      // time stamped by the camera
        *RS_READOUT, RT_PADDING, visualCorr->lj,             // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_CmToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_CmToBr, option)  // if opt time offset
    );

    if (!TimeInRangeForSplines(timePairI) || !TimeInRangeForSplines(timePairJ)) {
        return;
    }

    if (timePairI.first < timePairJ.first) {
        splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {timePairI, timePairJ},
                                        so3Meta);
        splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {timePairI, timePairJ},
                                       scaleMeta);
    } else {
        splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {timePairJ, timePairI},
                                        so3Meta);
        splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {timePairJ, timePairI},
                                       scaleMeta);
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
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

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
                           new ceres::HuberLoss(Configor::Prior::LossForReprojFactor * weight),
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

    if (!IsOptionWith(Opt::OPT_VISUAL_DEPTH, option)) {
        this->SetParameterBlockConstant(invDepth);
    } else {
        // set bound
        this->SetParameterLowerBound(invDepth, 0, 1E-3);
    }
}

/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_DnToBr | POS_DnInBr | TO_DnToBr |
 *   READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
 */
template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void Estimator::AddRGBDOpticalFlowConstraint(const OpticalFlowCorr::Ptr &ofCorr,
                                             const std::string &topic,
                                             Opt option,
                                             double weight) {
    auto &intri = parMagr->INTRI.RGBD.at(topic);
    // invalid depth
    if (intri->ActualDepth(ofCorr->depth) < 1E-3) {
        return;
    }

    const double TO_PADDING = Configor::Prior::TimeOffsetPadding;
    const double RT_PADDING = Configor::Prior::ReadoutTimePadding;
    double *RS_READOUT = &parMagr->TEMPORAL.RS_READOUT.at(topic);
    double *TO_DnToBr = &parMagr->TEMPORAL.TO_DnToBr.at(topic);

    if (ofCorr->MidPointVel(*RS_READOUT).norm() < Configor::Prior::LossForOpticalFlowFactor) {
        // small pixel velocity
        return;
    }

    std::pair<double, double> timePair = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(1),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(1),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_DnToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_DnToBr, option)  // if opt time offset
    );

    if (!TimeInRangeForSplines(timePair)) {
        return;
    }

    // prepare metas for splines
    SplineMetaType so3Meta, scaleMeta;

    splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {timePair}, so3Meta);
    splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {timePair}, scaleMeta);

    static constexpr int deriv = TimeDeriv::Deriv<type, TimeDeriv::LIN_VEL>();
    // create a cost function
    auto costFunc =
        VisualOpticalFlowFactor<Configor::Prior::SplineOrder, deriv, IsInvDepth>::Create(
            so3Meta, scaleMeta, ofCorr, weight);

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

    // depth
    costFunc->AddParameterBlock(1);

    costFunc->SetNumResiduals(2);

    // organize the param block vector
    std::vector<double *> paramBlockVec;

    // so3 knots param block
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

    // lin acce knots
    AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                   scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

    auto SO3_DnToBr = parMagr->EXTRI.SO3_DnToBr.at(topic).data();
    paramBlockVec.push_back(SO3_DnToBr);

    auto POS_DnInBr = parMagr->EXTRI.POS_DnInBr.at(topic).data();
    paramBlockVec.push_back(POS_DnInBr);

    paramBlockVec.push_back(TO_DnToBr);
    paramBlockVec.push_back(RS_READOUT);

    paramBlockVec.push_back(intri->intri->FXAddress());
    paramBlockVec.push_back(intri->intri->FYAddress());
    paramBlockVec.push_back(intri->intri->CXAddress());
    paramBlockVec.push_back(intri->intri->CYAddress());

    if constexpr (IsInvDepth) {
        paramBlockVec.push_back(&ofCorr->invDepth);
    } else {
        paramBlockVec.push_back(&ofCorr->depth);
    }

    // pass to problem
    this->AddResidualBlock(costFunc,
                           new ceres::HuberLoss(Configor::Prior::LossForOpticalFlowFactor * weight),
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

    // two cases we do not estimate the depth:
    // 1. without 'Opt::OPT_VISUAL_DEPTH' option
    // 2. the rgbd camera moves too slow, and without depth observability
    if (!IsOptionWith(Opt::OPT_VISUAL_DEPTH, option) || !ofCorr->withDepthObservability) {
        if constexpr (IsInvDepth) {
            this->SetParameterBlockConstant(&ofCorr->invDepth);
        } else {
            this->SetParameterBlockConstant(&ofCorr->depth);
        }
    }
}

/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
 *   READOUT_TIME | FX | FY | CX | CY | DEPTH ]
 */
template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void Estimator::AddVisualOpticalFlowConstraint(const OpticalFlowCorr::Ptr &ofCorr,
                                               const std::string &topic,
                                               Opt option,
                                               double weight) {
    // invalid depth
    if (ofCorr->depth < 1E-3) {
        return;
    }

    auto &intri = parMagr->INTRI.Camera.at(topic);
    const double TO_PADDING = Configor::Prior::TimeOffsetPadding;
    const double RT_PADDING = Configor::Prior::ReadoutTimePadding;
    double *RS_READOUT = &parMagr->TEMPORAL.RS_READOUT.at(topic);
    double *TO_CmToBr = &parMagr->TEMPORAL.TO_CmToBr.at(topic);

    if (ofCorr->MidPointVel(*RS_READOUT).norm() < Configor::Prior::LossForOpticalFlowFactor) {
        // small pixel velocity
        return;
    }

    std::pair<double, double> timePair = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(1),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(1),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_CmToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_CmToBr, option)  // if opt time offset
    );

    if (!TimeInRangeForSplines(timePair)) {
        return;
    }

    // prepare metas for splines
    SplineMetaType so3Meta, scaleMeta;

    splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {timePair}, so3Meta);
    splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE, {timePair}, scaleMeta);

    static constexpr int deriv = TimeDeriv::Deriv<type, TimeDeriv::LIN_VEL>();
    // create a cost function
    auto costFunc =
        VisualOpticalFlowFactor<Configor::Prior::SplineOrder, deriv, IsInvDepth>::Create(
            so3Meta, scaleMeta, ofCorr, weight);

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

    costFunc->SetNumResiduals(2);

    // organize the param block vector
    std::vector<double *> paramBlockVec;

    // so3 knots param block
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

    // lin acce knots
    AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                   scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

    auto SO3_CmToBr = parMagr->EXTRI.SO3_CmToBr.at(topic).data();
    paramBlockVec.push_back(SO3_CmToBr);

    auto POS_CmInBr = parMagr->EXTRI.POS_CmInBr.at(topic).data();
    paramBlockVec.push_back(POS_CmInBr);

    paramBlockVec.push_back(TO_CmToBr);
    paramBlockVec.push_back(RS_READOUT);

    paramBlockVec.push_back(intri->FXAddress());
    paramBlockVec.push_back(intri->FYAddress());
    paramBlockVec.push_back(intri->CXAddress());
    paramBlockVec.push_back(intri->CYAddress());

    if constexpr (IsInvDepth) {
        paramBlockVec.push_back(&ofCorr->invDepth);
    } else {
        paramBlockVec.push_back(&ofCorr->depth);
    }

    // pass to problem
    this->AddResidualBlock(costFunc,
                           new ceres::HuberLoss(Configor::Prior::LossForOpticalFlowFactor * weight),
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

    if (!IsOptionWith(Opt::OPT_VISUAL_DEPTH, option)) {
        if constexpr (IsInvDepth) {
            this->SetParameterBlockConstant(&ofCorr->invDepth);
        } else {
            this->SetParameterBlockConstant(&ofCorr->depth);
        }
    }
}

/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
 * READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
 */
template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void Estimator::AddVisualOpticalFlowReprojConstraint(const OpticalFlowCorrPtr &ofCorr,
                                                     const std::string &topic,
                                                     Opt option,
                                                     double weight) {
    // invalid depth
    if (ofCorr->depth < 1E-3) {
        return;
    }

    auto &intri = parMagr->INTRI.Camera.at(topic);
    const double TO_PADDING = Configor::Prior::TimeOffsetPadding;
    const double RT_PADDING = Configor::Prior::ReadoutTimePadding;
    double *RS_READOUT = &parMagr->TEMPORAL.RS_READOUT.at(topic);
    double *TO_CmToBr = &parMagr->TEMPORAL.TO_CmToBr.at(topic);

    if (ofCorr->MidPointVel(*RS_READOUT).norm() < Configor::Prior::LossForOpticalFlowFactor) {
        // small pixel velocity
        return;
    }

    std::pair<double, double> timePairFir = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(0),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(0),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_CmToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_CmToBr, option)  // if opt time offset
    );
    std::pair<double, double> timePairMid = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(1),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(1),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_CmToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_CmToBr, option)  // if opt time offset
    );
    std::pair<double, double> timePairLast = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(2),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(2),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_CmToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_CmToBr, option)  // if opt time offset
    );
    if (!TimeInRangeForSplines(timePairFir) || !TimeInRangeForSplines(timePairMid) ||
        !TimeInRangeForSplines(timePairLast)) {
        return;
    }

    // prepare metas for splines
    SplineMetaType so3Meta, scaleMeta;

    splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                    {timePairFir, timePairMid, timePairLast}, so3Meta);
    splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE,
                                   {timePairFir, timePairMid, timePairLast}, scaleMeta);

    static constexpr int deriv = TimeDeriv::Deriv<type, TimeDeriv::LIN_POS>();
    // create a cost function
    auto costFunc =
        VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, deriv, IsInvDepth>::Create(
            so3Meta, scaleMeta, ofCorr, weight);

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

    costFunc->SetNumResiduals(4);

    // organize the param block vector
    std::vector<double *> paramBlockVec;

    // so3 knots param block
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

    // lin acce knots
    AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                   scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

    auto SO3_CmToBr = parMagr->EXTRI.SO3_CmToBr.at(topic).data();
    paramBlockVec.push_back(SO3_CmToBr);

    auto POS_CmInBr = parMagr->EXTRI.POS_CmInBr.at(topic).data();
    paramBlockVec.push_back(POS_CmInBr);

    paramBlockVec.push_back(TO_CmToBr);
    paramBlockVec.push_back(RS_READOUT);

    paramBlockVec.push_back(intri->FXAddress());
    paramBlockVec.push_back(intri->FYAddress());
    paramBlockVec.push_back(intri->CXAddress());
    paramBlockVec.push_back(intri->CYAddress());

    if constexpr (IsInvDepth) {
        paramBlockVec.push_back(&ofCorr->invDepth);
    } else {
        paramBlockVec.push_back(&ofCorr->depth);
    }

    // pass to problem
    this->AddResidualBlock(costFunc,
                           new ceres::HuberLoss(Configor::Prior::LossForReprojFactor * weight),
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

    if (!IsOptionWith(Opt::OPT_VISUAL_DEPTH, option)) {
        if constexpr (IsInvDepth) {
            this->SetParameterBlockConstant(&ofCorr->invDepth);
        } else {
            this->SetParameterBlockConstant(&ofCorr->depth);
        }
    }
}

/**
 * param blocks:
 * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_DnToBr | POS_DnInBr | TO_DnToBr |
 *   READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
 */
template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void Estimator::AddRGBDOpticalFlowReprojConstraint(const OpticalFlowCorrPtr &ofCorr,
                                                   const std::string &topic,
                                                   Opt option,
                                                   double weight) {
    // invalid depth
    if (ofCorr->depth < 1E-3) {
        return;
    }

    auto &intri = parMagr->INTRI.RGBD.at(topic)->intri;
    const double TO_PADDING = Configor::Prior::TimeOffsetPadding;
    const double RT_PADDING = Configor::Prior::ReadoutTimePadding;
    double *RS_READOUT = &parMagr->TEMPORAL.RS_READOUT.at(topic);
    double *TO_DnToBr = &parMagr->TEMPORAL.TO_DnToBr.at(topic);

    if (ofCorr->MidPointVel(*RS_READOUT).norm() < Configor::Prior::LossForOpticalFlowFactor) {
        // small pixel velocity
        return;
    }

    std::pair<double, double> timePairFir = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(0),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(0),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_DnToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_DnToBr, option)  // if opt time offset
    );
    std::pair<double, double> timePairMid = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(1),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(1),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_DnToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_DnToBr, option)  // if opt time offset
    );
    std::pair<double, double> timePairLast = ConsideredTimeRangeForCameraStamp(
        ofCorr->timeAry.at(2),                               // time stamped by the camera
        *RS_READOUT, RT_PADDING, ofCorr->rdFactorAry.at(2),  // the readout factor
        IsOptionWith(Opt::OPT_RS_CAM_READOUT_TIME, option),  // if optimize rs readout time
        *TO_DnToBr, TO_PADDING, IsOptionWith(Opt::OPT_TO_DnToBr, option)  // if opt time offset
    );

    if (!TimeInRangeForSplines(timePairFir) || !TimeInRangeForSplines(timePairMid) ||
        !TimeInRangeForSplines(timePairLast)) {
        return;
    }

    // prepare metas for splines
    SplineMetaType so3Meta, scaleMeta;

    splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                    {timePairFir, timePairMid, timePairLast}, so3Meta);
    splines->CalculateRdSplineMeta(Configor::Preference::SCALE_SPLINE,
                                   {timePairFir, timePairMid, timePairLast}, scaleMeta);

    static constexpr int deriv = TimeDeriv::Deriv<type, TimeDeriv::LIN_POS>();
    // create a cost function
    auto costFunc =
        VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, deriv, IsInvDepth>::Create(
            so3Meta, scaleMeta, ofCorr, weight);

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

    costFunc->SetNumResiduals(4);

    // organize the param block vector
    std::vector<double *> paramBlockVec;

    // so3 knots param block
    AddSo3KnotsData(paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                    !IsOptionWith(Opt::OPT_SO3_SPLINE, option));

    // lin acce knots
    AddRdKnotsData(paramBlockVec, splines->GetRdSpline(Configor::Preference::SCALE_SPLINE),
                   scaleMeta, !IsOptionWith(Opt::OPT_SCALE_SPLINE, option));

    auto SO3_DnToBr = parMagr->EXTRI.SO3_DnToBr.at(topic).data();
    paramBlockVec.push_back(SO3_DnToBr);

    auto POS_DnInBr = parMagr->EXTRI.POS_DnInBr.at(topic).data();
    paramBlockVec.push_back(POS_DnInBr);

    paramBlockVec.push_back(TO_DnToBr);
    paramBlockVec.push_back(RS_READOUT);

    paramBlockVec.push_back(intri->FXAddress());
    paramBlockVec.push_back(intri->FYAddress());
    paramBlockVec.push_back(intri->CXAddress());
    paramBlockVec.push_back(intri->CYAddress());

    if constexpr (IsInvDepth) {
        paramBlockVec.push_back(&ofCorr->invDepth);
    } else {
        paramBlockVec.push_back(&ofCorr->depth);
    }

    // pass to problem
    this->AddResidualBlock(costFunc,
                           new ceres::HuberLoss(Configor::Prior::LossForReprojFactor * weight),
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
        this->SetParameterBlockConstant(intri->FXAddress());
        this->SetParameterBlockConstant(intri->FYAddress());
    }

    if (!IsOptionWith(Opt::OPT_CAM_PRINCIPAL_POINT, option)) {
        this->SetParameterBlockConstant(intri->CXAddress());
        this->SetParameterBlockConstant(intri->CYAddress());
    }

    if (!IsOptionWith(Opt::OPT_VISUAL_DEPTH, option)) {
        if constexpr (IsInvDepth) {
            this->SetParameterBlockConstant(&ofCorr->invDepth);
        } else {
            this->SetParameterBlockConstant(&ofCorr->depth);
        }
    }
}
}  // namespace ns_ikalibr

#endif  // IKALIBR_ESTIMATOR_TPL_HPP
