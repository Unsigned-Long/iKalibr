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

#include "solver/calib_solver_tpl.hpp"
#include "core/lidar_odometer.h"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void CalibSolver::InitScaleSpline() const {
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);

    spdlog::info("performing scale spline recovery...");

    auto estimator = Estimator::Create(_splines, _parMagr);
    auto optOption = OptOption::OPT_SCALE_SPLINE;

    switch (GetScaleType()) {
        case TimeDeriv::LIN_ACCE_SPLINE: {
            // only multiple imus are involved
            this->AddAcceFactor<TimeDeriv::LIN_ACCE_SPLINE>(
                estimator, Configor::DataStream::ReferIMU, optOption);
        } break;
        case TimeDeriv::LIN_VEL_SPLINE: {
            // only multiple radars and imus are involved
            constexpr int VelDeriv =
                TimeDeriv::Deriv<TimeDeriv::LIN_VEL_SPLINE, TimeDeriv::LIN_VEL>();
            optOption |= OptOption::OPT_SO3_SPLINE;

            // radar-derived velocities
            for (const auto &[topic, data] : _dataMagr->GetRadarMeasurements()) {
                const double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(topic);
                const auto &SO3_RjToBr = _parMagr->EXTRI.SO3_RjToBr.at(topic);
                const Eigen::Vector3d &POS_RjInBr = _parMagr->EXTRI.POS_RjInBr.at(topic);
                // const double weight = Configor::DataStream::RadarTopics.at(topic).Weight;
                constexpr double weight = 10.0;

                for (const auto &ary : data) {
                    double timeByBr = ary->GetTimestamp() + TO_RjToBr;
                    if (ary->GetTargets().size() < 10 || !so3Spline.TimeStampInRange(timeByBr)) {
                        continue;
                    }

                    auto SO3_BrToW = so3Spline.Evaluate(timeByBr);
                    auto ANG_VEL_BrToWInW = SO3_BrToW * so3Spline.VelocityBody(timeByBr);
                    Eigen::Vector3d LIN_VEL_BrToWInW =
                        SO3_BrToW * SO3_RjToBr * ary->RadarVelocityFromStaticTargetArray() -
                        Sophus::SO3d::hat(ANG_VEL_BrToWInW) * (SO3_BrToW * POS_RjInBr);

                    estimator->AddLinearScaleConstraint<VelDeriv>(
                        timeByBr,          // time stamped the reference imu
                        LIN_VEL_BrToWInW,  // the linear velocity
                        optOption,         // the optimization option
                        weight);           // the weigh
                }
            }

            // rgbd-derived velocities
            for (const auto &[rgbdTopic, bodyFrameVels] : _initAsset->rgbdBodyFrameVels) {
                // const double weight = Configor::DataStream::RGBDTopics.at(rgbdTopic).Weight;
                constexpr double weight = 10.0;
                double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(rgbdTopic);
                const auto &SO3_DnToBr = _parMagr->EXTRI.SO3_DnToBr.at(rgbdTopic);
                const Eigen::Vector3d &POS_DnInBr = _parMagr->EXTRI.POS_DnInBr.at(rgbdTopic);

                for (const auto &[frame, vel] : bodyFrameVels) {
                    double timeByBr = frame->GetTimestamp() + TO_DnToBr;
                    if (!so3Spline.TimeStampInRange(timeByBr)) {
                        continue;
                    }

                    auto SO3_BrToW = so3Spline.Evaluate(timeByBr);
                    auto ANG_VEL_BrToWInW = SO3_BrToW * so3Spline.VelocityBody(timeByBr);
                    Eigen::Vector3d LIN_VEL_BrToWInW =
                        SO3_BrToW * SO3_DnToBr * vel -
                        Sophus::SO3d::hat(ANG_VEL_BrToWInW) * (SO3_BrToW * POS_DnInBr);

                    estimator->AddLinearScaleConstraint<VelDeriv>(
                        timeByBr,          // time stamped the reference imu
                        LIN_VEL_BrToWInW,  // the linear velocity
                        optOption,         // the optimization option
                        weight);           // the weigh
                }
            }

            // camera-derived velocities
            for (const auto &[camTopic, bodyFrameVels] : _initAsset->velCamBodyFrameVels) {
                constexpr double weight = 10.0;
                double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                const auto &SO3_CmToBr = _parMagr->EXTRI.SO3_CmToBr.at(camTopic);
                const Eigen::Vector3d &POS_CmInBr = _parMagr->EXTRI.POS_CmInBr.at(camTopic);

                for (const auto &[frame, vel] : bodyFrameVels) {
                    double timeByBr = frame->GetTimestamp() + TO_CmToBr;
                    if (!so3Spline.TimeStampInRange(timeByBr)) {
                        continue;
                    }

                    auto SO3_BrToW = so3Spline.Evaluate(timeByBr);
                    auto ANG_VEL_BrToWInW = SO3_BrToW * so3Spline.VelocityBody(timeByBr);
                    Eigen::Vector3d LIN_VEL_BrToWInW =
                        SO3_BrToW * SO3_CmToBr * vel -
                        Sophus::SO3d::hat(ANG_VEL_BrToWInW) * (SO3_BrToW * POS_CmInBr);

                    estimator->AddLinearScaleConstraint<VelDeriv>(
                        timeByBr,          // time stamped the reference imu
                        LIN_VEL_BrToWInW,  // the linear velocity
                        optOption,         // the optimization option
                        weight);           // the weigh
                }
            }

            // add acceleration factor using inertial measurements only from reference IMU
            this->AddAcceFactor<TimeDeriv::LIN_VEL_SPLINE>(
                estimator, Configor::DataStream::ReferIMU, optOption);
            this->AddGyroFactor(estimator, Configor::DataStream::ReferIMU, optOption);

            /**
             * if optimize time offsets, we first recover the scale spline, then construct a ls
             * problem to recover time offsets, only for radars
             */
            if (Configor::IsRadarIntegrated() && Configor::Prior::OptTemporalParams) {
                // we don't want to output the solving information
                auto solveOpt = Estimator::DefaultSolverOptions(
                    // the thread count
                    Configor::Preference::AvailableThreads(),
                    // do not output information
                    false, Configor::Preference::UseCudaInSolving);
                estimator->Solve(solveOpt);

                /**
                 * in the new constructed ls problem, we estimate the time offsets of radars
                 */
                estimator = Estimator::Create(_splines, _parMagr);
                optOption = OptOption::OPT_TO_RjToBr;
                for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
                    this->AddRadarFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
                }
            }

        } break;
        case TimeDeriv::LIN_POS_SPLINE: {
            // this value equals to zero
            constexpr int PosDeriv =
                TimeDeriv::Deriv<TimeDeriv::LIN_POS_SPLINE, TimeDeriv::LIN_POS>();
            optOption |= OptOption::OPT_SO3_SPLINE;

            /**
             * we do not directly operate the origin spline, we construct new splines with large
             * knot distance and recover them use the lidar/camera-derived poses, then based on
             * these rough spline, we recover the origin splines
             */
            std::vector<double> headTimeVec, tailTimeVec;
            for (const auto &[lidarTopic, odometer] : _initAsset->lidarOdometers) {
                const double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
                headTimeVec.push_back(odometer->GetOdomPoseVec().front().timeStamp + TO_LkToBr);
                tailTimeVec.push_back(odometer->GetOdomPoseVec().back().timeStamp + TO_LkToBr);
            }
            for (const auto &[camTopic, poseSeq] : _initAsset->sfmPoseSeq) {
                const double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                headTimeVec.push_back(poseSeq.front().timeStamp + TO_CmToBr);
                tailTimeVec.push_back(poseSeq.back().timeStamp + TO_CmToBr);
            }

            const double minTime = _dataMagr->GetCalibStartTimestamp();
            const double maxTime = _dataMagr->GetCalibEndTimestamp();

            // create rough splines (the time distance is larger than that from configure as our
            // poses are not dense)
            double dtRoughSum = 0.0;
            int count = 0;
            if (Configor::IsLiDARIntegrated()) {
                dtRoughSum += 2.0 / _dataMagr->GetLiDARAvgFrequency();
                ++count;
            }
            if (Configor::IsPosCameraIntegrated()) {
                dtRoughSum += 2.0 / _dataMagr->GetCameraAvgFrequency();
                ++count;
            }
            const double dtRough = dtRoughSum / count;
            // create the rough splines
            auto roughSplines = CreateSplineBundle(minTime, maxTime, dtRough, dtRough);
            const auto &rSo3Spline = roughSplines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
            const auto &rScaleSpline =
                roughSplines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

            estimator = Estimator::Create(roughSplines, _parMagr);

            for (const auto &[lidarTopic, odometer] : _initAsset->lidarOdometers) {
                double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
                // const double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;
                constexpr double weight = 10.0;
                auto SE3_BrToLk = _parMagr->EXTRI.SE3_LkToBr(lidarTopic).inverse();

                // the translation part is not exactly rigorous as the translation spline is not
                // recovered yet this assumes that the origin of map frame of {Lk} is the same as
                // that of the world frame {br0}
                auto SE3_Lk0ToBr0 = this->CurLkToW(odometer->GetMapTime(), lidarTopic);

                if (SE3_Lk0ToBr0 == std::nullopt) {
                    throw Status(Status::CRITICAL,
                                 "map time of '{}' is out of time range of splines!", lidarTopic);
                }

                for (const auto &SE3_LkToLk0 : odometer->GetOdomPoseVec()) {
                    auto SE3_BrToBr0 = *SE3_Lk0ToBr0 * SE3_LkToLk0.se3() * SE3_BrToLk;
                    estimator->AddSO3Constraint(
                        SE3_LkToLk0.timeStamp + TO_LkToBr,  // the time stamped the reference imu
                        SE3_BrToBr0.so3(),                  // the rotation
                        optOption,                          // the optimization option
                        weight);                            // the weight

                    estimator->AddLinearScaleConstraint<PosDeriv>(
                        SE3_LkToLk0.timeStamp + TO_LkToBr,  // the time stamped the reference imu
                        SE3_BrToBr0.translation(),          // the translation
                        optOption,                          // the optimization option
                        weight);                            // the weight
                }
            }

            for (const auto &[camTopic, poseSeq] : _initAsset->sfmPoseSeq) {
                double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                // const double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;
                constexpr double weight = 10.0;
                auto SE3_BrToCm = _parMagr->EXTRI.SE3_CmToBr(camTopic).inverse();

                /**
                 * the translation part is not exactly rigorous as the translation spline is not
                 * recovered yet this assumes that the origin of map frame of {Cm} is the same as
                 * that of the world frame {br0}
                 */
                auto SE3_Cm0ToBr0 = this->CurCmToW(poseSeq.front().timeStamp, camTopic);

                if (SE3_Cm0ToBr0 == std::nullopt) {
                    throw Status(Status::CRITICAL,
                                 "map time of '{}' is out of time range of splines!", camTopic);
                }

                for (const auto &SE3_CmToCm0 : poseSeq) {
                    auto SE3_BrToBr0 = *SE3_Cm0ToBr0 * SE3_CmToCm0.se3() * SE3_BrToCm;
                    estimator->AddSO3Constraint(
                        SE3_CmToCm0.timeStamp + TO_CmToBr,  // the time stamped the reference imu
                        SE3_BrToBr0.so3(),                  // the rotation
                        optOption,                          // the optimization option
                        weight);                            // the weight
                    estimator->AddLinearScaleConstraint<PosDeriv>(
                        SE3_CmToCm0.timeStamp + TO_CmToBr,  // the time stamped the reference imu
                        SE3_BrToBr0.translation(),          // the translation
                        optOption,                          // the optimization option
                        weight);                            // the weight
                }
            }

            // add tail factors (constraints) to maintain enough observability
            estimator->AddLinScaleTailConstraint(optOption, 1.0);
            estimator->AddSO3TailConstraint(optOption, 1.0);

            // add inertial factor using inertial measurements only from reference IMU
            this->AddAcceFactor<TimeDeriv::LIN_POS_SPLINE>(
                estimator, Configor::DataStream::ReferIMU, optOption);
            this->AddGyroFactor(estimator, Configor::DataStream::ReferIMU, optOption);

            spdlog::info("fitting rough splines using sensor-derived pose sequence...");
            // estimator->PrintUninvolvedKnots();

            // we don't want to output the solving information
            auto solveOpt = Estimator::DefaultSolverOptions(
                // the thread count
                Configor::Preference::AvailableThreads(),
                // do not output information
                false, Configor::Preference::UseCudaInSolving);
            estimator->Solve(solveOpt, _priori);

            spdlog::info("fitting rough splines finished.");

            estimator = Estimator::Create(_splines, _parMagr);
            for (double t = minTime; t < maxTime;) {
                estimator->AddSO3Constraint(t,  // the time stamped the reference imu
                                            rSo3Spline.Evaluate(t),  // the rotation
                                            optOption,               // the optimization option
                                            1.0);                    // the weight
                estimator->AddLinearScaleConstraint<PosDeriv>(
                    t,                         // the time stamped the reference imu
                    rScaleSpline.Evaluate(t),  // the translation
                    optOption,                 // the optimization option
                    1.0);                      // the weight
                t += 0.01;
            }
            // add tail factors (constraints) to maintain enough observability
            estimator->AddLinScaleTailConstraint(optOption, 1.0);
            estimator->AddSO3TailConstraint(optOption, 1.0);
            // estimator->PrintUninvolvedKnots();
        } break;
    }
    auto sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

    if (GetScaleType() == TimeDeriv::LIN_POS_SPLINE && Configor::IsRadarIntegrated() &&
        Configor::Prior::OptTemporalParams) {
        // in this case, the time offsets of radars have not been recovered
        estimator = Estimator::Create(_splines, _parMagr);
        for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
            this->AddRadarFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic,
                                                            OptOption::OPT_TO_RjToBr);
        }
        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }
}
}  // namespace ns_ikalibr