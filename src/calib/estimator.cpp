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

#include "calib/estimator.h"
#include "factor/imu_gyro_factor.hpp"
#include "factor/inertial_align_factor.hpp"
#include "factor/lidar_inertial_align_factor.hpp"
#include "factor/visual_inertial_align_factor.hpp"
#include "factor/radar_inertial_align_factor.hpp"
#include "factor/hand_eye_rot_align_factor.hpp"
#include "factor/so3_factor.hpp"
#include "factor/visual_proj_factor.hpp"
#include "factor/linear_knots_factor.hpp"
#include "ctraj/core/trajectory_estimator.h"

_3_

namespace ns_ikalibr {

    std::shared_ptr<ceres::EigenQuaternionManifold> Estimator::QUATER_MANIFOLD(new ceres::EigenQuaternionManifold());
    std::shared_ptr<ceres::SphereManifold<3>> Estimator::GRAVITY_MANIFOLD(new ceres::SphereManifold<3>());

    ceres::Problem::Options Estimator::DefaultProblemOptions() {
        return ns_ctraj::TrajectoryEstimator<Configor::Prior::SplineOrder>::DefaultProblemOptions();
    }

    ceres::Solver::Options Estimator::DefaultSolverOptions(int threadNum, bool toStdout, bool useCUDA) {
        auto defaultSolverOptions = ns_ctraj::TrajectoryEstimator<Configor::Prior::SplineOrder>::DefaultSolverOptions(
                threadNum, toStdout, useCUDA
        );
        if (!useCUDA) {
            defaultSolverOptions.linear_solver_type = ceres::DENSE_SCHUR;
        }
        defaultSolverOptions.trust_region_strategy_type = ceres::DOGLEG;
        return defaultSolverOptions;
    }

    Estimator::Estimator(SplineBundleType::Ptr splines, CalibParamManager::Ptr calibParamManager)
            : ceres::Problem(DefaultProblemOptions()), splines(std::move(splines)),
              parMagr(std::move(calibParamManager)) {}

    Estimator::Ptr
    Estimator::Create(const SplineBundleType::Ptr &splines, const CalibParamManager::Ptr &calibParamManager) {
        return std::make_shared<Estimator>(splines, calibParamManager);
    }

    ceres::Solver::Summary Estimator::Solve(const ceres::Solver::Options &options) {
        ceres::Solver::Summary summary;
        ceres::Solve(options, this, &summary);
        return summary;
    }

    void Estimator::AddRdKnotsData(std::vector<double *> &paramBlockVec,
                                   const Estimator::SplineBundleType::RdSplineType &spline,
                                   const Estimator::SplineMetaType &splineMeta, bool setToConst) {
        // for each segment
        for (const auto &seg: splineMeta.segments) {
            // the factor 'seg.dt * 0.5' is the treatment for numerical accuracy
            auto idxMaster = spline.ComputeTIndex(seg.t0 + seg.dt * 0.5).second;

            // from the first control point to the last control point
            for (std::size_t i = idxMaster; i < idxMaster + seg.NumParameters(); ++i) {
                auto *data = const_cast<double *>(spline.GetKnot(static_cast<int>(i)).data());

                this->AddParameterBlock(data, 3);

                paramBlockVec.push_back(data);
                // set this param block to be constant
                if (setToConst) { this->SetParameterBlockConstant(data); }
            }
        }
    }

    void Estimator::AddSo3KnotsData(std::vector<double *> &paramBlockVec,
                                    const Estimator::SplineBundleType::So3SplineType &spline,
                                    const Estimator::SplineMetaType &splineMeta, bool setToConst) {
        // for each segment
        for (const auto &seg: splineMeta.segments) {
            // the factor 'seg.dt * 0.5' is the treatment for numerical accuracy
            auto idxMaster = spline.ComputeTIndex(seg.t0 + seg.dt * 0.5).second;

            // from the first control point to the last control point
            for (std::size_t i = idxMaster; i < idxMaster + seg.NumParameters(); ++i) {
                auto *data = const_cast<double *>(spline.GetKnot(static_cast<int>(i)).data());
                // the local parameterization is very important!!!
                this->AddParameterBlock(data, 4, QUATER_MANIFOLD.get());

                paramBlockVec.push_back(data);
                // set this param block to be constant
                if (setToConst) { this->SetParameterBlockConstant(data); }
            }
        }
    }

    Eigen::MatrixXd Estimator::CRSMatrix2EigenMatrix(ceres::CRSMatrix *jacobian_crs_matrix) {
        Eigen::MatrixXd J(jacobian_crs_matrix->num_rows, jacobian_crs_matrix->num_cols);
        J.setZero();

        std::vector<int> jacobian_crs_matrix_rows, jacobian_crs_matrix_cols;
        std::vector<double> jacobian_crs_matrix_values;
        jacobian_crs_matrix_rows = jacobian_crs_matrix->rows;
        jacobian_crs_matrix_cols = jacobian_crs_matrix->cols;
        jacobian_crs_matrix_values = jacobian_crs_matrix->values;

        int cur_index_in_cols_and_values = 0;
        // rows is a num_rows + 1 sized array
        int row_size = static_cast<int>(jacobian_crs_matrix_rows.size()) - 1;
        // outer loop traverse rows, inner loop traverse cols and values
        for (int row_index = 0; row_index < row_size; ++row_index) {
            while (cur_index_in_cols_and_values < jacobian_crs_matrix_rows[row_index + 1]) {
                J(row_index, jacobian_crs_matrix_cols[cur_index_in_cols_and_values]) =
                        jacobian_crs_matrix_values[cur_index_in_cols_and_values];
                cur_index_in_cols_and_values++;
            }
        }
        return J;
    }

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | GYRO_BIAS | GYRO_MAP_COEFF | SO3_AtoG | SO3_BiToBr | TO_BiToBr ]
     */
    void Estimator::AddIMUGyroMeasurement(const IMUFrame::Ptr &imuFrame, const std::string &topic,
                                          Opt option, double gyroWeight) {
        // prepare metas for splines
        SplineMetaType so3Meta;

        // different relative control points finding [single vs. range]
        // for the inertial measurements from the reference IMU, there is no need to consider a time padding,
        // as its time offsets would be fixed as identity
        if (IsOptionWith(Opt::OPT_TO_BiToBr, option) && Configor::DataStream::ReferIMU != topic) {
            double minTime = imuFrame->GetTimestamp() - Configor::Prior::TimeOffsetPadding;
            double maxTime = imuFrame->GetTimestamp() + Configor::Prior::TimeOffsetPadding;
            // invalid time stamp
            if (!splines->TimeInRangeForSo3(minTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(maxTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{minTime, maxTime}}, so3Meta);
        } else {
            double curTime = imuFrame->GetTimestamp() + parMagr->TEMPORAL.TO_BiToBr.at(topic);

            // check point time stamp
            if (!splines->TimeInRangeForSo3(curTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{curTime, curTime}}, so3Meta);
        }

        // create a cost function
        auto costFunc = IMUGyroFactor<Configor::Prior::SplineOrder>::Create(so3Meta, imuFrame, gyroWeight);

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }

        // GYRO gyroBias
        costFunc->AddParameterBlock(3);
        // GYRO map coeff
        costFunc->AddParameterBlock(6);
        // SO3_AtoG
        costFunc->AddParameterBlock(4);
        // SO3_BiToBr
        costFunc->AddParameterBlock(4);
        // TIME_OFFSET_BiToBc
        costFunc->AddParameterBlock(1);

        // set Residuals
        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(
                paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                !IsOptionWith(Opt::OPT_SO3_SPLINE, option)
        );

        // GYRO gyroBias
        auto gyroBias = parMagr->INTRI.IMU.at(topic)->GYRO.BIAS.data();
        paramBlockVec.push_back(gyroBias);
        // GYRO map coeff
        auto gyroMapCoeff = parMagr->INTRI.IMU.at(topic)->GYRO.MAP_COEFF.data();
        paramBlockVec.push_back(gyroMapCoeff);
        // SO3_AtoG
        auto SO3_AtoG = parMagr->INTRI.IMU.at(topic)->SO3_AtoG.data();
        paramBlockVec.push_back(SO3_AtoG);
        // SO3_BiToBr
        auto SO3_BiToBr = parMagr->EXTRI.SO3_BiToBr.at(topic).data();
        paramBlockVec.push_back(SO3_BiToBr);
        // TIME_OFFSET_BiToBc
        auto TIME_OFFSET_BiToBc = &parMagr->TEMPORAL.TO_BiToBr.at(topic);
        paramBlockVec.push_back(TIME_OFFSET_BiToBc);

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

        this->SetManifold(SO3_AtoG, QUATER_MANIFOLD.get());
        this->SetManifold(SO3_BiToBr, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_GYRO_BIAS, option)) {
            this->SetParameterBlockConstant(gyroBias);
        }

        if (!IsOptionWith(Opt::OPT_GYRO_MAP_COEFF, option)) {
            this->SetParameterBlockConstant(gyroMapCoeff);
        }

        if (!IsOptionWith(Opt::OPT_SO3_AtoG, option)) {
            this->SetParameterBlockConstant(SO3_AtoG);
        }

        if (!IsOptionWith(Opt::OPT_SO3_BiToBr, option)) {
            this->SetParameterBlockConstant(SO3_BiToBr);
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
     * [ SO3_LkToBr | POS_LkInBr | POS_BiInBr | S_VEL | E_VEL | GRAVITY ]
     */
    void Estimator::AddLiDARInertialAlignment(const std::vector<IMUFrame::Ptr> &data,
                                              const std::string &lidarTopic, const std::string &imuTopic,
                                              const ns_ctraj::Posed &sPose, const ns_ctraj::Posed &ePose,
                                              double mapTime, Eigen::Vector3d *sVel, Eigen::Vector3d *eVel,
                                              Estimator::Opt option, double weight) {
        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        double st = sPose.timeStamp, et = ePose.timeStamp, TO_LkToBr = parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);

        if (!so3Spline.TimeStampInRange(st + TO_LkToBr) || !so3Spline.TimeStampInRange(et + TO_LkToBr)) { return; }

        if (!so3Spline.TimeStampInRange(mapTime + TO_LkToBr)) {
            throw Status(Status::CRITICAL, "the map time is not in spline time range!");
        }

        double TO_LkToBi = TO_LkToBr - parMagr->TEMPORAL.TO_BiToBr.at(imuTopic);
        auto [velVecMat, posVecMat] = InertialPosIntegration(data, imuTopic, st + TO_LkToBi, et + TO_LkToBi);

        // create a cost function
        auto helper = LiDARInertialAlignHelper<Configor::Prior::SplineOrder>(
                so3Spline, sPose, ePose, mapTime, TO_LkToBr, velVecMat, posVecMat
        );
        auto costFunc = LiDARInertialAlignFactor<Configor::Prior::SplineOrder>::Create(helper, weight);

        costFunc->AddParameterBlock(4);
        costFunc->AddParameterBlock(3);

        costFunc->AddParameterBlock(3);

        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(3);

        costFunc->AddParameterBlock(3);

        costFunc->SetNumResiduals(6);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        auto SO3_LkToBr = parMagr->EXTRI.SO3_LkToBr.at(lidarTopic).data();
        paramBlockVec.push_back(SO3_LkToBr);

        auto POS_LkInBr = parMagr->EXTRI.POS_LkInBr.at(lidarTopic).data();
        paramBlockVec.push_back(POS_LkInBr);

        auto POS_BiInBr = parMagr->EXTRI.POS_BiInBr.at(imuTopic).data();
        paramBlockVec.push_back(POS_BiInBr);

        paramBlockVec.push_back(sVel->data());
        paramBlockVec.push_back(eVel->data());

        auto GRAVITY = parMagr->GRAVITY.data();
        paramBlockVec.push_back(GRAVITY);

        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);
        this->SetManifold(SO3_LkToBr, QUATER_MANIFOLD.get());
        this->SetManifold(GRAVITY, GRAVITY_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_SO3_LkToBr, option)) {
            this->SetParameterBlockConstant(SO3_LkToBr);
        }
        if (!IsOptionWith(Opt::OPT_POS_LkInBr, option)) {
            this->SetParameterBlockConstant(POS_LkInBr);
        }
        if (!IsOptionWith(Opt::OPT_POS_BiInBr, option)) {
            this->SetParameterBlockConstant(POS_BiInBr);
        }
        if (!IsOptionWith(Opt::OPT_GRAVITY, option)) {
            this->SetParameterBlockConstant(GRAVITY);
        }
    }

    /**
     * param blocks:
     * [ SO3_CmToBr | POS_CmInBr | POS_BiInBr | S_VEL | E_VEL | GRAVITY | SCALE ]
     */
    void Estimator::AddVisualInertialAlignment(const std::vector<IMUFrame::Ptr> &data, const std::string &camTopic,
                                               const std::string &imuTopic, const ns_ctraj::Posed &sPose,
                                               const ns_ctraj::Posed &ePose, double mapTime, Eigen::Vector3d *sVel,
                                               Eigen::Vector3d *eVel, double *SCALE, Estimator::Opt option,
                                               double weight) {
        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        double st = sPose.timeStamp, et = ePose.timeStamp, TO_CmToBr = parMagr->TEMPORAL.TO_CmToBr.at(camTopic);

        if (!so3Spline.TimeStampInRange(st + TO_CmToBr) || !so3Spline.TimeStampInRange(et + TO_CmToBr)) { return; }

        if (!so3Spline.TimeStampInRange(mapTime + TO_CmToBr)) {
            throw Status(Status::CRITICAL, "the map time is not in spline time range!");
        }

        double TO_CmToBi = TO_CmToBr - parMagr->TEMPORAL.TO_BiToBr.at(imuTopic);
        auto [velVecMat, posVecMat] = InertialPosIntegration(data, imuTopic, st + TO_CmToBi, et + TO_CmToBi);

        // create a cost function
        auto helper = VisualInertialAlignHelper<Configor::Prior::SplineOrder>(
                so3Spline, sPose, ePose, mapTime, TO_CmToBr, velVecMat, posVecMat
        );
        auto costFunc = VisualInertialAlignFactor<Configor::Prior::SplineOrder>::Create(helper, weight);

        costFunc->AddParameterBlock(4);
        costFunc->AddParameterBlock(3);

        costFunc->AddParameterBlock(3);

        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(3);

        costFunc->AddParameterBlock(3);

        costFunc->AddParameterBlock(1);

        costFunc->SetNumResiduals(6);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        auto SO3_CmToBr = parMagr->EXTRI.SO3_CmToBr.at(camTopic).data();
        paramBlockVec.push_back(SO3_CmToBr);

        auto POS_CmInBr = parMagr->EXTRI.POS_CmInBr.at(camTopic).data();
        paramBlockVec.push_back(POS_CmInBr);

        auto POS_BiInBr = parMagr->EXTRI.POS_BiInBr.at(imuTopic).data();
        paramBlockVec.push_back(POS_BiInBr);

        paramBlockVec.push_back(sVel->data());
        paramBlockVec.push_back(eVel->data());

        auto GRAVITY = parMagr->GRAVITY.data();
        paramBlockVec.push_back(GRAVITY);

        paramBlockVec.push_back(SCALE);

        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);
        this->SetManifold(SO3_CmToBr, QUATER_MANIFOLD.get());
        this->SetManifold(GRAVITY, GRAVITY_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_SO3_CmToBr, option)) {
            this->SetParameterBlockConstant(SO3_CmToBr);
        }
        if (!IsOptionWith(Opt::OPT_POS_CmInBr, option)) {
            this->SetParameterBlockConstant(POS_CmInBr);
        }
        if (!IsOptionWith(Opt::OPT_POS_BiInBr, option)) {
            this->SetParameterBlockConstant(POS_BiInBr);
        }
        if (!IsOptionWith(Opt::OPT_GRAVITY, option)) {
            this->SetParameterBlockConstant(GRAVITY);
        }
        if (!IsOptionWith(Opt::OPT_VISUAL_GLOBAL_SCALE, option)) {
            this->SetParameterBlockConstant(SCALE);
        } else {
            this->SetParameterLowerBound(SCALE, 0, 1E-3);
        }
    }

    /**
     * param blocks:
     * [ SO3_CmToBr | POS_CmInBr | POS_BiInBr | S_VEL | E_VEL | GRAVITY | SCALE ]
     */
    /**
     * param blocks:
     * [ POS_BiInBr | START_VEL | END_VEL | GRAVITY ]
     */
    void Estimator::AddInertialAlignment(const std::vector<IMUFrame::Ptr> &data, const std::string &imuTopic,
                                         double sTimeByBr, double eTimeByBr, Eigen::Vector3d *sVel,
                                         Eigen::Vector3d *eVel, Estimator::Opt option, double weight) {
        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        if (!so3Spline.TimeStampInRange(sTimeByBr) || !so3Spline.TimeStampInRange(eTimeByBr)) { return; }

        double TO_BrToBi = -parMagr->TEMPORAL.TO_BiToBr.at(imuTopic);
        auto velVecMat = InertialVelIntegration(data, imuTopic, sTimeByBr + TO_BrToBi, eTimeByBr + TO_BrToBi);

        auto helper = InertialAlignHelper(eTimeByBr - sTimeByBr, velVecMat);
        auto costFunc = InertialAlignFactor::Create(helper, weight);

        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(3);

        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        auto POS_BiInBr = parMagr->EXTRI.POS_BiInBr.at(imuTopic).data();
        paramBlockVec.push_back(POS_BiInBr);

        auto GRAVITY = parMagr->GRAVITY.data();
        paramBlockVec.push_back(GRAVITY);

        paramBlockVec.push_back(sVel->data());
        paramBlockVec.push_back(eVel->data());

        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);
        this->SetManifold(GRAVITY, GRAVITY_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_POS_BiInBr, option)) {
            this->SetParameterBlockConstant(POS_BiInBr);
        }
        if (!IsOptionWith(Opt::OPT_GRAVITY, option)) {
            this->SetParameterBlockConstant(GRAVITY);
        }
    }

    /**
     * param blocks:
     * [ POS_BiInBr | SO3_RjToBr | POS_RjInBr | GRAVITY ]
     */
    void Estimator::AddRadarInertialAlignment(const std::vector<IMUFrame::Ptr> &data, const std::string &imuTopic,
                                              const std::string &radarTopic, const RadarTargetArray::Ptr &sRadarAry,
                                              const RadarTargetArray::Ptr &eRadarAry, Estimator::Opt option,
                                              double weight) {
        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        double st = sRadarAry->GetTimestamp(), et = eRadarAry->GetTimestamp();
        double TO_RjToBr = parMagr->TEMPORAL.TO_RjToBr.at(radarTopic);

        if (!so3Spline.TimeStampInRange(st + TO_RjToBr) || !so3Spline.TimeStampInRange(et + TO_RjToBr)) { return; }

        double TO_RjToBi = TO_RjToBr - parMagr->TEMPORAL.TO_BiToBr.at(imuTopic);
        auto velVecMat = InertialVelIntegration(data, imuTopic, st + TO_RjToBi, et + TO_RjToBi);

        // create a cost function
        auto helper = RadarInertialAlignHelper<Configor::Prior::SplineOrder>(
                so3Spline, sRadarAry, eRadarAry, TO_RjToBr, velVecMat
        );
        auto costFunc = RadarInertialAlignFactor<Configor::Prior::SplineOrder>::Create(helper, weight);

        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(4);
        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(3);

        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // POS_BiInBc
        auto POS_BiInBr = parMagr->EXTRI.POS_BiInBr.at(imuTopic).data();
        paramBlockVec.push_back(POS_BiInBr);
        // SO3_RjToBc
        auto SO3_RjToBr = parMagr->EXTRI.SO3_RjToBr.at(radarTopic).data();
        paramBlockVec.push_back(SO3_RjToBr);
        // POS_RjInBc
        auto POS_RjInBr = parMagr->EXTRI.POS_RjInBr.at(radarTopic).data();
        paramBlockVec.push_back(POS_RjInBr);
        // GRAVITY
        auto gravity = parMagr->GRAVITY.data();

        paramBlockVec.push_back(gravity);

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);
        this->SetManifold(gravity, GRAVITY_MANIFOLD.get());
        this->SetManifold(SO3_RjToBr, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_GRAVITY, option)) {
            this->SetParameterBlockConstant(gravity);
        }
        if (!IsOptionWith(Opt::OPT_POS_BiInBr, option)) {
            this->SetParameterBlockConstant(POS_BiInBr);
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
     * [ SO3 | ... | SO3 | SO3_LkToBr | TO_LkToBr ]
     */
    void Estimator::AddHandEyeRotationAlignmentForLiDAR(const std::string &lidarTopic, double tLastByLk,
                                                        double tCurByLk, const Sophus::SO3d &so3LastLkToM,
                                                        const Sophus::SO3d &so3CurLkToM, Estimator::Opt option,
                                                        double weight) {
        // prepare metas for splines
        SplineMetaType so3Meta;

        // different relative control points finding [single vs. range]
        if (IsOptionWith(Opt::OPT_TO_LkToBr, option)) {

            double lastMinTime = tLastByLk - Configor::Prior::TimeOffsetPadding;
            double lastMaxTime = tLastByLk + Configor::Prior::TimeOffsetPadding;
            // invalid time stamp
            if (!splines->TimeInRangeForSo3(lastMinTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(lastMaxTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }

            double curMinTime = tCurByLk - Configor::Prior::TimeOffsetPadding;
            double curMaxTime = tCurByLk + Configor::Prior::TimeOffsetPadding;
            // invalid time stamp
            if (!splines->TimeInRangeForSo3(curMinTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(curMaxTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }

            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                            {{lastMinTime, lastMaxTime},
                                             {curMinTime,  curMaxTime}}, so3Meta);
        } else {
            double lastTime = tLastByLk + parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
            double curTime = tCurByLk + parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);

            // check point time stamp
            if (!splines->TimeInRangeForSo3(lastTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(curTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                            {{lastTime, lastTime},
                                             {curTime,  curTime}}, so3Meta);
        }

        // create a cost function
        auto costFunc = HandEyeRotationAlignFactor<Configor::Prior::SplineOrder>::Create(
                so3Meta, tLastByLk, tCurByLk, so3LastLkToM.inverse() * so3CurLkToM, weight
        );

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }

        // SO3_LkToBr
        costFunc->AddParameterBlock(4);
        // TO_LkToBr
        costFunc->AddParameterBlock(1);

        // set Residuals
        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(
                paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                !IsOptionWith(Opt::OPT_SO3_SPLINE, option)
        );

        auto SO3_LkToBr = parMagr->EXTRI.SO3_LkToBr.at(lidarTopic).data();
        paramBlockVec.push_back(SO3_LkToBr);

        auto TO_LkToBr = &parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
        paramBlockVec.push_back(TO_LkToBr);

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

        this->SetManifold(SO3_LkToBr, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_SO3_LkToBr, option)) {
            this->SetParameterBlockConstant(SO3_LkToBr);
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
     * [ SO3 | ... | SO3 | SO3_CmToBr | TO_CmToBr ]
     */
    void Estimator::AddHandEyeRotationAlignmentForCamera(const std::string &camTopic, double tLastByCm, double tCurByCm,
                                                         const Sophus::SO3d &so3LastCmToW,
                                                         const Sophus::SO3d &so3CurCmToW, Estimator::Opt option,
                                                         double weight) {
        // prepare metas for splines
        SplineMetaType so3Meta;

        // different relative control points finding [single vs. range]
        if (IsOptionWith(Opt::OPT_TO_CmToBr, option)) {

            double lastMinTime = tLastByCm - Configor::Prior::TimeOffsetPadding;
            double lastMaxTime = tLastByCm + Configor::Prior::TimeOffsetPadding;
            // invalid time stamp
            if (!splines->TimeInRangeForSo3(lastMinTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(lastMaxTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }

            double curMinTime = tCurByCm - Configor::Prior::TimeOffsetPadding;
            double curMaxTime = tCurByCm + Configor::Prior::TimeOffsetPadding;
            // invalid time stamp
            if (!splines->TimeInRangeForSo3(curMinTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(curMaxTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }

            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                            {{lastMinTime, lastMaxTime},
                                             {curMinTime,  curMaxTime}}, so3Meta);
        } else {
            double lastTime = tLastByCm + parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
            double curTime = tCurByCm + parMagr->TEMPORAL.TO_CmToBr.at(camTopic);

            // check point time stamp
            if (!splines->TimeInRangeForSo3(lastTime, Configor::Preference::SO3_SPLINE) ||
                !splines->TimeInRangeForSo3(curTime, Configor::Preference::SO3_SPLINE)) {
                return;
            }
            splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE,
                                            {{lastTime, lastTime},
                                             {curTime,  curTime}}, so3Meta);
        }

        // create a cost function
        auto costFunc = HandEyeRotationAlignFactor<Configor::Prior::SplineOrder>::Create(
                so3Meta, tLastByCm, tCurByCm, so3LastCmToW.inverse() * so3CurCmToW, weight
        );

        // so3 knots param block [each has four sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }

        // SO3_CmToBr
        costFunc->AddParameterBlock(4);
        // TO_CmToBr
        costFunc->AddParameterBlock(1);

        // set Residuals
        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(
                paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                !IsOptionWith(Opt::OPT_SO3_SPLINE, option)
        );

        auto SO3_CmToBr = parMagr->EXTRI.SO3_CmToBr.at(camTopic).data();
        paramBlockVec.push_back(SO3_CmToBr);

        auto TO_CmToBr = &parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
        paramBlockVec.push_back(TO_CmToBr);

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

        this->SetManifold(SO3_CmToBr, QUATER_MANIFOLD.get());

        if (!IsOptionWith(Opt::OPT_SO3_CmToBr, option)) {
            this->SetParameterBlockConstant(SO3_CmToBr);
        }

        if (!IsOptionWith(Opt::OPT_TO_CmToBr, option)) {
            this->SetParameterBlockConstant(TO_CmToBr);
        } else {
            // set bound
            this->SetParameterLowerBound(TO_CmToBr, 0, -Configor::Prior::TimeOffsetPadding);
            this->SetParameterUpperBound(TO_CmToBr, 0, Configor::Prior::TimeOffsetPadding);
        }
    }

    /**
     * param blocks:
     * [ SO3 | ... | SO3 ]
     */
    void Estimator::AddSO3Constraint(double timeByBr, const Sophus::SO3d &so3, Estimator::Opt option,
                                     double weight) {
        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        // check point time stamp
        if (!so3Spline.TimeStampInRange(timeByBr) || !so3Spline.TimeStampInRange(timeByBr)) { return; }

        SplineMetaType so3Meta;
        splines->CalculateSo3SplineMeta(Configor::Preference::SO3_SPLINE, {{timeByBr, timeByBr}}, so3Meta);

        // create a cost function
        auto costFunc = SO3Factor<Configor::Prior::SplineOrder>::Create(
                so3Meta, timeByBr, so3, weight
        );

        // pos knots param block [each has three sub params]
        for (int i = 0; i < static_cast<int>(so3Meta.NumParameters()); ++i) {
            costFunc->AddParameterBlock(4);
        }

        // the Residual
        costFunc->SetNumResiduals(3);

        // organize the param block vector
        std::vector<double *> paramBlockVec;

        // so3 knots param block
        AddSo3KnotsData(
                paramBlockVec, splines->GetSo3Spline(Configor::Preference::SO3_SPLINE), so3Meta,
                !IsOptionWith(Opt::OPT_SO3_SPLINE, option)
        );

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);
    }

    /**
     * param blocks:
     * [ SO3_CurCToW | POS_CurCInW | POS_LMInW ]
     */
    void Estimator::AddVisualProjectionFactor(ns_veta::Posed *T_CurCToW, Eigen::Vector3d *POS_LMInW,
                                              const ns_veta::PinholeIntrinsic::Ptr &intri,
                                              const Eigen::Vector2d &feat, double weight) {

        // create a cost function
        auto costFunc = VisualProjFactor::Create(feat, intri, weight);

        // par blocks
        costFunc->AddParameterBlock(4);
        costFunc->AddParameterBlock(3);
        costFunc->AddParameterBlock(3);

        // the Residual
        costFunc->SetNumResiduals(2);

        // organize the param block vector
        std::vector<double *> paramBlockVec;
        paramBlockVec.push_back(T_CurCToW->Rotation().data());
        paramBlockVec.push_back(T_CurCToW->Translation().data());
        paramBlockVec.push_back(POS_LMInW->data());

        // pass to problem
        this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

        this->SetManifold(T_CurCToW->Rotation().data(), QUATER_MANIFOLD.get());
    }

    std::pair<Eigen::Vector3d, Eigen::Matrix3d>
    Estimator::InertialVelIntegration(const std::vector<IMUFrame::Ptr> &data,
                                      const std::string &imuTopic, double sTimeByBi, double eTimeByBi) {
        auto vecMatSeq = InertialIntegrationBase(data, imuTopic, sTimeByBi, eTimeByBi);
        return {TrapIntegrationOnce(vecMatSeq.first), TrapIntegrationOnce(vecMatSeq.second)};
    }

    std::pair<std::pair<Eigen::Vector3d, Eigen::Matrix3d>, std::pair<Eigen::Vector3d, Eigen::Matrix3d>>
    Estimator::InertialPosIntegration(const std::vector<IMUFrame::Ptr> &data,
                                      const std::string &imuTopic, double sTimeByBi, double eTimeByBi) {
        auto vecMatSeq = InertialIntegrationBase(data, imuTopic, sTimeByBi, eTimeByBi);
        return {{TrapIntegrationOnce(vecMatSeq.first),  TrapIntegrationOnce(vecMatSeq.second)},
                {TrapIntegrationTwice(vecMatSeq.first), TrapIntegrationTwice(vecMatSeq.second)}};
    }

    std::pair<std::vector<std::pair<double, Eigen::Vector3d>>, std::vector<std::pair<double, Eigen::Matrix3d>>>
    Estimator::InertialIntegrationBase(const std::vector<IMUFrame::Ptr> &data,
                                       const std::string &imuTopic, double sTimeByBi, double eTimeByBi) {
        // vector and matrix sequence for integration
        std::vector<std::pair<double, Eigen::Vector3d>> vecSeq;
        std::vector<std::pair<double, Eigen::Matrix3d>> matSeq;

        double timeOffset = parMagr->TEMPORAL.TO_BiToBr.at(imuTopic);
        auto SO3_BiToBr = parMagr->EXTRI.SO3_BiToBr.at(imuTopic);

        // extract data by considering the initialized time offsets
        auto [sIter, eIter] = CalibDataManager::ExtractIMUDataPiece(data, sTimeByBi, eTimeByBi);

        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);

        for (auto iter = sIter; iter != eIter; ++iter) {
            const auto &frame = *iter;
            double curTimeByBr = frame->GetTimestamp() + timeOffset;

            if (!splines->TimeInRangeForSo3(curTimeByBr, Configor::Preference::SO3_SPLINE)) { continue; }

            auto SO3_BrToBr0 = so3Spline.Evaluate(curTimeByBr);

            // angular velocity in world
            auto SO3_VEL_BrToBr0InBr0 = SO3_BrToBr0 * so3Spline.VelocityBody(curTimeByBr);
            Eigen::Matrix3d SO3_VEL_MAT = Sophus::SO3d::hat(SO3_VEL_BrToBr0InBr0);

            // angular acceleration in world
            auto SO3_ACCE_BrToBr0InBr0 = SO3_BrToBr0 * so3Spline.AccelerationBody(curTimeByBr);
            Eigen::Matrix3d SO3_ACCE_MAT = Sophus::SO3d::hat(SO3_ACCE_BrToBr0InBr0);

            // store
            vecSeq.emplace_back(curTimeByBr, SO3_BrToBr0 * SO3_BiToBr * frame->GetAcce());
            matSeq.emplace_back(curTimeByBr, (SO3_ACCE_MAT + SO3_VEL_MAT * SO3_VEL_MAT) * SO3_BrToBr0.matrix());
        }

        return {vecSeq, matSeq};
    }

    Eigen::MatrixXd Estimator::GetHessianMatrix(const std::vector<double *> &consideredParBlocks, int numThread) {
        // remove params that are not involved
        ceres::Problem::EvaluateOptions evalOpt;
        evalOpt.parameter_blocks = consideredParBlocks;
        evalOpt.num_threads = numThread;

        // evaluate
        ceres::CRSMatrix jacobianCRSMatrix;
        this->Evaluate(evalOpt, nullptr, nullptr, nullptr, &jacobianCRSMatrix);

        // obtain hessian matrix and residual vector
        Eigen::MatrixXd JMat = CRSMatrix2EigenMatrix(&jacobianCRSMatrix);

        Eigen::MatrixXd HMat = JMat.transpose() * JMat;

        return HMat;
    }

    void Estimator::SetRefIMUParamsConstant() {
        auto SO3_BiToBr = parMagr->EXTRI.SO3_BiToBr.at(Configor::DataStream::ReferIMU).data();
        if (this->HasParameterBlock(SO3_BiToBr)) {
            this->SetParameterBlockConstant(SO3_BiToBr);
        }

        auto POS_BiInBr = parMagr->EXTRI.POS_BiInBr.at(Configor::DataStream::ReferIMU).data();
        if (this->HasParameterBlock(POS_BiInBr)) {
            this->SetParameterBlockConstant(POS_BiInBr);
        }

        auto TO_BiToBr = &parMagr->TEMPORAL.TO_BiToBr.at(Configor::DataStream::ReferIMU);
        if (this->HasParameterBlock(TO_BiToBr)) {
            this->SetParameterBlockConstant(TO_BiToBr);
        }
    }

    void Estimator::FixFirSO3ControlPoint() {
        const auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        for (int i = 0; i < static_cast<int>(so3Spline.GetKnots().size()); ++i) {
            auto data = so3Spline.GetKnot(i).data();
            if (this->HasParameterBlock(data)) {
                this->SetParameterBlockConstant(data);
                break;
            }
        }
    }

    /**
     * param blocks:
     * [ example for four-order spline: VEL | VEL | VEL | VEL ]
     */
    void Estimator::AddLinScaleTailConstraint(Opt option, double weight, int count) {
        auto &velSpline = splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
        for (int j = 0; j < count - 2; ++j) {
            auto costFunc = RdLinearKnotsFactor::Create(weight);
            costFunc->AddParameterBlock(3);
            costFunc->AddParameterBlock(3);
            costFunc->AddParameterBlock(3);
            costFunc->SetNumResiduals(3);

            // organize the param block vector
            std::vector<double *> paramBlockVec(3);
            for (int i = 0; i < 3; ++i) {
                paramBlockVec.at(i) = velSpline.GetKnot(
                        j + i + static_cast<int>(velSpline.GetKnots().size()) - Configor::Prior::SplineOrder
                ).data();
            }

            this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

            if (!IsOptionWith(Opt::OPT_SCALE_SPLINE, option)) {
                for (auto &knot: paramBlockVec) { this->SetParameterBlockConstant(knot); }
            }
        }
    }

    /**
     * param blocks:
     * [ example for four-order spline: SO3 | SO3 | SO3 | SO3 ]
     */
    void Estimator::AddSO3TailConstraint(Opt option, double weight, int count) {
        auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        for (int j = 0; j < count - 2; ++j) {
            auto costFunc = So3LinearKnotsFactor::Create(weight);
            costFunc->AddParameterBlock(4);
            costFunc->AddParameterBlock(4);
            costFunc->AddParameterBlock(4);
            costFunc->SetNumResiduals(3);

            // organize the param block vector
            std::vector<double *> paramBlockVec(3);
            for (int i = 0; i < 3; ++i) {
                paramBlockVec.at(i) = so3Spline.GetKnot(
                        j + i + static_cast<int>(so3Spline.GetKnots().size()) - Configor::Prior::SplineOrder
                ).data();
            }

            this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

            for (const auto &item: paramBlockVec) { this->SetManifold(item, QUATER_MANIFOLD.get()); }

            if (!IsOptionWith(Opt::OPT_SO3_SPLINE, option)) {
                for (auto &knot: paramBlockVec) { this->SetParameterBlockConstant(knot); }
            }
        }
    }

    void Estimator::AddLinScaleHeadConstraint(Estimator::Opt option, double weight, int count) {
        auto &velSpline = splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
        for (int j = 0; j < count - 2; ++j) {
            auto costFunc = RdLinearKnotsFactor::Create(weight);
            costFunc->AddParameterBlock(3);
            costFunc->AddParameterBlock(3);
            costFunc->AddParameterBlock(3);
            costFunc->SetNumResiduals(3);

            // organize the param block vector
            std::vector<double *> paramBlockVec(3);
            for (int i = 0; i < 3; ++i) {
                paramBlockVec.at(i) = velSpline.GetKnot(j + i).data();
            }

            this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

            if (!IsOptionWith(Opt::OPT_SCALE_SPLINE, option)) {
                for (auto &knot: paramBlockVec) { this->SetParameterBlockConstant(knot); }
            }
        }
    }

    void Estimator::AddSO3HeadConstraint(Estimator::Opt option, double weight, int count) {
        auto &so3Spline = splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        for (int j = 0; j < count - 2; ++j) {
            auto costFunc = So3LinearKnotsFactor::Create(weight);
            costFunc->AddParameterBlock(4);
            costFunc->AddParameterBlock(4);
            costFunc->AddParameterBlock(4);
            costFunc->SetNumResiduals(3);

            // organize the param block vector
            std::vector<double *> paramBlockVec(3);
            for (int i = 0; i < 3; ++i) {
                paramBlockVec.at(i) = so3Spline.GetKnot(j + i).data();
            }

            this->AddResidualBlock(costFunc, nullptr, paramBlockVec);

            for (const auto &item: paramBlockVec) { this->SetManifold(item, QUATER_MANIFOLD.get()); }

            if (!IsOptionWith(Opt::OPT_SO3_SPLINE, option)) {
                for (auto &knot: paramBlockVec) { this->SetParameterBlockConstant(knot); }
            }
        }
    }
}