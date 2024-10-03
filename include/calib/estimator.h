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
#include "calib/time_deriv.hpp"
#include "ceres/ceres.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
using namespace magic_enum::bitwise_operators;
struct PointToSurfelCorr;
using PointToSurfelCorrPtr = std::shared_ptr<PointToSurfelCorr>;
struct VisualReProjCorr;
using VisualReProjCorrPtr = std::shared_ptr<VisualReProjCorr>;
struct OpticalFlowCorr;
using OpticalFlowCorrPtr = std::shared_ptr<OpticalFlowCorr>;

// myenumGenor Option OPT_SO3_SPLINE OPT_SCALE_SPLINE OPT_SO3_BiToBr OPT_POS_BiInBr
// OPT_SO3_RjToBr OPT_POS_RjInBr OPT_SO3_LkToBr OPT_POS_LkInBr OPT_SO3_CmToBr OPT_POS_CmInBr
// OPT_SO3_DnToBr OPT_POS_DnInBr OPT_TO_BiToBr OPT_TO_RjToBr OPT_TO_LkToBr OPT_TO_CmToBr
// OPT_TO_DnToBr OPT_GYRO_BIAS OPT_GYRO_MAP_COEFF OPT_ACCE_BIAS OPT_ACCE_MAP_COEFF OPT_SO3_AtoG
// OPT_GRAVITY OPT_VISUAL_GLOBAL_SCALE OPT_VISUAL_DEPTH OPT_RGBD_ALPHA
// OPT_RGBD_BETA OPT_CAM_FOCAL_LEN OPT_CAM_PRINCIPAL_POINT OPT_RS_CAM_READOUT_TIME
enum class OptOption : std::uint32_t {
    /**
     * @brief options
     */
    NONE = std::uint32_t(1) << 0,
    OPT_SO3_SPLINE = std::uint32_t(1) << 1,
    OPT_SCALE_SPLINE = std::uint32_t(1) << 2,

    OPT_SO3_BiToBr = std::uint32_t(1) << 3,
    OPT_POS_BiInBr = std::uint32_t(1) << 4,

    OPT_SO3_RjToBr = std::uint32_t(1) << 5,
    OPT_POS_RjInBr = std::uint32_t(1) << 6,

    OPT_SO3_LkToBr = std::uint32_t(1) << 7,
    OPT_POS_LkInBr = std::uint32_t(1) << 8,

    OPT_SO3_CmToBr = std::uint32_t(1) << 9,
    OPT_POS_CmInBr = std::uint32_t(1) << 10,

    OPT_SO3_DnToBr = std::uint32_t(1) << 11,
    OPT_POS_DnInBr = std::uint32_t(1) << 12,

    OPT_TO_BiToBr = std::uint32_t(1) << 13,
    OPT_TO_RjToBr = std::uint32_t(1) << 14,
    OPT_TO_LkToBr = std::uint32_t(1) << 15,
    OPT_TO_CmToBr = std::uint32_t(1) << 16,
    OPT_TO_DnToBr = std::uint32_t(1) << 17,

    OPT_GYRO_BIAS = std::uint32_t(1) << 18,
    OPT_GYRO_MAP_COEFF = std::uint32_t(1) << 19,
    OPT_ACCE_BIAS = std::uint32_t(1) << 20,
    OPT_ACCE_MAP_COEFF = std::uint32_t(1) << 21,
    OPT_SO3_AtoG = std::uint32_t(1) << 22,

    OPT_GRAVITY = std::uint32_t(1) << 23,

    OPT_VISUAL_GLOBAL_SCALE = std::uint32_t(1) << 24,
    OPT_VISUAL_DEPTH = std::uint32_t(1) << 25,

    OPT_RGBD_ALPHA = std::uint32_t(1) << 26,
    OPT_RGBD_BETA = std::uint32_t(1) << 27,

    OPT_CAM_FOCAL_LEN = std::uint32_t(1) << 28,
    OPT_CAM_PRINCIPAL_POINT = std::uint32_t(1) << 29,

    OPT_RS_CAM_READOUT_TIME = std::uint32_t(1) << 30,

    ALL = OPT_SO3_SPLINE | OPT_SCALE_SPLINE | OPT_SO3_BiToBr | OPT_POS_BiInBr | OPT_SO3_RjToBr |
          OPT_POS_RjInBr | OPT_SO3_LkToBr | OPT_POS_LkInBr | OPT_SO3_CmToBr | OPT_POS_CmInBr |
          OPT_SO3_DnToBr | OPT_POS_DnInBr | OPT_TO_BiToBr | OPT_TO_RjToBr | OPT_TO_LkToBr |
          OPT_TO_CmToBr | OPT_TO_DnToBr | OPT_GYRO_BIAS | OPT_GYRO_MAP_COEFF | OPT_ACCE_BIAS |
          OPT_ACCE_MAP_COEFF | OPT_SO3_AtoG | OPT_GRAVITY | OPT_VISUAL_GLOBAL_SCALE |
          OPT_VISUAL_DEPTH | OPT_RGBD_ALPHA | OPT_RGBD_BETA | OPT_CAM_FOCAL_LEN |
          OPT_CAM_PRINCIPAL_POINT | OPT_RS_CAM_READOUT_TIME
};

struct SpatialTemporalPriori;
using SpatialTemporalPrioriPtr = std::shared_ptr<SpatialTemporalPriori>;

class Estimator : public ceres::Problem {
public:
    using Ptr = std::shared_ptr<Estimator>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;
    using SplineMetaType = ns_ctraj::SplineMeta<Configor::Prior::SplineOrder>;
    using Opt = OptOption;

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
                               double acceWeight);

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

    void AddVelVisualInertialAlignment(const std::vector<IMUFrame::Ptr> &data,
                                       const std::string &imuTopic,
                                       const std::string &topic,
                                       const std::pair<CameraFrame::Ptr, Eigen::Vector3d> &sVelAry,
                                       double *sVelScale,
                                       const std::pair<CameraFrame::Ptr, Eigen::Vector3d> &eVelAry,
                                       double *eVelScale,
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
                             double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_RjToBr | POS_RjInBr | TO_RjToBr ]
     */
    template <int TimeDeriv>
    void AddLinearScaleConstraint(double timeByBr,
                                  const Eigen::Vector3d &linScaleOfDeriv,
                                  Estimator::Opt option,
                                  double weight);

    void AddSO3Constraint(double timeByBr,
                          const Sophus::SO3d &so3,
                          Estimator::Opt option,
                          double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_LkToBr | POS_LkInBr | TO_LkToBr ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddLiDARPointTiSurfelConstraint(const PointToSurfelCorrPtr &ptsCorr,
                                         const std::string &topic,
                                         Opt option,
                                         double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_DnToBr | POS_DnInBr | TO_DnToBr ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddRGBDPointTiSurfelConstraint(const PointToSurfelCorrPtr &ptsCorr,
                                        const std::string &topic,
                                        Opt option,
                                        double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
     * READOUT_TIME | FX | FY | CX | CY | GLOBAL_SCALE | INV_DEPTH ]
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddVisualReprojection(const VisualReProjCorrPtr &visualCorr,
                               const std::string &topic,
                               double *globalScale,
                               double *invDepth,
                               Opt option,
                               double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_DnToBr | POS_DnInBr | TO_DnToBr |
     *   READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    void AddRGBDOpticalFlowConstraint(const OpticalFlowCorrPtr &ofCorr,
                                      const std::string &topic,
                                      Opt option,
                                      double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
     *   READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    void AddVisualOpticalFlowConstraint(const OpticalFlowCorrPtr &ofCorr,
                                        const std::string &topic,
                                        Opt option,
                                        double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_CmToBr | POS_CmInBr | TO_CmToBr |
     * READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    void AddVisualOpticalFlowReprojConstraint(const OpticalFlowCorrPtr &velCorr,
                                              const std::string &topic,
                                              Opt option,
                                              double weight);

    /**
     * param blocks:
     * [ SO3 | ... | SO3 | LIN_SCALE | ... | LIN_SCALE | SO3_DnToBr | POS_DnInBr | TO_DnToBr |
     *   READOUT_TIME | FX | FY | CX | CY | DEPTH_INFO ]
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    void AddRGBDOpticalFlowReprojConstraint(const OpticalFlowCorrPtr &velCorr,
                                            const std::string &topic,
                                            Opt option,
                                            double weight);

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

    void PrintUninvolvedKnots() const;

    void AddVisualVelocityDepthFactor(Eigen::Vector3d *LIN_VEL_CmToWInCm,
                                      const OpticalFlowCorrPtr &corr,
                                      double TO_CamToBr,
                                      double readout,
                                      const Sophus::SO3d &SO3_CamToBr,
                                      const ns_veta::PinholeIntrinsic::Ptr &intri,
                                      double weight,
                                      bool estDepth,
                                      bool estVelDirOnly);

    void AddVisualVelocityDepthFactorForRGBD(Eigen::Vector3d *LIN_VEL_CmToWInCm,
                                             const OpticalFlowCorrPtr &corr,
                                             const std::string &rgbdTopic,
                                             double weight,
                                             bool estDepth,
                                             bool estVelDirOnly);

    void AddVisualVelocityDepthFactorForVelCam(Eigen::Vector3d *LIN_VEL_CmToWInCm,
                                               const OpticalFlowCorrPtr &corr,
                                               const std::string &topic,
                                               double weight,
                                               bool estDepth,
                                               bool estVelDirOnly);

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

    std::optional<std::pair<Eigen::Vector3d, Eigen::Matrix3d>> InertialVelIntegration(
        const std::vector<IMUFrame::Ptr> &data,
        const std::string &imuTopic,
        double sTimeByBi,
        double eTimeByBi);

    std::optional<std::pair<std::pair<Eigen::Vector3d, Eigen::Matrix3d>,
                            std::pair<Eigen::Vector3d, Eigen::Matrix3d>>>
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

    /**
     * compute the time range of knots to be considered in optimization based on given information
     * @param timeByCam the time stamped by the camera
     * @param RS_READOUT the readout time of rs camera
     * @param RT_PADDING the time padding of the readout time
     * @param rdFactor the readout factor
     * @param optReadout whether optimize the readout time
     * @param TO_CmToBr the time offset of the camera with respect to the reference IMU
     * @param TO_PADDING the time padding og the time offset
     * @param optTimeOffset whether optimize the time offset
     * @return the time range [min time, max time]
     */
    static std::pair<double, double> ConsideredTimeRangeForCameraStamp(double timeByCam,
                                                                       double RS_READOUT,
                                                                       double RT_PADDING,
                                                                       double rdFactor,
                                                                       bool optReadout,
                                                                       double TO_CmToBr,
                                                                       double TO_PADDING,
                                                                       bool optTimeOffset);

    /**
     * check whether time range (mainly from 'ConsideredTimeRangeForCameraStamp') is valid foe splines
     * @param timePair the time stamp pair
     * @return true: valid, false invalid
     */
    [[nodiscard]] bool TimeInRangeForSplines(const std::pair<double, double> &timePair) const;
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_ESTIMATOR_H
