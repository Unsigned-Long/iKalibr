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

#ifndef IKALIBR_CALIB_SOLVER_H
#define IKALIBR_CALIB_SOLVER_H

#include "calib/time_deriv.hpp"
#include "ceres/solver.h"
#include "config/configor.h"
#include "core/rot_only_vo.h"
#include "ctraj/core/pose.hpp"
#include "ctraj/core/spline_bundle.h"
#include "optional"
#include "pcl/point_cloud.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_veta {
class Veta;
using VetaPtr = std::shared_ptr<Veta>;
}  // namespace ns_veta

namespace pcl {
struct PointXYZRGBA;
using ColorPointCloudPtr = boost::shared_ptr<pcl::PointCloud<PointXYZRGBA>>;
}  // namespace pcl

struct PointXYZT;
using IKalibrPointCloudPtr = boost::shared_ptr<pcl::PointCloud<PointXYZT>>;

namespace ns_ikalibr {
struct CalibParamManager;
using CalibParamManagerPtr = std::shared_ptr<CalibParamManager>;
class CalibDataManager;
using CalibDataManagerPtr = std::shared_ptr<CalibDataManager>;
struct OpticalFlowTripleTrace;
using OpticalFlowTripleTracePtr = std::shared_ptr<OpticalFlowTripleTrace>;
struct SpatialTemporalPriori;
using SpatialTemporalPrioriPtr = std::shared_ptr<SpatialTemporalPriori>;
struct LiDAROdometer;
using LiDAROdometerPtr = std::shared_ptr<LiDAROdometer>;
struct VisualReProjCorrSeq;
using VisualReProjCorrSeqPtr = std::shared_ptr<VisualReProjCorrSeq>;
struct OpticalFlowCorr;
using OpticalFlowCorrPtr = std::shared_ptr<OpticalFlowCorr>;
struct PointToSurfelCorr;
using PointToSurfelCorrPtr = std::shared_ptr<PointToSurfelCorr>;
struct LiDARFrame;
using LiDARFramePtr = std::shared_ptr<LiDARFrame>;
class CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;
class Viewer;
using ViewerPtr = std::shared_ptr<Viewer>;
class Estimator;
using EstimatorPtr = std::shared_ptr<Estimator>;
enum class OptOption : std::uint32_t;

struct ImagesInfo {
public:
    std::string topic;
    std::string root_path;
    // image list [id, filename]
    std::map<ns_veta::IndexT, std::string> images;

    ImagesInfo(std::string topic,
               std::string rootPath,
               const std::map<ns_veta::IndexT, std::string> &images)
        : topic(std::move(topic)),
          root_path(std::move(rootPath)),
          images(images) {}

    [[nodiscard]] std::optional<std::string> GetImagePath(ns_veta::IndexT id) const;

    [[nodiscard]] std::optional<std::string> GetImageFilename(ns_veta::IndexT id) const;

    [[nodiscard]] std::map<ns_veta::IndexT, std::string> GetImagesIdxToName() const;

    [[nodiscard]] std::map<std::string, ns_veta::IndexT> GetImagesNameToIdx() const;

public:
    template <class Archive>
    void serialize(Archive &ar) {
        ar(CEREAL_NVP(topic), CEREAL_NVP(root_path), CEREAL_NVP(images));
    }
};

class CalibSolver {
public:
    using Ptr = std::shared_ptr<CalibSolver>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;

    friend class CalibSolverIO;

    struct BackUp {
    public:
        using Ptr = std::shared_ptr<BackUp>;

    public:
        // estimator
        EstimatorPtr estimator;
        // visual global scale
        std::shared_ptr<double> visualGlobalScale;
        // visual reprojection correspondences contains inverse depth parameters
        std::map<std::string, std::vector<VisualReProjCorrSeqPtr>> visualCorrs;
        // lidar global map
        IKalibrPointCloudPtr lidarMap;
        // lidar point-to-surfel correspondences
        std::map<std::string, std::vector<PointToSurfelCorrPtr>> lidarCorrs;
        // radar global map
        IKalibrPointCloudPtr radarMap;
        // visual optical flow correspondences, orienting to RGBDs and VelCameras
        std::map<std::string, std::vector<OpticalFlowCorrPtr>> ofCorrs;
    };

    struct InitAsset {
    public:
        using Ptr = std::shared_ptr<InitAsset>;

    public:
        // lidar-only odometers for each lidar
        std::map<std::string, LiDAROdometerPtr> lidarOdometers;
        // undistorted scans for each lidar
        std::map<std::string, std::vector<LiDARFramePtr>> undistFramesInScan;
        // rgbd-derived rgbd-frame velocities for each rgbd camera
        std::map<std::string, std::vector<std::pair<CameraFramePtr, Eigen::Vector3d>>>
            rgbdBodyFrameVels;
        // camera-derived rgbd-frame velocities for each vel-based camera
        std::map<std::string, std::vector<std::pair<CameraFramePtr, Eigen::Vector3d>>>
            velCamBodyFrameVelDirs;
        // camera-derived rgbd-frame velocities for each vel-based camera
        std::map<std::string, std::vector<std::pair<CameraFramePtr, Eigen::Vector3d>>>
            velCamBodyFrameVels;
        // SfM pose sequence for each camera
        std::map<std::string, std::vector<ns_ctraj::Posed>> sfmPoseSeq;
        // the global lidar map expressed in world frame
        IKalibrPointCloudPtr globalMap;
        // undistorted scans for each lidar expressed in the world frame
        std::map<std::string, std::vector<LiDARFramePtr>> undistFramesInMap;
    };

private:
    // the data manager for calibration
    CalibDataManagerPtr _dataMagr;
    // the parameter manager for calibration
    CalibParamManagerPtr _parMagr;
    // prior knowledge about spatiotemporal parameters used in optimization (if provided)
    SpatialTemporalPrioriPtr _priori;
    // rotation and linear scale splines (acceleration, velocity, or translation)
    SplineBundleType::Ptr _splines;
    // options used for ceres-related optimization
    ceres::Solver::Options _ceresOption;
    // viewer used to visualize entities in calibration
    ViewerPtr _viewer;
    // storge results from optimization for by-products-related output
    BackUp::Ptr _backup;
    // storge temporal results from initialization, which would be destroyed after initialization
    InitAsset::Ptr _initAsset;
    // indicates whether the solving is finished
    bool _solveFinished;

public:
    /**
     * create a solver for spatiotemporal calibration
     * @param calibDataManager the data manager
     * @param calibParamManager the parameter manager
     */
    explicit CalibSolver(CalibDataManagerPtr calibDataManager,
                         CalibParamManagerPtr calibParamManager);

    /**
     * create a solver shared pointer for spatiotemporal calibration
     * @param calibDataManager the data manager
     * @param calibParamManager the parameter manager
     * @return the shared pointer of this solver
     */
    static Ptr Create(const CalibDataManagerPtr &calibDataManager,
                      const CalibParamManagerPtr &calibParamManager);

    /**
     * perform the spatiotemporal calibration
     */
    void Process();

    /**
     * de-constructor
     */
    virtual ~CalibSolver();

protected:
    /**
     * transform an input veta using given transformation information, if scale is provide,
     * this veta would also ve scaled
     * @param veta the visual meta data to be transformed
     * @param curToNew the transformation information from current coordinate frame to new one
     * @param scale the scale factor, if no scaling is needed, set it to 1.0
     */
    static void PerformTransformForVeta(const ns_veta::VetaPtr &veta,
                                        const ns_veta::Posed &curToNew,
                                        double scale);

    /**
     * align vectors to a new coordinate frame where gravity pointing to negative z-axis.
     * the splines (both rotation and translation splines), as well as the gravity vector would be
     * transformed
     */
    void AlignStatesToGravity() const;

    /**
     * initialize (recover) the rotation spline using raw angular velocity measurements from
     * the gyroscope. If multiple gyroscopes (IMUs) are involved, the extrinsic rotations and
     * time offsets would be also recovered
     */
    void InitSO3Spline() const;

    /**
     * perform sensor-inertial alignment to recover the gravity vector and extrinsic translations
     */
    void InitSensorInertialAlign() const;

    /**
     * detailed sensor-inertial alignment for camera and IMU, this is the preparation for final
     * one-shot sensor-inertial alignment
     */
    void InitPrepPosCameraInertialAlign() const;

    /**
     * detailed sensor-inertial alignment for camera and IMU, this is the preparation for final
     * one-shot sensor-inertial alignment
     */
    void InitPrepVelCameraInertialAlign() const;

    /**
     * detailed sensor-inertial alignment for LiDAR and IMU, this is the preparation for final
     * one-shot sensor-inertial alignment
     */
    void InitPrepLiDARInertialAlign() const;

    /**
     * detailed sensor-inertial alignment for RGBD camera and IMU, this is the preparation for final
     * one-shot sensor-inertial alignment
     */
    void InitPrepRGBDInertialAlign() const;

    /**
     * detailed sensor-inertial alignment for radars and IMU, this is the preparation for final
     * one-shot sensor-inertial alignment
     */
    void InitPrepRadarInertialAlign();

    /**
     * detailed sensor-inertial alignment for IMUs, i.e., inertial-inertial alignment,
     * this is the preparation for final one-shot sensor-inertial alignment
     */
    void InitPrepInertialInertialAlign();

    /**
     * initialize the linear scale spline using by-products from sensor-inertial alignment
     * one of three kinds of linear scale spline, i.e., linear acceleration, linear velocity, and
     * translation splines, would be recovered
     */
    void InitScaleSpline() const;

    /**
     * preparation for the final batch optimization
     */
    void InitPrepBatchOpt() const;

    /**
     * build the global map for LiDARs
     * @return the final map and frames used for map construction.
     * these frames are also expressed in the map (global) frame, rather than local frames
     */
    std::tuple<IKalibrPointCloudPtr, std::map<std::string, std::vector<LiDARFramePtr>>>
    BuildGlobalMapOfLiDAR() const;

    /**
     * build the global map for radars, this is only for visualization if translation spline
     * is employed
     * @return the global map for radars expressed in the global coordinate frame
     */
    IKalibrPointCloudPtr BuildGlobalMapOfRadar() const;

    /**
     * build the global colorized map for RGBD cameras
     * @param topic ros topic for the RGBD
     * @return the global colorized map
     */
    pcl::ColorPointCloudPtr BuildGlobalColorMapOfRGBD(const std::string &topic) const;

    /**
     * build the global map for RGBDs
     * @return the global map, scans in global frame, and scans in local frame
     */
    std::tuple<IKalibrPointCloudPtr,
               // scans in global frame
               std::map<std::string, std::vector<IKalibrPointCloudPtr>>,
               // scans in local frame
               std::map<std::string, std::vector<IKalibrPointCloudPtr>>>
    BuildGlobalMapOfRGBD() const;

    /**
     * @param topic the ros topic of this vision sensor
     * @param ofVec the optical flow correspondence
     * @param intri the visual intrinsics
     * @param SE3_CurSenToW the function to obtain pose of the vision sensor given a timestamp and
     * the ros topic
     * @return the vision meta data structure
     */
    ns_veta::VetaPtr CreateVetaFromOpticalFlow(
        const std::string &topic,
        const std::vector<OpticalFlowCorrPtr> &ofVec,
        const ns_veta::PinholeIntrinsic::Ptr &intri,
        const std::function<std::optional<Sophus::SE3d>(
            const CalibSolver *, double, const std::string &)> &SE3_CurSenToW) const;

    /**
     * perform data association for LiDARs
     * @param map the global point cloud map
     * @param undistFrames the undistorted scans expressed in the global coordinate frame
     * @param ptsCountInEachScan construct how many correspondences in each scan
     * @return the point-to-surfel correspondences for each LiDAR
     */
    std::map<std::string, std::vector<PointToSurfelCorrPtr>> DataAssociationForLiDARs(
        const IKalibrPointCloudPtr &map,
        const std::map<std::string, std::vector<LiDARFramePtr>> &undistFrames,
        int ptsCountInEachScan) const;

    /**
     * perform data association for pos-derived cameras
     * @return the visual reprojection correspondences for each optical camera
     */
    std::map<std::string, std::vector<VisualReProjCorrSeqPtr>> DataAssociationForPosCameras() const;

    /**
     * perform data association for vel-derived cameras
     * @return the optical flow correspondences for each optical camera
     */
    std::map<std::string, std::vector<OpticalFlowCorrPtr>> DataAssociationForVelCameras() const;

    /**
     * perform data association for RGBD cameras
     * @param estDepth whether estimate depth of the constructed correspondences
     * @return the RGBD optical flow correspondence
     */
    std::map<std::string, std::vector<OpticalFlowCorrPtr>> DataAssociationForRGBDs(bool estDepth);

    /**
     * perform data associate for RGBD cameras, I mean the point-to-surfel data association.
     * If the depth camera is well-matched with the RGB camera, then use these correspondence
     * would help spatiotemporal calibration. Otherwise, this would introduce additional errors.
     * @param map the global point cloud map
     * @param scanInGFrame the scans in the global coordinate frame
     * @param scanInLFrame the scans in the local coordinate frame
     * @param ptsCountInEachScan construct how many correspondences in each scan
     * @return the point-to-surfel correspondences for each RGBD camera
     */
    std::map<std::string, std::vector<PointToSurfelCorrPtr>> DataAssociationForRGBDs(
        const IKalibrPointCloudPtr &map,
        const std::map<std::string, std::vector<IKalibrPointCloudPtr>> &scanInGFrame,
        const std::map<std::string, std::vector<IKalibrPointCloudPtr>> &scanInLFrame,
        int ptsCountInEachScan) const;

    /**
     * the final continuous-time-based batch optimization
     * @param optOption the option for optimization, deciding which variable (state) would be
     * estimated
     * @param lidarPtsCorrs the point-to-surfel correspondences for LiDARs
     * @param visualReprojCorrs the visual reprojection correspondence for cameras
     * @param rgbdCorrs the optical flow correspondence for RGBDs
     * @param visualVelCorrs the optical flow correspondence for vel-powered cameras
     * @param rgbdPtsCorrs the point-to-surfel correspondences for RGBDs, its optional
     * @return the backup data from batch optimization
     */
    BackUp::Ptr BatchOptimization(
        OptOption optOption,
        const std::map<std::string, std::vector<PointToSurfelCorrPtr>> &lidarPtsCorrs,
        const std::map<std::string, std::vector<VisualReProjCorrSeqPtr>> &visualReprojCorrs,
        const std::map<std::string, std::vector<OpticalFlowCorrPtr>> &rgbdCorrs,
        const std::map<std::string, std::vector<OpticalFlowCorrPtr>> &visualVelCorrs,
        const std::optional<std::map<std::string, std::vector<PointToSurfelCorrPtr>>>
            &rgbdPtsCorrs = std::nullopt) const;

    /**
     * compute the pose of IMU in the global (world) coordinate frame
     * @param timeByBr the time stamped by the reference IMU
     * @return the reference IMU pose, if the timestamp is out of range, return 'std::nullopt'
     */
    std::optional<Sophus::SE3d> CurBrToW(double timeByBr) const;

    /**
     * compute the pose of LiDAR in the global (world) coordinate frame
     * @param timeByLk the time stamped by the LiDAR, i.e., the raw timestamp
     * @param topic the ros topic of this LiDAR
     * @return the LiDAR pose, if the timestamp is out of range, return 'std::nullopt'
     */
    std::optional<Sophus::SE3d> CurLkToW(double timeByLk, const std::string &topic) const;

    /**
     * compute the pose of camera in the global (world) coordinate frame
     * @param timeByCm the time stamped by the camera, i.e., the raw timestamp
     * @param topic the ros topic of this camera
     * @return the camera pose, if the timestamp is out of range, return 'std::nullopt'
     */
    std::optional<Sophus::SE3d> CurCmToW(double timeByCm, const std::string &topic) const;

    /**
     * compute the pose of RGBD camera in the global (world) coordinate frame
     * @param timeByDn the time stamped by the RGBD camera, i.e., the raw timestamp
     * @param topic the ros topic of this RGBD camera
     * @return the RGBD camera pose, if the timestamp is out of range, return 'std::nullopt'
     */
    std::optional<Sophus::SE3d> CurDnToW(double timeByDn, const std::string &topic) const;

    /**
     * compute the pose of radar in the global (world) coordinate frame
     * @param timeByRj the time stamped by the radar, i.e., the raw timestamp
     * @param topic the ros topic of this radar
     * @return the radar pose, if the timestamp is out of range, return 'std::nullopt'
     */
    std::optional<Sophus::SE3d> CurRjToW(double timeByRj, const std::string &topic) const;

    /**
     * create a spline bundle, including the rotation and linear scale splines
     * @param st the starr timestamp
     * @param et the end timestamp
     * @param so3Dt the time distance between two rotation control points
     * @param scaleDt the time distance between two linear scale control points
     * @return the created spline bundle
     */
    static SplineBundleType::Ptr CreateSplineBundle(double st,
                                                    double et,
                                                    double so3Dt,
                                                    double scaleDt);

    /**
     * get the type of the linear scale spline, it can be linear acceleration, linear velocity,
     * and translation spline, decided by the sensor suite to be calibrated
     * @return the linear scale spline type
     */
    static TimeDeriv::ScaleSplineType GetScaleType();

    /**
     * add radar Dopple velocity factors for the radar to the estimator
     * @tparam type the linear scale spline type
     * @param estimator the estimator
     * @param radarTopic the ros topic of this radar
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddRadarFactor(EstimatorPtr &estimator,
                        const std::string &radarTopic,
                        OptOption option) const;

    /**
     * add accelerometer factors for the IMU to the estimator
     * @tparam type the linear scale spline type
     * @param estimator the estimator
     * @param imuTopic the ros topic of this IMU
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type>
    void AddAcceFactor(EstimatorPtr &estimator,
                       const std::string &imuTopic,
                       OptOption option) const;

    /**
     * add gyroscope factors for the IMU to the estimator
     * @param estimator the estimator
     * @param imuTopic the ros topic of this IMU
     * @param option the option for the optimization
     */
    void AddGyroFactor(EstimatorPtr &estimator,
                       const std::string &imuTopic,
                       OptOption option) const;

    /**
     * add point-to-surfel factors for the LiDAR to the estimator
     * @tparam type the linear scale spline type
     * @param estimator the estimator
     * @param lidarTopic the ros topic of this LiDAR
     * @param corrs the point-to-surfel correspondences
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type>
    static void AddLiDARPointToSurfelFactor(EstimatorPtr &estimator,
                                            const std::string &lidarTopic,
                                            const std::vector<PointToSurfelCorrPtr> &corrs,
                                            OptOption option);

    /**
     * add point-to-surfel factors for the RGBD to the estimator
     * @tparam type the linear scale spline type
     * @param estimator the estimator
     * @param rgbdTopic the ros topic of this RGBD
     * @param corrs the point-to-surfel correspondences
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type>
    static void AddRGBDPointToSurfelFactor(EstimatorPtr &estimator,
                                           const std::string &rgbdTopic,
                                           const std::vector<PointToSurfelCorrPtr> &corrs,
                                           OptOption option);

    /**
     * add visual reprojection factors for the camera to the estimator
     * @tparam type the linear scale spline type
     * @param estimator the estimator
     * @param camTopic the ros topic of this camera
     * @param corrs the visual reprojection correspondences
     * @param globalScale the global scale factor  of the visual structure of this camera
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type>
    static void AddVisualReprojectionFactor(EstimatorPtr &estimator,
                                            const std::string &camTopic,
                                            const std::vector<VisualReProjCorrSeqPtr> &corrs,
                                            double *globalScale,
                                            OptOption option);

    /**
     * add optical flow factors for the RGBD camera to the estimator
     * @tparam type the linear scale spline type
     * @tparam IsInvDepth estimate the depth or the inverse depth
     * @param estimator the estimator
     * @param rgbdTopic the ros topic of this RGBD camera
     * @param corrs the optical flow correspondences
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    static void AddRGBDOpticalFlowFactor(EstimatorPtr &estimator,
                                         const std::string &rgbdTopic,
                                         const std::vector<OpticalFlowCorrPtr> &corrs,
                                         OptOption option);

    /**
     * add optical flow factors for the optical camera to the estimator
     * @tparam type the linear scale spline type
     * @tparam IsInvDepth estimate the depth or the inverse depth
     * @param estimator the estimator
     * @param camTopic the ros topic of this camera
     * @param corrs the optical flow correspondences
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    static void AddVisualOpticalFlowFactor(EstimatorPtr &estimator,
                                           const std::string &camTopic,
                                           const std::vector<OpticalFlowCorrPtr> &corrs,
                                           OptOption option);

    /**
     * add optical flow factors for the optical camera to the estimator
     * @tparam type the linear scale spline type
     * @tparam IsInvDepth estimate the depth or the inverse depth
     * @param estimator the estimator
     * @param camTopic the ros topic of this camera
     * @param corrs the optical flow correspondences
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    static void AddVisualOpticalFlowReprojFactor(EstimatorPtr &estimator,
                                                 const std::string &camTopic,
                                                 const std::vector<OpticalFlowCorrPtr> &corrs,
                                                 OptOption option);

    /**
     * add optical flow factors for the optical camera to the estimator
     * @tparam type the linear scale spline type
     * @tparam IsInvDepth estimate the depth or the inverse depth
     * @param estimator the estimator
     * @param rgbdTopic the ros topic of this camera
     * @param corrs the optical flow correspondences
     * @param option the option for the optimization
     */
    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    static void AddRGBDOpticalFlowReprojFactor(EstimatorPtr &estimator,
                                               const std::string &rgbdTopic,
                                               const std::vector<OpticalFlowCorrPtr> &corrs,
                                               OptOption option);

    /**
     * store images to the disk for structure from motion (SfM)
     * @param camTopic the ros topic of this camera
     * @param matchRes the match results of images of this camera
     */
    void StoreImagesForSfM(
        const std::string &camTopic,
        const std::set<std::pair<ns_veta::IndexT, ns_veta::IndexT>> &matchRes) const;

    /**
     * once the SfM is performed, we load the results from the disk to the program
     * @param camTopic the ros topic of this camera
     * @param errorThd the reprojection error threshold
     * @param trackLenThd the track length threshold
     * @return the veta of this camera, storaging the visual structure
     */
    ns_veta::VetaPtr TryLoadSfMData(const std::string &camTopic,
                                    double errorThd,
                                    std::size_t trackLenThd) const;

    /**
     * downsample the landmarks for a veta
     * @param veta the visual meta data
     * @param lmNumThd the landmark number threshold
     * @param obvNumThd the observation number threshold
     */
    static void DownsampleVeta(const ns_veta::VetaPtr &veta,
                               std::size_t lmNumThd,
                               std::size_t obvNumThd);

    /**
     * is a camera a RS camera
     * @param camTopic the ros topic of this camera
     * @return if this camera is a RS camera
     */
    static bool IsRSCamera(const std::string &camTopic);

    /**
     * save spatiotemporal calibration results to disk
     * @param par the spatiotemporal parameter manager
     * @param desc the description of this spatiotemporal parameter
     */
    static void SaveStageCalibParam(const CalibParamManagerPtr &par, const std::string &desc);

    /**
     * create optical flow triple trace for RGBD cameras
     * @param trackInfoList the optical tracking information
     * @param trackThd the threshold of track length
     * @return the optical flow trace of this RGBD camera
     */
    static std::vector<OpticalFlowTripleTracePtr> CreateOpticalFlowTrace(
        const std::list<RotOnlyVisualOdometer::FeatTrackingInfo> &trackInfoList, int trackThd);
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_CALIB_SOLVER_H
