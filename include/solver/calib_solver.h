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

#include "calib/calib_data_manager.h"
#include "calib/calib_param_manager.h"
#include "calib_solver_io.h"
#include "calib/estimator.h"
#include "core/rot_only_vo.h"
#include "viewer/viewer.h"
#include "optional"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct OpticalFlowTripleTrace;
using OpticalFlowTripleTracePtr = std::shared_ptr<OpticalFlowTripleTrace>;
struct SpatialTemporalPriori;
using SpatialTemporalPrioriPtr = std::shared_ptr<SpatialTemporalPriori>;
struct LiDAROdometer;
using LiDAROdometerPtr = std::shared_ptr<LiDAROdometer>;
struct VisualReProjCorrSeq;
using VisualReProjCorrSeqPtr = std::shared_ptr<VisualReProjCorrSeq>;

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
        Estimator::Ptr estimator;
        // visual global scale
        std::shared_ptr<double> visualGlobalScale;
        // visual reprojection correspondences contains inverse depth parameters
        std::map<std::string, std::vector<VisualReProjCorrSeqPtr>> visualCorrs;
        // lidar global map
        IKalibrPointCloud::Ptr lidarMap;
        // lidar point-to-surfel correspondences
        std::map<std::string, std::vector<PointToSurfelCorrPtr>> lidarCorrs;
        // radar global map
        IKalibrPointCloud::Ptr radarMap;
        // rgbd velocity correspondences
        std::map<std::string, std::vector<RGBDVelocityCorrPtr>> rgbdCorrs;
    };

    struct InitAsset {
    public:
        using Ptr = std::shared_ptr<InitAsset>;

    public:
        // lidar-only odometers for each lidar
        std::map<std::string, LiDAROdometerPtr> lidarOdometers;
        // undistorted scans for each lidar
        std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFramesInScan;
        // rgbd-derived rgbd-frame velocities for each rgbd camera
        std::map<std::string, std::vector<std::pair<CameraFrame::Ptr, Eigen::Vector3d>>>
            rgbdBodyFrameVels;
        // SfM pose sequence for each camera
        std::map<std::string, std::vector<ns_ctraj::Posed>> sfmPoseSeq;
        // the global lidar map expressed in world frame
        IKalibrPointCloud::Ptr globalMap;
        // undistorted scans for each lidar expressed in the world frame
        std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFramesInMap;
    };

private:
    CalibDataManager::Ptr _dataMagr;
    CalibParamManager::Ptr _parMagr;

    SpatialTemporalPrioriPtr _priori;

    SplineBundleType::Ptr _splines;
    ceres::Solver::Options _ceresOption;

    Viewer::Ptr _viewer;
    BackUp::Ptr _backup;
    InitAsset::Ptr _initAsset;

    bool _solveFinished;

public:
    explicit CalibSolver(CalibDataManager::Ptr calibDataManager,
                         CalibParamManager::Ptr calibParamManager);

    static CalibSolver::Ptr Create(const CalibDataManager::Ptr &calibDataManager,
                                   const CalibParamManager::Ptr &calibParamManager);

    void Process();

    virtual ~CalibSolver();

protected:
    static void PerformTransformForVeta(const ns_veta::Veta::Ptr &veta,
                                        const ns_veta::Posed &curToNew,
                                        double scale);

    void AlignStatesToGravity();

    void InitSO3Spline();

    void InitSensorInertialAlign();

    void InitPrepCameraInertialAlign();

    void InitPrepLiDARInertialAlign();

    void InitPrepRGBDInertialAlign();

    void InitPrepRadarInertialAlign();

    void InitScaleSpline();

    void InitPrepBatchOpt();

    std::tuple<IKalibrPointCloud::Ptr, std::map<std::string, std::vector<LiDARFrame::Ptr>>>
    BuildGlobalMapOfLiDAR();

    IKalibrPointCloud::Ptr BuildGlobalMapOfRadar();

    ColorPointCloud::Ptr BuildGlobalColorMapOfRGBD(const std::string &topic);

    std::tuple<IKalibrPointCloud::Ptr,
               // scans in global frame
               std::map<std::string, std::vector<IKalibrPointCloud::Ptr>>,
               // scans in local frame
               std::map<std::string, std::vector<IKalibrPointCloud::Ptr>>>
    BuildGlobalMapOfRGBD();

    ns_veta::Veta::Ptr CreateVetaFromRGBD(const std::string &topic);

    std::map<std::string, std::vector<PointToSurfelCorrPtr>> DataAssociationForLiDARs(
        const IKalibrPointCloud::Ptr &map,
        const std::map<std::string, std::vector<LiDARFrame::Ptr>> &undistFrames,
        int ptsCountInEachScan);

    std::map<std::string, std::vector<VisualReProjCorrSeqPtr>> DataAssociationForCameras();

    std::map<std::string, std::vector<RGBDVelocityCorrPtr>> DataAssociationForRGBDs(bool estDepth);

    std::map<std::string, std::vector<PointToSurfelCorrPtr>> DataAssociationForRGBDs(
        const IKalibrPointCloud::Ptr &map,
        const std::map<std::string, std::vector<IKalibrPointCloud::Ptr>> &scanInGFrame,
        const std::map<std::string, std::vector<IKalibrPointCloud::Ptr>> &scanInLFrame,
        int ptsCountInEachScan);

    BackUp::Ptr BatchOptimization(
        OptOption::Option optOption,
        const std::map<std::string, std::vector<PointToSurfelCorrPtr>> &lidarPtsCorrs,
        const std::map<std::string, std::vector<VisualReProjCorrSeqPtr>> &visualCorrs,
        const std::map<std::string, std::vector<RGBDVelocityCorrPtr>> &rgbdCorrs,
        const std::optional<std::map<std::string, std::vector<PointToSurfelCorrPtr>>>
            &rgbdPtsCorrs = std::nullopt);

    std::optional<Sophus::SE3d> CurBrToW(double timeByBr);

    std::optional<Sophus::SE3d> CurLkToW(double timeByLk, const std::string &topic);

    std::optional<Sophus::SE3d> CurCmToW(double timeByCm, const std::string &topic);

    std::optional<Sophus::SE3d> CurDnToW(double timeByDn, const std::string &topic);

    std::optional<Sophus::SE3d> CurRjToW(double timeByRj, const std::string &topic);

    static SplineBundleType::Ptr CreateSplineBundle(double st,
                                                    double et,
                                                    double so3Dt,
                                                    double scaleDt);

    static TimeDeriv::ScaleSplineType GetScaleType();

    template <TimeDeriv::ScaleSplineType type>
    void AddRadarFactor(Estimator::Ptr &estimator,
                        const std::string &radarTopic,
                        Estimator::Opt option);

    template <TimeDeriv::ScaleSplineType type>
    void AddAcceFactor(Estimator::Ptr &estimator,
                       const std::string &imuTopic,
                       Estimator::Opt option);

    void AddGyroFactor(Estimator::Ptr &estimator,
                       const std::string &imuTopic,
                       Estimator::Opt option);

    template <TimeDeriv::ScaleSplineType type>
    void AddLiDARPointToSurfelFactor(Estimator::Ptr &estimator,
                                     const std::string &lidarTopic,
                                     const std::vector<PointToSurfelCorrPtr> &corrs,
                                     Estimator::Opt option);

    template <TimeDeriv::ScaleSplineType type>
    void AddRGBDPointToSurfelFactor(Estimator::Ptr &estimator,
                                    const std::string &rgbdTopic,
                                    const std::vector<PointToSurfelCorrPtr> &corrs,
                                    Estimator::Opt option);

    template <TimeDeriv::ScaleSplineType type>
    void AddVisualReprojectionFactor(Estimator::Ptr &estimator,
                                     const std::string &camTopic,
                                     const std::vector<VisualReProjCorrSeqPtr> &corrs,
                                     double *globalScale,
                                     Estimator::Opt option);

    template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
    void AddRGBDVelocityFactor(Estimator::Ptr &estimator,
                               const std::string &rgbdTopic,
                               const std::vector<RGBDVelocityCorrPtr> &corrs,
                               Estimator::Opt option);

    void StoreImagesForSfM(const std::string &camTopic,
                           const std::set<std::pair<ns_veta::IndexT, ns_veta::IndexT>> &matchRes);

    ns_veta::Veta::Ptr TryLoadSfMData(const std::string &camTopic,
                                      double errorThd,
                                      std::size_t trackLenThd);

    static void DownsampleVeta(const ns_veta::Veta::Ptr &veta,
                               std::size_t lmNumThd,
                               std::size_t obvNumThd);

    static bool IsRSCamera(const std::string &camTopic);

    static void SaveStageCalibParam(const CalibParamManager::Ptr &par, const std::string &desc);

    static std::vector<OpticalFlowTripleTracePtr> CreateOpticalFlowTraceForRGBD(
        const std::list<RotOnlyVisualOdometer::FeatTrackingInfo> &trackInfoList,
        const std::string &topic);
};

struct CeresDebugCallBack : public ceres::IterationCallback {
private:
    CalibParamManager::Ptr _parMagr;
    const std::string _outputDir;
    std::ofstream _iterInfoFile;
    int _idx;

public:
    explicit CeresDebugCallBack(CalibParamManager::Ptr calibParamManager);

    ~CeresDebugCallBack() override;

    ceres::CallbackReturnType operator()(const ceres::IterationSummary &summary) override;
};

struct CeresViewerCallBack : public ceres::IterationCallback {
private:
    Viewer::Ptr _viewer;

public:
    explicit CeresViewerCallBack(Viewer::Ptr viewer);

    ceres::CallbackReturnType operator()(const ceres::IterationSummary &summary) override;
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_CALIB_SOLVER_H
