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
#include "factor/point_to_surfel_factor.hpp"
#include "solver/batch_opt_option.hpp"
#include "solver/calib_solver.h"
#include "util/utils_tpl.hpp"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::Process() {
    auto outputParams = IsOptionWith(OutputOption::ParamInEachIter, Configor::Preference::Outputs);
    if (outputParams) {
        SaveStageCalibParam(_parMagr, "stage_0_init");
    }
    /**
     * perform the initialization procedure to recover spatiotemporal parameters and
     * interested states in the estimator, such as the splines and the gravity vector
     */

    /* initialize (recover) the rotation spline using raw angular velocity measurements from the
     * gyroscope. If multiple gyroscopes (IMUs) are involved, the extrinsic rotations and time
     * offsets would be also recovered
     */
    this->InitSO3Spline();
    if (outputParams) {
        SaveStageCalibParam(_parMagr, "stage_1_rot_fit");
    }

    /**
     * perform sensor-inertial alignment to recover the gravity vector and extrinsic translations.
     * for each type of sensor, the preparation is first performed to obtain necessary quantities
     * for one-shot sensor-inertial alignment
     */
    this->InitPrepPosCameraInertialAlign();  // visual-inertial (pos scale spline based)
    this->InitPrepVelCameraInertialAlign();  // visual-inertial (vel scale spline based)
    this->InitPrepRGBDInertialAlign();       // rgbd-inertial
    this->InitPrepLiDARInertialAlign();      // lidar-inertial
    this->InitPrepRadarInertialAlign();      // radar-inertial
    this->InitPrepInertialInertialAlign();   // inertial-inertial
    this->InitSensorInertialAlign();         // one-shot sensor-inertial alignment

    if (outputParams) {
        SaveStageCalibParam(_parMagr, "stage_2_align");
    }

    /**
     * recover the linear scale spline using quantities from the one-shot sensor-inertial alignment
     */
    this->InitScaleSpline();

    if (outputParams) {
        SaveStageCalibParam(_parMagr, "stage_3_scale_fit");
    }

    /**
     * this is the end of the initialization procedure, also the start of the batch optimizaion.
     * some preparation operatiors would be performed here
     */
    this->InitPrepBatchOpt();

    /**
     * once the initialization procedure is finished, we print the recovered spatiotemporal
     * parameters to users.
     */
    _parMagr->ShowParamStatus();

    /**
     * perform multi-stage batch optimizations to refine all initialized states to global optimal
     * ones, which is descided by the optimization options from the 'BatchOptOption'
     */
    const int ptsCountInEachScan = Configor::Prior::LiDARDataAssociate::PointToSurfelCountInScan;
    auto options = BatchOptOption::GetOptions();

    for (int i = 0; i < static_cast<int>(options.size()); ++i) {
        spdlog::info("perform '{}-th' batch optimization...", i);
        /**
         * the preparation visualization tasks before the batch optimization.
         */
        _viewer->ClearViewer(Viewer::VIEW_MAP);
        if (Configor::IsRadarIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
            // add radar cloud if radars and pose spline is maintained
            auto color = ns_viewer::Colour::Black().WithAlpha(0.2f);
            _viewer->AddCloud(BuildGlobalMapOfRadar(), Viewer::VIEW_MAP, color, 2.0f);
        }
        std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> lidarPtsCorr;

        /**
         * for the first batch optimization, the lidar global map and frames in the global frames
         * are from the odometry performed in the initialization procedure, to save the computation
         * consumption
         */
        if (i == 0) {
            lidarPtsCorr = DataAssociationForLiDARs(
                // the global lidar map
                _initAsset->globalMap,
                // undistorted frame expressed in the global map
                _initAsset->undistFramesInMap, ptsCountInEachScan);
            _initAsset = nullptr;  // deconstruct data from initialization
        } else {
            auto [curGlobalMap, curUndistFramesInMap] = BuildGlobalMapOfLiDAR();
            lidarPtsCorr = DataAssociationForLiDARs(
                // the global lidar map
                curGlobalMap,
                // undistorted frame expressed in the global map
                curUndistFramesInMap, ptsCountInEachScan);
            // 'curGlobalMap' and 'curUndistFramesInMap' would be deconstructed here
        }
        /**
         * perform batch optimization, association correspondences of cameras, rgbds, and lidars are
         * from addition constructed, while for imus and radars, raw measurements can be directly
         * fused into the estimator.
         * the results would be storaged in '_backup' for by-products output
         */
        _backup = this->BatchOptimization(
            // optimization option
            options.at(i),
            // point to surfel data association for LiDARs
            lidarPtsCorr,
            // visual reprojection data association for cameras
            DataAssociationForPosCameras(),
            // visual velocity creation for rgbd cameras
            DataAssociationForRGBDs(IsOptionWith(OptOption::OPT_VISUAL_DEPTH, options.at(i))),
            // visual velocity creation for optical cameras
            DataAssociationForVelCameras());

        /**
         * update the viewer and output the spatiotemporal parameters after this batch optimization
         * if output is needed, output the stage parameters to the disk
         */
        _viewer->UpdateSplineViewer();
        _parMagr->ShowParamStatus();
        if (outputParams) {
            SaveStageCalibParam(_parMagr, "stage_4_bo_" + std::to_string(i));
        }
    }

/**
 * currently, the data association is performed based on single-sensor data streams, however,
 * cross-model data association can be performed between different sensors, such as point-to-surfel
 * association between (camera, rgbd, and radar) sparse point cloud and the lidar global map.
 * However, this module has not been finished and tested
 */
#define USE_CROSS_MODEL_REFINEMENT 0
#if USE_CROSS_MODEL_REFINEMENT
    for (int i = 0; i < 3; ++i) {
        spdlog::info("perform '{}-th' cross-model batch optimization...", i);
        _viewer->ClearViewer(Viewer::VIEW_MAP);
        // add radar cloud if radars and pose spline is maintained
        if (Configor::IsRadarIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
            auto color = ns_viewer::Colour::Black().WithAlpha(0.2f);
            _viewer->AddCloud(BuildGlobalMapOfRadar(), Viewer::VIEW_MAP, color, 2.0f);
        }

        std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> lidarPtsCorr;
        {
            auto [curGlobalMap, curUndistFramesInMap] = BuildGlobalMapOfLiDAR();
            lidarPtsCorr = DataAssociationForLiDARs(curGlobalMap,
                                                    // frames in map
                                                    curUndistFramesInMap, ptsCountInEachScan);
            // 'curGlobalMap' and 'curUndistFramesInMap' would be deconstructed here
        }
        std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> rgbdPtsCorr;
        {
            auto [rgbdGlobalMap, rgbdScansInGFrame, rgbdScansInLFrame] = BuildGlobalMapOfRGBD();
            rgbdPtsCorr = DataAssociationForRGBDs(rgbdGlobalMap,
                                                  // scans in global frame
                                                  rgbdScansInGFrame,
                                                  // scans in local frame
                                                  rgbdScansInLFrame, ptsCountInEachScan);
            // 'rgbdGlobalMap', 'rgbdScansInGFrame', and 'rgbdScansInLFrame' would be deconstructed
        }
        _backup = this->BatchOptimization(
            // optimization option, here we use the final one
            options.back(),
            // point to surfel data association for LiDARs
            lidarPtsCorr,
            // visual reprojection data association for cameras
            DataAssociationForCameras(),
            // visual velocity creation for rgbd cameras
            {},  // DataAssociationForRGBDs(IsOptionWith(OptOption::OPT_VISUAL_DEPTH,
                 // options.back())),
            // point to surfel data association for RGBDs
            // attention!!! if depth images are not well-matched with rgb images
            // the spatiotemporal parameters between rgb camera and depth camera would be conflict
            rgbdPtsCorr);

        _viewer->UpdateSplineViewer();
        _parMagr->ShowParamStatus();

        if (IsOptionWith(OutputOption::ParamInEachIter, Configor::Preference::Outputs)) {
            SaveStageCalibParam(_parMagr, "stage_5_bo_" + std::to_string(i));
        }
    }
#endif
#undef USE_CROSS_MODEL_REFINEMENT

    /**
     * some tasks after batch optimization
     */
    _viewer->ClearViewer(Viewer::VIEW_MAP);
    if (Configor::IsLiDARIntegrated()) {
        spdlog::info("build final lidar map and point-to-surfel correspondences...");
        // aligned map
        const auto final = BuildGlobalMapOfLiDAR();
        _backup->lidarMap = std::get<0>(final);
        // use large 'ptsCountInEachScan' to keep all point-to-surfel corrs
        // lidar map and corr map would be added to the viewer in this function
        _backup->lidarCorrs =
            DataAssociationForLiDARs(std::get<0>(final), std::get<1>(final), 100000);
    }
    if (Configor::IsRadarIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
        spdlog::info("build final radar map...");
        // radar map would be added to the viewer in this function
        _backup->radarMap = BuildGlobalMapOfRadar();
    }
    if (Configor::IsPosCameraIntegrated()) {
        for (const auto &[topic, sfmData] : _dataMagr->GetSfMData()) {
            _viewer->AddVeta(sfmData, Viewer::VIEW_MAP);
        }
    }
    if (Configor::IsRGBDIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
        // add veta from pixel dynamics
        for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
            const auto &veta = CreateVetaFromOpticalFlow(topic, _backup->ofCorrs.at(topic),
                                                         _parMagr->INTRI.RGBD.at(topic)->intri,
                                                         &CalibSolver::CurDnToW);

            if (veta != nullptr) {
                const auto lenThd = Configor::DataStream::RGBDTopics.at(topic).TrackLengthMin;
                DownsampleVeta(veta, 10000, lenThd);
                // we do not show the pose
                _viewer->AddVeta(veta, Viewer::VIEW_MAP, {}, ns_viewer::Entity::GetUniqueColour());
            }
        }
    }
    if (Configor::IsVelCameraIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
        // add veta from pixel dynamics
        for (const auto &[topic, _] : Configor::DataStream::VelCameraTopics()) {
            const auto &intri = _parMagr->INTRI.Camera.at(topic);

            auto veta = CreateVetaFromOpticalFlow(topic, _backup->ofCorrs.at(topic), intri,
                                                  &CalibSolver::CurCmToW);
            if (veta != nullptr) {
                DownsampleVeta(veta, 10000,
                               Configor::DataStream::CameraTopics.at(topic).TrackLengthMin);
                // we do not show the pose
                _viewer->AddVeta(veta, Viewer::VIEW_MAP, {}, ns_viewer::Entity::GetUniqueColour());
            }
        }
    }

    _solveFinished = true;

    spdlog::info(
        "Solving is finished! Focus on the viewer and press [ctrl+'s'] to save the current scene!");
    spdlog::info("Focus on the viewer and press ['w', 's', 'a', 'd'] to zoom spline viewer!");
}
}  // namespace ns_ikalibr