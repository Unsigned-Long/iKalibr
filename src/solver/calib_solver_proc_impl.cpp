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

#include "solver/calib_solver.h"
#include "solver/batch_opt_option.hpp"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::Process() {
    if (IsOptionWith(OutputOption::ParamInEachIter, Configor::Preference::Outputs)) {
        SaveStageCalibParam(_parMagr, "stage_0_init");
    }
    // --------------
    // initialization
    // --------------
    spdlog::info("initialization...");

    // recover so3 spline
    this->InitSO3Spline();

    // prepare for sensor-alignment alignment
    this->InitPrepCameraInertialAlign();
    this->InitPrepRGBDInertialAlign();
    this->InitPrepLiDARInertialAlign();
    this->InitPrepRadarInertialAlign();

    // sensor-alignment alignment
    this->InitSensorInertialAlign();

    // recover scale spline
    this->InitScaleSpline();

    // prepare for final batch optimizations
    this->InitPrepBatchOpt();

    _parMagr->ShowParamStatus();

    // ------------------
    // batch optimization
    // ------------------
    const int ptsCountInEachScan = Configor::Prior::LiDARDataAssociate::PointToSurfelCountInScan;

    auto options = BatchOptOption::GetOptions();
    for (int i = 0; i < static_cast<int>(options.size()); ++i) {
        spdlog::info("perform '{}-th' batch optimization...", i);
        _viewer->ClearViewer(Viewer::VIEW_MAP);
        // add radar cloud if radars and pose spline is maintained
        if (Configor::IsRadarIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
            auto color = ns_viewer::Colour::Black().WithAlpha(0.2f);
            _viewer->AddCloud(BuildGlobalMapOfRadar(), Viewer::VIEW_MAP, color, 2.0f);
        }
        if (i == 0) {
            // here we store the estimator just for output the hessian matrix
            _backup = this->BatchOptimization(
                // optimization option
                options.at(i),
                // point to surfel data association for LiDARs
                DataAssociationForLiDARs(_initAsset->globalMap, _initAsset->undistFramesInMap,
                                         ptsCountInEachScan),
                // visual reprojection data association for cameras
                DataAssociationForCameras(),
                // visual velocity creation for rgbd cameras
                DataAssociationForRGBDs(
                    IsOptionWith(OptOption::Option::OPT_RGBD_DEPTH, options.at(i))));
            // deconstruct data from initialization
            _initAsset = nullptr;
        } else {
            auto [curGlobalMap, curUndistFramesInMap] = BuildGlobalMapOfLiDAR();
            _backup = this->BatchOptimization(
                // optimization option
                options.at(i),
                // point to surfel data association for LiDARs
                DataAssociationForLiDARs(curGlobalMap, curUndistFramesInMap, ptsCountInEachScan),
                // visual reprojection data association for cameras
                DataAssociationForCameras(),
                // visual velocity creation for rgbd cameras
                DataAssociationForRGBDs(
                    IsOptionWith(OptOption::Option::OPT_RGBD_DEPTH, options.at(i))));
            // 'curGlobalMap' and 'curUndistFramesInMap' would be deconstructed here
        }

        _viewer->UpdateSplineViewer();
        _parMagr->ShowParamStatus();

        if (IsOptionWith(OutputOption::ParamInEachIter, Configor::Preference::Outputs)) {
            SaveStageCalibParam(_parMagr, "stage_4_bo_" + std::to_string(i));
        }
    }
    // backup data and update the viewer
    _viewer->ClearViewer(Viewer::VIEW_MAP);
    if (Configor::IsLiDARIntegrated()) {
        spdlog::info("build final lidar map and point-to-surfel correspondences...");
        // aligned map
        auto final = BuildGlobalMapOfLiDAR();
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
    if (Configor::IsCameraIntegrated()) {
        for (const auto &[topic, sfmData] : _dataMagr->GetSfMData()) {
            _viewer->AddVeta(sfmData, Viewer::VIEW_MAP);
        }
    }
    if (Configor::IsRGBDIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
        // add veta from pixel dynamics
        for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
            auto veta = CreateVetaFromRGBD(topic);
            if (veta != nullptr) {
                DownsampleVeta(veta, 10000,
                               Configor::DataStream::RGBDTopics.at(topic).TrackLengthMin);
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