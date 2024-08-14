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

#include "calib/calib_solver.h"
#include "calib/batch_opt_option.hpp"
#include "core/rotation_estimator.h"
#include "core/lidar_odometer.h"
#include "core/scan_undistortion.h"
#include "core/pts_association.h"
#include "core/rot_only_vo.h"
#include "core/visual_reproj_association.h"
#include "opencv2/highgui.hpp"
#include "util/tqdm.h"
#include "magic_enum_flags.hpp"
#include "core/visual_pixel_dynamic.h"
#include "core/visual_velocity_estimator.h"
#include "core/visual_velocity_sac.h"

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
    auto [globalMap, undistFramesInMap] = Initialization();
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
                DataAssociationForLiDARs(globalMap, undistFramesInMap, ptsCountInEachScan),
                // visual reprojection data association for cameras
                DataAssociationForCameras());
            // deconstruct data
            globalMap.reset(), undistFramesInMap.clear();
        } else {
            auto [curGlobalMap, curUndistFramesInMap] = BuildGlobalMapOfLiDAR();
            _backup = this->BatchOptimization(
                // optimization option
                options.at(i),
                // point to surfel data association for LiDARs
                DataAssociationForLiDARs(curGlobalMap, curUndistFramesInMap, ptsCountInEachScan),
                // visual reprojection data association for cameras
                DataAssociationForCameras());
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

    _solveFinished = true;

    spdlog::info(
        "Solving is finished! Focus on the viewer and press [ctrl+'s'] to save the current scene!");
    spdlog::info("Focus on the viewer and press ['w', 's', 'a', 'd'] to zoom spline viewer!");
}

std::tuple<IKalibrPointCloud::Ptr, std::map<std::string, std::vector<LiDARFrame::Ptr>>>
CalibSolver::Initialization() {
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are
    // poor
    const double st =
        std::max(so3Spline.MinTime(), scaleSpline.MinTime()) + Configor::Prior::TimeOffsetPadding;
    const double et =
        std::min(so3Spline.MaxTime(), scaleSpline.MaxTime()) - Configor::Prior::TimeOffsetPadding;

    // ----------------------------------------------------------------------------------------------------
    // Step 1: initialize (i) so3 spline of the reference IMU, (ii) extrinsic rotations, (iii) time
    // offsets
    // ----------------------------------------------------------------------------------------------------
    spdlog::info("fitting rotation b-spline...");

    // here we  recover the so3 spline first using only the angular velocities from the reference
    // IMU, then estimates other quantities by involving angular velocity measurements from other
    // IMUs. For better readability, we could also optimize them together (not current version)
    auto estimator = Estimator::Create(_splines, _parMagr);
    auto optOption = OptOption::Option::OPT_SO3_SPLINE;
    this->AddGyroFactor(estimator, Configor::DataStream::ReferIMU, optOption);
    auto sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

    if (Configor::DataStream::IMUTopics.size() > 1) {
        estimator = Estimator::Create(_splines, _parMagr);
        if (Configor::Prior::OptTemporalParams) {
            optOption |= OptOption::Option::OPT_TO_BiToBr;
        }
        for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
            this->AddGyroFactor(estimator, topic, optOption);
        }
        // make this problem full rank
        estimator->SetRefIMUParamsConstant();
        // estimator->FixFirSO3ControlPoint();

        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    if (IsOptionWith(OutputOption::ParamInEachIter, Configor::Preference::Outputs)) {
        SaveStageCalibParam(_parMagr, "stage_1_rot_fit");
    }

    // --------------------------------------------------------------------
    // Step 2.1: perform rotation-only visual odometer to recover rotations
    // --------------------------------------------------------------------
    std::shared_ptr<tqdm> bar;
    std::map<std::string, RotOnlyVisualOdometer::Ptr> rotOnlyOdom;
    if (Configor::IsCameraIntegrated()) {
        // how many features to maintain in each image
        constexpr int featNumPerImg = 300;
        // the min distance between two features (to ensure features are distributed uniformly)
        constexpr int minDist = 25;
        for (const auto &[topic, frameVec] : _dataMagr->GetCameraMeasurements()) {
            spdlog::info(
                "perform rotation-only visual odometer to recover extrinsic rotations for '{}'...",
                topic);

            // estimates rotations
            auto odometer = RotOnlyVisualOdometer::Create(featNumPerImg, minDist,
                                                          _parMagr->INTRI.Camera.at(topic));
            // estimates extrinsic rotation between the camera and reference IMU using the estimated
            // rotations
            auto rotEstimator = RotationEstimator::Create();

            bar = std::make_shared<tqdm>();
            for (int i = 0; i < static_cast<int>(frameVec.size()); ++i) {
                bar->progress(i, static_cast<int>(frameVec.size()));

                // if tracking current frame failed, the rotation-only odometer would re-initialize
                if (!odometer->GrabFrame(frameVec.at(i))) {
                    spdlog::warn(
                        "tracking failed when grab the '{}' image frame!!! try to reinitialize", i);
                }

                // we do not want to try to recover the extrinsic rotation too frequent
                if ((odometer->GetRotations().size() < 50) ||
                    (odometer->GetRotations().size() % 5 != 0)) {
                    continue;
                }

                // estimate the extrinsic rotation
                rotEstimator->Estimate(so3Spline, odometer->GetRotations());

                // check solver status
                if (rotEstimator->SolveStatus()) {
                    // assign the estimated extrinsic rotation
                    _parMagr->EXTRI.SO3_CmToBr.at(topic) = rotEstimator->GetSO3SensorToSpline();
                    // once we solve the rotation successfully, break this for loop
                    bar->finish();
                    break;
                }
            }
            if (!rotEstimator->SolveStatus()) {
                throw Status(Status::ERROR,
                             "initialize rotation 'SO3_CmToBr' failed, this may be related to "
                             "insufficiently excited motion or bad images.");
            } else {
                spdlog::info("extrinsic rotation of '{}' is recovered using '{:06}' frames", topic,
                             odometer->GetRotations().size());
            }
            _viewer->UpdateSensorViewer();

            rotOnlyOdom.insert({topic, odometer});
        }
        cv::destroyAllWindows();
    }

    // ---------------------------------------------------------------------------
    // Step 2.2: initialize time offsets of Cameras by hand eye rotation alignment
    // ---------------------------------------------------------------------------
    if (Configor::IsCameraIntegrated() && Configor::Prior::OptTemporalParams) {
        // in last step 2.1, we use a discrete-time rotation-only hand-eye alignment,
        // where only the extrinsic rotations are considered, but here we use a continuous-time
        // rotation-only hand-eye alignment where both extrinsic rotations and temporal parameters
        // are considered.
        spdlog::info("perform rotation alignment to initialize time offset for each camera...");
        estimator = Estimator::Create(_splines, _parMagr);
        optOption = OptOption::Option::OPT_SO3_CmToBr;
        if (Configor::Prior::OptTemporalParams) {
            optOption |= OptOption::Option::OPT_TO_CmToBr;
        }

        // perform time offset estimation and extrinsic rotation refinement
        for (const auto &[topic, _] : Configor::DataStream::CameraTopics) {
            const auto &rotations = rotOnlyOdom.at(topic)->GetRotations();
            // this field should be zero
            double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(topic);
            double weight = Configor::DataStream::CameraTopics.at(topic).Weight;

            for (int i = 0; i < static_cast<int>(rotations.size()) - 1; ++i) {
                const auto &sPose = rotations.at(i), ePose = rotations.at(i + 1);
                // we throw the head and tail data as the rotations from the fitted SO3 Spline in
                // that range are poor
                if (sPose.first + TO_CmToBr < st || ePose.first + TO_CmToBr > et) {
                    continue;
                }

                estimator->AddHandEyeRotationAlignmentForCamera(
                    topic, sPose.first, ePose.first, sPose.second, ePose.second, optOption, weight);
            }
        }

        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    // ---------------------------------
    // Step 2.3: perform SfM for cameras
    // ---------------------------------
    if (Configor::IsCameraIntegrated()) {
        spdlog::info("perform SfM for each camera...");
        int needSfMCount = 0;
        for (const auto &[topic, data] : _dataMagr->GetCameraMeasurements()) {
            // load data if SfM has been performed
            auto veta = TryLoadSfMData(topic,
                                       // for rs camera, as the rs effect is not considered in SfM,
                                       // we relax the landmark selection condition
                                       IsRSCamera(topic) ? 2.0 : 1.0,
                                       Configor::DataStream::CameraTopics.at(topic).TrackLengthMin);
            if (veta != nullptr) {
                spdlog::info("SfM data for camera '{}' is valid!", topic);
                spdlog::info("down sample SfM data for camera '{}'", topic);
                // keep too many landmarks and features in estimator is not always good
                DownsampleVeta(veta, 10000, 10);
                // we store SfM datas in '_dataMagr' as it is the 'calibration data manager'
                _dataMagr->SetSfMData(topic, veta);
                // just for visualization
                _viewer->AddVeta(veta, Viewer::VIEW_MAP);
                spdlog::info(
                    "SfM info for topic '{}' after filtering: view count: {}, landmark count: {}",
                    topic, veta->views.size(), veta->structure.size());
                continue;
            } else {
                spdlog::warn(
                    "SfM data for camera '{}' is not valid! Perform SfM using colmap first!",
                    topic);
            }

            // if the camera is integrated, output images for SfM using thirdparty
            // note that the rotation priors of each frame from the extrinsic rotation and so3
            // spline are utilized in this process to accelerate the feature matching
            auto sfm = VisionOnlySfM::Create(topic, data, _parMagr, so3Spline, _viewer);
            // we use the thirdparty library, i.e., colmap, to solve the SfM problem, rather than
            // ours sfm->PreProcess(); sfm->StructureFromMotion();
            spdlog::info("store images of '{}' for SfM...", topic);
            // output image frames and their corresponding information, which would be re-load after
            // SfM are performed
            StoreImagesForSfM(topic, sfm->FindCovisibility(0.1));
            ++needSfMCount;
        }
        if (needSfMCount != 0) {
            throw ns_ikalibr::Status(Status::FINE,
                                     "images have been output to '{}/images'. "
                                     "Perform SfM for each camera using command lines in file "
                                     "'/sfm_ws/sfm-command-line.txt'!",
                                     ns_ikalibr::Configor::DataStream::OutputPath);
        }
    }

    // -----------------------------------------------------------------------
    // Step 2.4: refine time offsets of Cameras by hand eye rotation alignment
    // -----------------------------------------------------------------------
    if (Configor::IsCameraIntegrated()) {
        // perhaps this is not too necessary. The difference from the last continuous-time
        // rotation-only hand-eye alignment is that we use the rotations recovered by colmap-derived
        // ones, rather than rotation-only odometer-derived ones.
        spdlog::info(
            "perform rotation alignment to refine time offset for each camera using SfM "
            "results...");
        estimator = Estimator::Create(_splines, _parMagr);
        optOption = OptOption::Option::OPT_SO3_CmToBr;
        if (Configor::Prior::OptTemporalParams) {
            optOption |= OptOption::Option::OPT_TO_CmToBr;
        }

        for (const auto &[camTopic, veta] : _dataMagr->GetSfMData()) {
            double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
            double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;

            const auto &frames = _dataMagr->GetCameraMeasurements(camTopic);
            std::vector<ns_ctraj::Posed> constructedFrames;
            constructedFrames.reserve(frames.size());

            // find constructed frames by SfM
            for (const auto &frame : frames) {
                auto viewIter = veta->views.find(frame->GetId());
                if (viewIter == veta->views.cend()) {
                    continue;
                }
                auto poseIter = veta->poses.find(viewIter->second->poseId);
                if (poseIter == veta->poses.cend()) {
                    continue;
                }
                const auto &pose = poseIter->second;
                constructedFrames.emplace_back(pose.Rotation(), pose.Translation(),
                                               frame->GetTimestamp());
            }

            for (int i = 0; i < static_cast<int>(constructedFrames.size()) - 1; ++i) {
                const auto &sPose = constructedFrames.at(i), ePose = constructedFrames.at(i + 1);
                // we throw the head and tail data as the rotations from the fitted SO3 Spline in
                // that range are poor
                if (sPose.timeStamp + TO_CmToBr < st || ePose.timeStamp + TO_CmToBr > et) {
                    continue;
                }

                estimator->AddHandEyeRotationAlignmentForCamera(camTopic, sPose.timeStamp,
                                                                ePose.timeStamp, sPose.so3,
                                                                ePose.so3, optOption, weight);
            }
        }

        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    // -------------------------------------------------------------------------------------
    // Step 3.1: perform rotation-only visual odometer to recover rotations and time offsets
    // -------------------------------------------------------------------------------------
    // 'FeatTrackingInfo' is a list for each rgbd camera, as fail tracking leads to multiple pieces
    std::map<std::string, std::list<RotOnlyVisualOdometer::FeatTrackingInfo>> RGBDTrackingInfo;
    if (Configor::IsRGBDIntegrated()) {
        // how many features to maintain in each image
        constexpr int featNumPerImg = 300;
        // the min distance between two features (to ensure features are distributed uniformly)
        constexpr int minDist = 25;
        for (const auto &[topic, frameVec] : _dataMagr->GetRGBDMeasurements()) {
            spdlog::info(
                "perform rotation-only visual odometer to recover extrinsic rotations for RGBD "
                "camera '{}'...",
                topic);

            // estimates rotations
            auto odometer = RotOnlyVisualOdometer::Create(featNumPerImg, minDist,
                                                          _parMagr->INTRI.RGBD.at(topic)->intri);
            // estimates extrinsic rotation between the rgbd and reference IMU using the estimated
            // rotations
            auto rotEstimator = RotationEstimator::Create();

            bar = std::make_shared<tqdm>();
            for (int i = 0; i < static_cast<int>(frameVec.size()); ++i) {
                bar->progress(i, static_cast<int>(frameVec.size()));

                // if tracking current frame failed, the rotation-only odometer would re-initialize
                if (!odometer->GrabFrame(frameVec.at(i))) {
                    spdlog::warn(
                        "tracking failed when grab the '{}' image frame!!! try to reinitialize", i);
                    // save the tracking information
                    RGBDTrackingInfo[topic].push_back(odometer->GetLmTrackInfo());
                    // clear workspace
                    odometer->ResetWorkspace();
                }

                // we do not want to try to recover the extrinsic rotation too frequent (or has been
                // recovered)
                if (rotEstimator->SolveStatus() || (odometer->GetRotations().size() < 50) ||
                    (odometer->GetRotations().size() % 5 != 0)) {
                    continue;
                }

                // estimate the extrinsic rotation
                rotEstimator->Estimate(so3Spline, odometer->GetRotations());

                // check solver status
                if (rotEstimator->SolveStatus()) {
                    // assign the estimated extrinsic rotation
                    _parMagr->EXTRI.SO3_DnToBr.at(topic) = rotEstimator->GetSO3SensorToSpline();

                    // perform rotation alignment to estimate time offset
                    if (Configor::Prior::OptTemporalParams) {
                        estimator = Estimator::Create(_splines, _parMagr);

                        optOption =
                            OptOption::Option::OPT_SO3_DnToBr | OptOption::Option::OPT_TO_DnToBr;
                        double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(topic);
                        double weight = Configor::DataStream::RGBDTopics.at(topic).Weight;

                        const auto &rotations = odometer->GetRotations();
                        for (int i = 0; i < static_cast<int>(rotations.size()) - 1; ++i) {
                            const auto &sRot = rotations.at(i), eRot = rotations.at(i + 1);
                            // we throw the head and tail data as the rotations from the fitted SO3
                            // Spline in that range are poor
                            if (sRot.first + TO_DnToBr < st || eRot.first + TO_DnToBr > et) {
                                continue;
                            }

                            estimator->AddHandEyeRotationAlignmentForRGBD(
                                topic, sRot.first, eRot.first, sRot.second, eRot.second, optOption,
                                weight);
                        }

                        // we don't want to output the solving information
                        estimator->Solve(Estimator::DefaultSolverOptions(
                                             Configor::Preference::AvailableThreads(), false,
                                             Configor::Preference::UseCudaInSolving),
                                         _priori);
                    }
                    _viewer->UpdateSensorViewer();
                }
            }
            bar->finish();

            // add tracking info
            RGBDTrackingInfo[topic].push_back(odometer->GetLmTrackInfo());

            // check solver status
            if (!rotEstimator->SolveStatus()) {
                throw Status(Status::ERROR,
                             "initialize rotation 'SO3_DnToBr' failed, this may be related to "
                             "insufficiently excited motion or bad images.");
            }
        }
        cv::destroyAllWindows();
    }

    // ------------------------------------------------------------
    // Step 3.2: estimate rgbd-derived body-frame linear velocities
    // ------------------------------------------------------------
    for (const auto &[topic, trackInfoList] : RGBDTrackingInfo) {
        // store
        _dataMagr->SetRGBDPixelDynamics(topic, CreateVisualPixelDynamicForRGBD(trackInfoList));
    }
    // topic, camera frame, body-frame velocity
    std::map<std::string, std::vector<std::pair<CameraFrame::Ptr, Eigen::Vector3d>>>
        rgbdBodyFrameVels;
    for (const auto &[topic, dynamics] : _dataMagr->GetRGBDPixelDynamics()) {
        spdlog::info("estimate RGBD-derived linear velocities for '{}'...", topic);
        // reorganize rgbd-dynamics, store them by frame index
        // camera frame, dynamics in this frame (timestamp, pixel, velocity, depth)
        std::map<CameraFrame::Ptr,
                 std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>>>
            dynamicsInFrame;
        const auto &intri = _parMagr->INTRI.RGBD.at(topic);
        const auto &cameraType = EnumCast::stringToEnum<CameraModelType>(
            Configor::DataStream::RGBDTopics.at(topic).Type);
        const auto &readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);
        for (const auto &dynamic : dynamics) {
            // show the visual pixel dynamic image (tracking features, mid-point pixel velocity)
            // auto img = dynamic->CreatePixelDynamicMat(_parMagr->INTRI.RGBD.at(topic)->intri);
            auto midCamFrame = dynamic->GetMidCameraFrame();
            if (const auto &rgbdVelCorr = dynamic->CreateRGBDVelocityCorr(intri, cameraType);
                rgbdVelCorr->depth > 1E-2 /* 1 cm */) {
                // a valid depth
                dynamicsInFrame[midCamFrame].emplace_back(
                    rgbdVelCorr->MidPoint(), rgbdVelCorr->MidPointVel(readout), rgbdVelCorr->depth);
            }
        }

        // estimate rgbd-derived linear velocities for each frame
        const auto &rgbdIntri = _parMagr->INTRI.RGBD.at(topic);
        const double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(topic);
        const Sophus::SO3d &SO3_DnToBr = _parMagr->EXTRI.SO3_DnToBr.at(topic);
        for (const auto &[frame, dynamics] : dynamicsInFrame) {
            // at least two measurements are required, here we up the ante
            if (dynamics.size() < 5) {
                continue;
            }
            const double timeByBr = frame->GetTimestamp() + TO_DnToBr;
            if (timeByBr < st || timeByBr > et) {
                continue;
            }
            auto res = VisualVelocitySacProblem::VisualVelocityEstimationRANSAC(
                dynamics, rgbdIntri->intri, timeByBr, so3Spline, SO3_DnToBr);
            if (res) {
                rgbdBodyFrameVels[topic].emplace_back(frame, *res);
                // auto img = VisualVelocityEstimator::DrawVisualVelocityMat(
                // dynamics, rgbdIntri->intri, timeByBr, so3Spline, SO3_DnToBr, *res, frame, 0.25);
            }
        }
        // sort timestamps
        auto &curRGBDVels = rgbdBodyFrameVels.at(topic);
        std::sort(curRGBDVels.begin(), curRGBDVels.end(), [](const auto &p1, const auto &p2) {
            return p1.first->GetTimestamp() < p2.first->GetTimestamp();
        });
    }

    // -----------------------------------------------------------------------
    // Step 4.1: initialize extrinsic rotations of LiDARs if they are involved
    // -----------------------------------------------------------------------
    if (Configor::IsLiDARIntegrated()) {
        spdlog::info("LiDARs are integrated, initializing extrinsic rotations of LiDARs...");
        for (const auto &[topic, data] : _dataMagr->GetLiDARMeasurements()) {
            spdlog::info(
                "performing ndt odometer for '{}' for extrinsic rotation initialization...", topic);

            auto lidarOdometer = LiDAROdometer::Create(
                static_cast<float>(Configor::Prior::NDTLiDAROdometer::Resolution),
                Configor::Preference::AvailableThreads());
            auto rotEstimator = RotationEstimator::Create();
            bar = std::make_shared<tqdm>();
            for (int i = 0; i < static_cast<int>(data.size()); ++i) {
                bar->progress(i, static_cast<int>(data.size()));
                // just for visualization
                _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
                _viewer->AddAlignedCloud(data.at(i)->GetScan(), Viewer::VIEW_ASSOCIATION);

                // run the lidar odometer(feed frame to ndt solver)
                lidarOdometer->FeedFrame(data.at(i));

                // we run rotation solver when frame size is 50, 55, 60, ...
                if (lidarOdometer->FrameSize() < 50 || (lidarOdometer->FrameSize() % 5 != 0)) {
                    continue;
                }

                // estimate the rotation
                rotEstimator->Estimate(so3Spline, lidarOdometer->GetOdomPoseVec());

                // check solver status
                if (rotEstimator->SolveStatus()) {
                    // update [rot lidar to imu]
                    _parMagr->EXTRI.SO3_LkToBr.at(topic) = rotEstimator->GetSO3SensorToSpline();
                    // once we solve the rotation successfully, just break
                    bar->finish();
                    break;
                }
            }
            if (!rotEstimator->SolveStatus()) {
                throw Status(Status::ERROR,
                             "initialize rotation 'SO3_LkToBr' failed, this may be related to the "
                             "'NDTResolution' of lidar odometer.");
            } else {
                spdlog::info("extrinsic rotation of '{}' is recovered using '{:06}' frames", topic,
                             lidarOdometer->GetOdomPoseVec().size());
            }
            _viewer->AddCloud(lidarOdometer->GetMap(), Viewer::VIEW_MAP,
                              ns_viewer::Entity::GetUniqueColour(), 2.0f);
            _viewer->UpdateSensorViewer();
        }
        _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
    }

    // --------------------------------------------------
    // Step 4.2: undistort lidar scans and rerun odometer
    // --------------------------------------------------
    std::map<std::string, LiDAROdometer::Ptr> lidarOdometers;
    std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFramesInScan;
    if (Configor::IsLiDARIntegrated()) {
        auto undistHelper = ScanUndistortion::Create(_splines, _parMagr);

        for (const auto &[topic, data] : _dataMagr->GetLiDARMeasurements()) {
            spdlog::info("undistort scans for lidar '{}'...", topic);

            // undistort rotation only using 'UNDIST_SO3' in initialization
            undistFramesInScan[topic] =
                undistHelper->UndistortToScan(data, topic, ScanUndistortion::Option::UNDIST_SO3);
            const auto &undistFrames = undistFramesInScan.at(topic);

            spdlog::info("rerun odometer for lidar '{}' using undistorted scans...", topic);
            lidarOdometers[topic] = LiDAROdometer::Create(
                static_cast<float>(Configor::Prior::NDTLiDAROdometer::Resolution),
                Configor::Preference::AvailableThreads());
            bar = std::make_shared<tqdm>();
            for (int i = 0; i < static_cast<int>(undistFrames.size()); ++i) {
                bar->progress(i, static_cast<int>(undistFrames.size()));

                _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
                _viewer->AddAlignedCloud(data.at(i)->GetScan(), Viewer::VIEW_ASSOCIATION);

                auto curUndistFrame = undistFrames.at(i);
                Eigen::Matrix4d predCurToLast;
                if (i == 0) {
                    predCurToLast = Eigen::Matrix4d::Identity();
                } else {
                    auto lastUndistFrame = undistFrames.at(i - 1);

                    if (curUndistFrame == nullptr || lastUndistFrame == nullptr) {
                        continue;
                    }

                    auto curLtoRef = this->CurLkToW(curUndistFrame->GetTimestamp(), topic);
                    auto lastLtoRef = this->CurLkToW(lastUndistFrame->GetTimestamp(), topic);

                    // if query pose successfully
                    if (curLtoRef && lastLtoRef) {
                        // the translation has not been initialized
                        predCurToLast =
                            ns_ctraj::Posed((*lastLtoRef).so3().inverse() * (*curLtoRef).so3(),
                                            Eigen::Vector3d::Zero())
                                .T();
                    } else {
                        predCurToLast = Eigen::Matrix4d::Identity();
                    }
                }
                lidarOdometers.at(topic)->FeedFrame(curUndistFrame, predCurToLast, i < 100);
            }
            bar->finish();

            _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
            _viewer->AddCloud(lidarOdometers.at(topic)->GetMap(), Viewer::VIEW_MAP,
                              ns_viewer::Entity::GetUniqueColour(), 2.0f);
        }
    }

    // --------------------------------------------------------------------------
    // Step 4.3: initialize time offsets of LiDARs by hand eye rotation alignment
    // --------------------------------------------------------------------------
    if (Configor::IsLiDARIntegrated()) {
        spdlog::info("performing  hand eye rotation alignment for LiDARs...");
        estimator = Estimator::Create(_splines, _parMagr);
        optOption = OptOption::Option::OPT_SO3_LkToBr;
        if (Configor::Prior::OptTemporalParams) {
            optOption |= OptOption::Option::OPT_TO_LkToBr;
        }

        for (const auto &[lidarTopic, odometer] : lidarOdometers) {
            const auto &poseSeq = odometer->GetOdomPoseVec();
            double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
            double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;

            for (int i = 0; i < static_cast<int>(poseSeq.size()) - 1; ++i) {
                const auto &sPose = poseSeq.at(i), ePose = poseSeq.at(i + 1);
                // we throw the head and tail data as the rotations from the fitted SO3 Spline in
                // that range are poor
                if (sPose.timeStamp + TO_LkToBr < st || ePose.timeStamp + TO_LkToBr > et) {
                    continue;
                }

                estimator->AddHandEyeRotationAlignmentForLiDAR(lidarTopic, sPose.timeStamp,
                                                               ePose.timeStamp, sPose.so3,
                                                               ePose.so3, optOption, weight);
            }
        }

        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    // ---------------------------------------------------------------------------
    // Step 5: initialize other spatial parameters using inertial-sensor alignment
    // ---------------------------------------------------------------------------
    spdlog::info("performing inertial alignment to initialize other spatial parameters...");

    // assign the gravity roughly, f = a - g, g = a - f
    Eigen::Vector3d firRefAcce =
        _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU).front()->GetAcce();
    _parMagr->GRAVITY = -firRefAcce.normalized() * Configor::Prior::GravityNorm;
    spdlog::info("rough assigned gravity in world frame: ['{:.3f}', '{:.3f}', '{:.3f}']",
                 _parMagr->GRAVITY(0), _parMagr->GRAVITY(1), _parMagr->GRAVITY(2));

    estimator = Estimator::Create(_splines, _parMagr);

    // we do not optimization the already initialized extrinsic rotations (IMUs', Cameras', and
    // LiDARs') here
    optOption =
        // lidar extrinsic translations
        OptOption::Option::OPT_POS_LkInBr |
        // camera extrinsic translations and visual scale
        OptOption::Option::OPT_POS_CmInBr | OptOption::Option::OPT_VISUAL_GLOBAL_SCALE |
        // radar extrinsics
        OptOption::Option::OPT_SO3_RjToBr | OptOption::Option::OPT_POS_RjInBr |
        // rgbd extrinsics
        OptOption::Option::OPT_POS_DnInBr |
        // imu extrinsic translations and gravity
        OptOption::Option::OPT_POS_BiInBr | OptOption::Option::OPT_GRAVITY;

    // lidar-inertial alignment
    std::map<std::string, std::vector<Eigen::Vector3d>> linVelSeqLk;
    for (const auto &[lidarTopic, odometer] : lidarOdometers) {
        const auto &poseSeq = odometer->GetOdomPoseVec();
        double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);

        // create linear velocity sequence
        linVelSeqLk[lidarTopic] =
            std::vector<Eigen::Vector3d>(poseSeq.size(), Eigen::Vector3d::Zero());
        auto &curLidarLinVelSeq = linVelSeqLk.at(lidarTopic);
        double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;

        const auto &imuFrames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);
        spdlog::info("add lidar-inertial alignment factors for '{}' and '{}'...", lidarTopic,
                     Configor::DataStream::ReferIMU);

        static constexpr int STEP = 5;
        for (int i = 0; i < static_cast<int>(poseSeq.size()) - STEP; ++i) {
            const auto &sPose = poseSeq.at(i), ePose = poseSeq.at(i + STEP);
            // we throw the head and tail data as the rotations from the fitted SO3 Spline in
            // that range are poor
            if (sPose.timeStamp + TO_LkToBr < st || ePose.timeStamp + TO_LkToBr > et) {
                continue;
            }

            if (ePose.timeStamp - sPose.timeStamp < 1E-3 ||
                ePose.timeStamp - sPose.timeStamp > 0.5) {
                continue;
            }

            estimator->AddLiDARInertialAlignment(
                imuFrames, lidarTopic, Configor::DataStream::ReferIMU, sPose, ePose,
                odometer->GetMapTime(), &curLidarLinVelSeq.at(i), &curLidarLinVelSeq.at(i + STEP),
                optOption, weight);
        }
    }

    // visual-inertial alignment
    std::map<std::string, std::vector<Eigen::Vector3d>> linVelSeqCm;
    std::map<std::string, double> visualScaleSeq;
    std::map<std::string, std::vector<ns_ctraj::Posed>> sfmPoseSeq;
    for (const auto &[camTopic, veta] : _dataMagr->GetSfMData()) {
        double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);

        const auto &frames = _dataMagr->GetCameraMeasurements(camTopic);
        // first frame to world
        ns_veta::Posed FirCtoW;
        double firCTime = 0.0;
        for (const auto &frame : frames) {
            auto viewIter = veta->views.find(frame->GetId());
            if (viewIter == veta->views.cend()) {
                continue;
            }
            auto poseIter = veta->poses.find(viewIter->second->poseId);
            if (poseIter == veta->poses.cend()) {
                continue;
            }
            FirCtoW = poseIter->second, firCTime = frame->GetTimestamp();
            break;
        }
        // set first valid camera frame as world frame (transform from world to first frame)
        PerformTransformForVeta(veta, FirCtoW.Inverse(), 1.0);

        // load poses
        sfmPoseSeq[camTopic] = {};
        auto &constructedFrames = sfmPoseSeq.at(camTopic);
        constructedFrames.reserve(frames.size());
        for (const auto &frame : frames) {
            auto viewIter = veta->views.find(frame->GetId());
            if (viewIter == veta->views.cend()) {
                continue;
            }
            auto poseIter = veta->poses.find(viewIter->second->poseId);
            if (poseIter == veta->poses.cend()) {
                continue;
            }
            const auto &pose = poseIter->second;
            constructedFrames.emplace_back(pose.Rotation(), pose.Translation(),
                                           frame->GetTimestamp());
        }

        // create linear velocity sequence
        linVelSeqCm[camTopic] =
            std::vector<Eigen::Vector3d>(constructedFrames.size(), Eigen::Vector3d::Zero());
        auto &curCamLinVelSeq = linVelSeqCm.at(camTopic);
        double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;
        // create scale
        visualScaleSeq[camTopic] = 1.0;
        auto &scale = visualScaleSeq.at(camTopic);

        const auto &imuFrames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);
        spdlog::info("add visual-inertial alignment factors for '{}' and '{}'...", camTopic,
                     Configor::DataStream::ReferIMU);

        static constexpr int STEP = 5;
        for (int i = 0; i < static_cast<int>(constructedFrames.size()) - STEP; ++i) {
            const auto &sPose = constructedFrames.at(i), ePose = constructedFrames.at(i + STEP);
            // we throw the head and tail data as the rotations from the fitted SO3 Spline in
            // that range are poor
            if (sPose.timeStamp + TO_CmToBr < st || ePose.timeStamp + TO_CmToBr > et) {
                continue;
            }

            if (ePose.timeStamp - sPose.timeStamp < 1E-3 ||
                ePose.timeStamp - sPose.timeStamp > 0.5) {
                continue;
            }

            estimator->AddVisualInertialAlignment(
                imuFrames, camTopic, Configor::DataStream::ReferIMU, sPose, ePose, firCTime,
                &curCamLinVelSeq.at(i), &curCamLinVelSeq.at(i + STEP), &scale, optOption, weight);
        }
    }

    // inertial alignment (only when more than or equal to 2 num IMUs are involved)
    constexpr double dt = 0.1;
    std::vector<Eigen::Vector3d> linVelSeqBr(std::floor((et - st) / dt), Eigen::Vector3d::Zero());
    if (Configor::DataStream::IMUTopics.size() >= 2) {
        for (const auto &[topic, frames] : _dataMagr->GetIMUMeasurements()) {
            spdlog::info("add inertial alignment factors for '{}'...", topic);

            for (int i = 0; i < static_cast<int>(linVelSeqBr.size()) - 1; ++i) {
                int sIdx = i, eIdx = i + 1;
                double sTimeByBr = sIdx * dt + st, eTimeByBr = eIdx * dt + st;
                Eigen::Vector3d *sVel = &linVelSeqBr.at(sIdx), *eVel = &linVelSeqBr.at(eIdx);

                estimator->AddInertialAlignment(
                    frames, topic, sTimeByBr, eTimeByBr, sVel, eVel, optOption,
                    Configor::DataStream::IMUTopics.at(topic).AcceWeight);
            }
        }
    }

    // radar-inertial alignment
    for (const auto &[radarTopic, radarMes] : _dataMagr->GetRadarMeasurements()) {
        double weight = Configor::DataStream::RadarTopics.at(radarTopic).Weight;
        double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(radarTopic);

        const auto &frames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);
        spdlog::info("add radar-inertial alignment factors for '{}' and '{}'...", radarTopic,
                     Configor::DataStream::ReferIMU);

        for (int i = 0; i < static_cast<int>(radarMes.size()) - 1; ++i) {
            const auto &sArray = radarMes.at(i), eArray = radarMes.at(i + 1);

            // spdlog::info("sAry count: {}, eAry count: {}",
            //              sArray->GetTargets().size(), eArray->GetTargets().size());

            // to estimate the radar velocity by linear least-squares solver
            // the minim targets number required is 3
            if (sArray->GetTargets().size() < 10 || eArray->GetTargets().size() < 10) {
                continue;
            }

            if (sArray->GetTimestamp() + TO_RjToBr < st ||
                eArray->GetTimestamp() + TO_RjToBr > et) {
                continue;
            }

            if (eArray->GetTimestamp() - sArray->GetTimestamp() < 1E-3 ||
                eArray->GetTimestamp() - sArray->GetTimestamp() > 0.5) {
                continue;
            }

            estimator->AddRadarInertialRotRoughAlignment(frames, Configor::DataStream::ReferIMU,
                                                         radarTopic, sArray, eArray, optOption,
                                                         weight);
        }
    }

    // rgbd-inertial alignment
    for (const auto &[rgbdTopic, bodyFrameVels] : rgbdBodyFrameVels) {
        double weight = Configor::DataStream::RGBDTopics.at(rgbdTopic).Weight;
        double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(rgbdTopic);

        const auto &frames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);
        spdlog::info("add rgbd-inertial alignment factors for '{}' and '{}'...", rgbdTopic,
                     Configor::DataStream::ReferIMU);

        for (int i = 0; i < static_cast<int>(bodyFrameVels.size()) - 1; ++i) {
            const auto &[sFrame, sVel] = bodyFrameVels.at(i);
            const auto &[eFrame, eVel] = bodyFrameVels.at(i + 1);

            if (sFrame->GetTimestamp() + TO_DnToBr < st ||
                eFrame->GetTimestamp() + TO_DnToBr > et) {
                continue;
            }

            if (eFrame->GetTimestamp() - sFrame->GetTimestamp() < 1E-3 ||
                eFrame->GetTimestamp() - sFrame->GetTimestamp() > 0.5) {
                continue;
            }

            // we don't consider the readout time here (i.e., the camera frame time is treated as
            // the velocity of the midpoint)
            estimator->AddRGBDInertialAlignment(frames, Configor::DataStream::ReferIMU, rgbdTopic,
                                                bodyFrameVels.at(i), bodyFrameVels.at(i + 1),
                                                optOption, weight);
        }
    }

    // fix spatiotemporal parameters of reference sensor
    // make this problem full rank
    estimator->SetRefIMUParamsConstant();

    sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

    if (Configor::IsRadarIntegrated()) {
        estimator = Estimator::Create(_splines, _parMagr);
        // radar-inertial alignment
        for (const auto &[radarTopic, radarMes] : _dataMagr->GetRadarMeasurements()) {
            double weight = Configor::DataStream::RadarTopics.at(radarTopic).Weight;
            double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(radarTopic);

            const auto &frames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);
            spdlog::info("add radar-inertial alignment factors for '{}' and '{}'...", radarTopic,
                         Configor::DataStream::ReferIMU);

            for (int i = 0; i < static_cast<int>(radarMes.size()) - 1; ++i) {
                const auto &sArray = radarMes.at(i), eArray = radarMes.at(i + 1);

                // spdlog::info("sAry count: {}, eAry count: {}",
                //              sArray->GetTargets().size(), eArray->GetTargets().size());

                // to estimate the radar velocity by linear least-squares solver
                // the minim targets number required is 3
                if (sArray->GetTargets().size() < 10 || eArray->GetTargets().size() < 10) {
                    continue;
                }

                if (sArray->GetTimestamp() + TO_RjToBr < st ||
                    eArray->GetTimestamp() + TO_RjToBr > et) {
                    continue;
                }

                estimator->AddRadarInertialAlignment(frames, Configor::DataStream::ReferIMU,
                                                     radarTopic, sArray, eArray, optOption, weight);
            }
        }
        estimator->SetRefIMUParamsConstant();
        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    if (IsOptionWith(OutputOption::ParamInEachIter, Configor::Preference::Outputs)) {
        SaveStageCalibParam(_parMagr, "stage_2_align");
    }

    for (const auto &[topic, lidarOdom] : lidarOdometers) {
        _viewer->AddCloud(lidarOdom->GetMap(), Viewer::VIEW_MAP,
                          ns_viewer::Entity::GetUniqueColour(), 2.0f);
    }
    // recover scale fro veta
    for (const auto &[camTopic, veta] : _dataMagr->GetSfMData()) {
        spdlog::info("visual global scale for camera '{}': {:.3f}", camTopic,
                     visualScaleSeq.at(camTopic));
        PerformTransformForVeta(veta, ns_veta::Posed(), visualScaleSeq.at(camTopic));
        _viewer->AddVeta(veta, Viewer::VIEW_MAP);
    }
    // perform scale for 'sfmPoseSeq', which would used for scale spline recovery
    for (auto &[camTopic, poseSeq] : sfmPoseSeq) {
        const double scale = visualScaleSeq.at(camTopic);
        for (auto &pose : poseSeq) {
            pose.t *= scale;
        }
    }

    // deconstruct data
    linVelSeqLk.clear(), linVelSeqBr.clear(), linVelSeqCm.clear(), visualScaleSeq.clear();

    // ----------------------------------------------------------------------------
    // Step 6: recover scale spline (linear acceleration, velocity, or translation)
    // ----------------------------------------------------------------------------
    spdlog::info("performing scale spline recovery...");

    estimator = Estimator::Create(_splines, _parMagr);
    optOption = OptOption::Option::OPT_SCALE_SPLINE;

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
            optOption |= OptOption::Option::OPT_SO3_SPLINE;
            // radar-derived velocities
            for (const auto &[topic, data] : _dataMagr->GetRadarMeasurements()) {
                const double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(topic);
                const auto &SO3_RjToBr = _parMagr->EXTRI.SO3_RjToBr.at(topic);
                const Eigen::Vector3d &POS_RjInBr = _parMagr->EXTRI.POS_RjInBr.at(topic);
                const double weight = Configor::DataStream::RadarTopics.at(topic).Weight;

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

                    estimator->AddLinearScaleConstraint<VelDeriv>(timeByBr, LIN_VEL_BrToWInW,
                                                                  optOption, weight);
                }
            }

            // rgbd-derived velocities
            for (const auto &[rgbdTopic, bodyFrameVels] : rgbdBodyFrameVels) {
                double weight = Configor::DataStream::RGBDTopics.at(rgbdTopic).Weight;
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

                    estimator->AddLinearScaleConstraint<VelDeriv>(timeByBr, LIN_VEL_BrToWInW,
                                                                  optOption, weight);
                }
            }

            // add acceleration factor using inertial measurements only from reference IMU
            this->AddAcceFactor<TimeDeriv::LIN_VEL_SPLINE>(
                estimator, Configor::DataStream::ReferIMU, optOption);
            this->AddGyroFactor(estimator, Configor::DataStream::ReferIMU, optOption);
            // if optimize time offsets, we first recover the scale spline, then construct a ls
            // problem to recover time offsets, only for radars
            if (Configor::IsRadarIntegrated() && Configor::Prior::OptTemporalParams) {
                // we don't want to output the solving information
                estimator->Solve(
                    Estimator::DefaultSolverOptions(Configor::Preference::AvailableThreads(), false,
                                                    Configor::Preference::UseCudaInSolving));

                estimator = Estimator::Create(_splines, _parMagr);
                optOption = OptOption::Option::OPT_TO_RjToBr;
                for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
                    this->AddRadarFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
                }
            }

        } break;
        case TimeDeriv::LIN_POS_SPLINE: {
            CalibSolver::SplineBundleType::Ptr roughSplines;
            // this value equals to zero
            constexpr int PosDeriv =
                TimeDeriv::Deriv<TimeDeriv::LIN_POS_SPLINE, TimeDeriv::LIN_POS>();
            optOption |= OptOption::Option::OPT_SO3_SPLINE;

            std::vector<double> headTimeVec, tailTimeVec;
            for (const auto &[lidarTopic, odometer] : lidarOdometers) {
                const double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
                headTimeVec.push_back(odometer->GetOdomPoseVec().front().timeStamp + TO_LkToBr);
                tailTimeVec.push_back(odometer->GetOdomPoseVec().back().timeStamp + TO_LkToBr);
            }
            for (const auto &[camTopic, poseSeq] : sfmPoseSeq) {
                const double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                headTimeVec.push_back(poseSeq.front().timeStamp + TO_CmToBr);
                tailTimeVec.push_back(poseSeq.back().timeStamp + TO_CmToBr);
            }

            const double minTime = *std::min_element(headTimeVec.begin(), headTimeVec.end()) - 1E-3;
            const double maxTime = *std::max_element(tailTimeVec.begin(), tailTimeVec.end()) + 1E-3;

            // create rough splines (the time distance is larger than that from configure as our
            // poses are not dense)
            double dtRoughSum = 0.0;
            int count = 0;
            if (Configor::IsLiDARIntegrated()) {
                dtRoughSum += 2.0 / _dataMagr->GetLiDARAvgFrequency();
                ++count;
            }
            if (Configor::IsCameraIntegrated()) {
                dtRoughSum += 2.0 / _dataMagr->GetCameraAvgFrequency();
                ++count;
            }
            const double dtRough = dtRoughSum / count;
            roughSplines = CreateSplineBundle(minTime, maxTime, dtRough, dtRough);

            estimator = Estimator::Create(roughSplines, nullptr);

            for (const auto &[lidarTopic, odometer] : lidarOdometers) {
                double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
                double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;
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
                    estimator->AddSO3Constraint(SE3_LkToLk0.timeStamp + TO_LkToBr,
                                                SE3_BrToBr0.so3(), optOption, weight);
                    estimator->AddLinearScaleConstraint<PosDeriv>(SE3_LkToLk0.timeStamp + TO_LkToBr,
                                                                  SE3_BrToBr0.translation(),
                                                                  optOption, weight);
                }
            }

            for (const auto &[camTopic, poseSeq] : sfmPoseSeq) {
                double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;
                auto SE3_BrToCm = _parMagr->EXTRI.SE3_CmToBr(camTopic).inverse();

                // the translation part is not exactly rigorous as the translation spline is not
                // recovered yet this assumes that the origin of map frame of {Cm} is the same as
                // that of the world frame {br0}
                auto SE3_Cm0ToBr0 = this->CurCmToW(poseSeq.front().timeStamp, camTopic);

                if (SE3_Cm0ToBr0 == std::nullopt) {
                    throw Status(Status::CRITICAL,
                                 "map time of '{}' is out of time range of splines!", camTopic);
                }

                for (const auto &SE3_CmToCm0 : poseSeq) {
                    auto SE3_BrToBr0 = *SE3_Cm0ToBr0 * SE3_CmToCm0.se3() * SE3_BrToCm;
                    estimator->AddSO3Constraint(SE3_CmToCm0.timeStamp + TO_CmToBr,
                                                SE3_BrToBr0.so3(), optOption, weight);
                    estimator->AddLinearScaleConstraint<PosDeriv>(SE3_CmToCm0.timeStamp + TO_CmToBr,
                                                                  SE3_BrToBr0.translation(),
                                                                  optOption, weight);
                }
            }

            // add tail factors (constraints) to maintain enough observability
            estimator->AddLinScaleTailConstraint(optOption, 1.0);
            estimator->AddSO3TailConstraint(optOption, 1.0);

            // we don't want to output the solving information
            estimator->Solve(
                Estimator::DefaultSolverOptions(Configor::Preference::AvailableThreads(), false,
                                                Configor::Preference::UseCudaInSolving),
                _priori);

            const auto &rSo3Spline = roughSplines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
            const auto &rScaleSpline =
                roughSplines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

            estimator = Estimator::Create(_splines, _parMagr);

            for (double t = minTime; t < maxTime;) {
                estimator->AddSO3Constraint(t, rSo3Spline.Evaluate(t), optOption, 1.0);
                estimator->AddLinearScaleConstraint<PosDeriv>(t, rScaleSpline.Evaluate(t),
                                                              optOption, 1.0);
                t += 0.01;
            }
        } break;
    }
    sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

    if (GetScaleType() == TimeDeriv::LIN_POS_SPLINE && Configor::IsRadarIntegrated() &&
        Configor::Prior::OptTemporalParams) {
        // in this case, the time offsets of radars have not been recovered
        estimator = Estimator::Create(_splines, _parMagr);
        for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
            this->AddRadarFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic,
                                                            OptOption::Option::OPT_TO_RjToBr);
        }
        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    if (IsOptionWith(OutputOption::ParamInEachIter, Configor::Preference::Outputs)) {
        SaveStageCalibParam(_parMagr, "stage_3_scale_fit");
    }

    _parMagr->ShowParamStatus();
    IKALIBR_DEBUG
    std::cin.get();

    //------------------------------------------------------
    // Step 7: align initialized states to gravity direction
    //------------------------------------------------------
    spdlog::info("aligning all states to gravity direction...");
    AlignStatesToGravity();
    _viewer->UpdateSplineViewer();

    // --------------------------------------------------------------
    // Step 8.1 : trans veta to world frame if Cameras are integrated
    // --------------------------------------------------------------
    for (const auto &[camTopic, poseSeq] : sfmPoseSeq) {
        auto SE3_Cm0ToBr0 = this->CurCmToW(poseSeq.front().timeStamp, camTopic);

        if (SE3_Cm0ToBr0 == std::nullopt) {
            throw Status(Status::CRITICAL, "map time of '{}' is out of time range of splines!",
                         camTopic);
        }
        PerformTransformForVeta(_dataMagr->GetSfMData(camTopic),
                                ns_veta::Posed(SE3_Cm0ToBr0->so3(), SE3_Cm0ToBr0->translation()),
                                1.0);
    }

    // ----------------------------------------------------------------
    // Step 8.2 : build map and undisto frames if LiDARs are integrated
    // ----------------------------------------------------------------
    spdlog::info("build global map and undisto lidar frames in world...");
    IKalibrPointCloud::Ptr globalMap(new IKalibrPointCloud);
    std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFramesInMap;
    for (const auto &[lidarTopic, odometer] : lidarOdometers) {
        auto SE3_Lk0ToBr0 = this->CurLkToW(odometer->GetMapTime(), lidarTopic);

        if (SE3_Lk0ToBr0 == std::nullopt) {
            throw Status(Status::CRITICAL, "map time of '{}' is out of time range of splines!",
                         lidarTopic);
        }

        const auto &curUndistFramesInScan = undistFramesInScan.at(lidarTopic);

        undistFramesInMap[lidarTopic] = std::vector<LiDARFrame::Ptr>(curUndistFramesInScan.size());
        auto &curUndistFramesInMap = undistFramesInMap.at(lidarTopic);

        const auto &poseSeq = odometer->GetOdomPoseVec();

        for (int i = 0; i < static_cast<int>(poseSeq.size()); ++i) {
            const auto &undistoScan = curUndistFramesInScan.at(i);
            if (undistoScan == nullptr) {
                curUndistFramesInMap.at(i) = nullptr;
            } else {
                const auto &SE3_ScanToLk0 = poseSeq.at(i);
                Sophus::SE3d SE3_ScanToBr0 = *SE3_Lk0ToBr0 * SE3_ScanToLk0.se3();
                IKalibrPointCloud::Ptr scanInBr0(new IKalibrPointCloud);
                pcl::transformPointCloud(*undistoScan->GetScan(), *scanInBr0,
                                         SE3_ScanToBr0.matrix().cast<float>());
                curUndistFramesInMap.at(i) =
                    LiDARFrame::Create(undistoScan->GetTimestamp(), scanInBr0);

                *globalMap += *scanInBr0;
            }
        }
    }
    return {globalMap, undistFramesInMap};
}

std::tuple<IKalibrPointCloud::Ptr, std::map<std::string, std::vector<LiDARFrame::Ptr>>>
CalibSolver::BuildGlobalMapOfLiDAR() {
    if (!Configor::IsLiDARIntegrated()) {
        return {};
    }

    // -------------------------------------------
    // undisto frames and build map in world frame
    // -------------------------------------------
    auto undistHelper = ScanUndistortion::Create(_splines, _parMagr);

    std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFrames;
    IKalibrPointCloud::Ptr mapCloud(new IKalibrPointCloud);
    for (const auto &[topic, data] : _dataMagr->GetLiDARMeasurements()) {
        spdlog::info("undistort scans for lidar '{}'...", topic);
        undistFrames[topic] =
            undistHelper->UndistortToRef(data, topic, ScanUndistortion::Option::ALL);

        spdlog::info("marge scans from lidar '{}' to map...", topic);
        for (auto &frame : undistFrames.at(topic)) {
            if (frame == nullptr) {
                continue;
            }
            *mapCloud += *frame->GetScan();
        }
    }

    IKalibrPointCloud::Ptr newMapCloud(new IKalibrPointCloud);
    std::vector<int> index;
    pcl::removeNaNFromPointCloud(*mapCloud, *newMapCloud, index);

    return {newMapCloud, undistFrames};
}

IKalibrPointCloud::Ptr CalibSolver::BuildGlobalMapOfRadar() {
    if (!Configor::IsRadarIntegrated() || GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
        return {};
    }

    IKalibrPointCloud::Ptr radarCloud(new IKalibrPointCloud);

    for (const auto &[topic, data] : _dataMagr->GetRadarMeasurements()) {
        const double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(topic);
        IKalibrPointCloud::Ptr curRadarCloud(new IKalibrPointCloud);
        curRadarCloud->reserve(data.size() * data.front()->GetTargets().size());
        for (const auto &ary : data) {
            for (const auto &frame : ary->GetTargets()) {
                auto SE3_CurRjToW = CurRjToW(frame->GetTimestamp(), topic);
                if (SE3_CurRjToW == std::nullopt) {
                    continue;
                }
                Eigen::Vector3d p = *SE3_CurRjToW * frame->GetTargetXYZ();
                IKalibrPoint p2;
                p2.timestamp = frame->GetTimestamp() + TO_RjToBr;
                p2.x = static_cast<float>(p(0));
                p2.y = static_cast<float>(p(1));
                p2.z = static_cast<float>(p(2));
                curRadarCloud->push_back(p2);
            }
        }
        *radarCloud += *curRadarCloud;
    }

    // ---------------------------
    // down sample the radar cloud
    // ---------------------------
    pcl::VoxelGrid<IKalibrPoint> filter;
    filter.setInputCloud(radarCloud);
    auto size = static_cast<float>(Configor::Prior::LiDARDataAssociate::MapDownSample * 2.0);
    filter.setLeafSize(size, size, size);

    IKalibrPointCloud::Ptr radarCloudSampled(new IKalibrPointCloud);
    filter.filter(*radarCloudSampled);
    _viewer->AddStarMarkCloud(radarCloudSampled, Viewer::VIEW_MAP);

    return radarCloud;
}

std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> CalibSolver::DataAssociationForLiDARs(
    const IKalibrPointCloud::Ptr &map,
    const std::map<std::string, std::vector<LiDARFrame::Ptr>> &undistFrames,
    int ptsCountInEachScan) {
    if (!Configor::IsLiDARIntegrated()) {
        return {};
    }

    // ---------------------------------
    // Step 1: down sample the map cloud
    // ---------------------------------
    pcl::VoxelGrid<IKalibrPoint> filter;
    filter.setInputCloud(map);
    auto size = static_cast<float>(Configor::Prior::LiDARDataAssociate::MapDownSample);
    filter.setLeafSize(size, size, size);

    IKalibrPointCloud::Ptr mapDownSampled(new IKalibrPointCloud);
    filter.filter(*mapDownSampled);

    _viewer->AddAlignedCloud(mapDownSampled, Viewer::VIEW_MAP, -_parMagr->GRAVITY.cast<float>(),
                             2.0f);

    // ------------------------------------------------
    // Step 2: perform data association for each frames
    // ------------------------------------------------
    auto associator = PointToSurfelAssociator::Create(
        // we use the dense map to create data associator for high-perform point-to-surfel search
        map, Configor::Prior::LiDARDataAssociate::MapResolution,
        Configor::Prior::LiDARDataAssociate::MapDepthLevels);
    auto condition = PointToSurfelCondition();
    _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
    _viewer->AddSurfelMap(associator->GetSurfelMap(), condition, Viewer::VIEW_ASSOCIATION);
    _viewer->AddCloud(mapDownSampled, Viewer::VIEW_ASSOCIATION,
                      ns_viewer::Colour::Black().WithAlpha(0.2f), 2.0f);

    // deconstruction
    mapDownSampled.reset();

    std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> pointToSurfel;

    std::size_t count = 0;
    std::shared_ptr<tqdm> bar;
    for (const auto &[topic, framesInMap] : undistFrames) {
        const auto &rawFrames = _dataMagr->GetLiDARMeasurements(topic);
        spdlog::info("perform point to surfel association for lidar '{}'...", topic);

        // for each scan, we keep 'ptsCountInEachScan' point to surfel corrs
        pointToSurfel[topic] = {};
        auto &curPointToSurfel = pointToSurfel.at(topic);
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(framesInMap.size()); ++i) {
            bar->progress(i, static_cast<int>(framesInMap.size()));

            if (framesInMap.at(i) == nullptr || rawFrames.at(i) == nullptr) {
                continue;
            }

            auto ptsVec = associator->Association(framesInMap.at(i)->GetScan(),
                                                  rawFrames.at(i)->GetScan(), condition);

            curPointToSurfel.insert(curPointToSurfel.end(), ptsVec.cbegin(), ptsVec.cend());
        }
        bar->finish();

        // downsample
        int expectCount = ptsCountInEachScan * static_cast<int>(rawFrames.size());
        if (static_cast<int>(curPointToSurfel.size()) > expectCount) {
            std::map<ufo::map::Node, std::vector<PointToSurfelCorr::Ptr>> nodes;
            for (const auto &corr : curPointToSurfel) {
                nodes[corr->node].push_back(corr);
            }
            std::size_t numEachNode = expectCount / nodes.size() + 1;
            curPointToSurfel.clear();

            // uniform sampling
            std::default_random_engine engine(
                std::chrono::steady_clock::now().time_since_epoch().count());
            for (const auto &[node, corrs] : nodes) {
                auto newCorrs =
                    SamplingWoutReplace2(engine, corrs, std::min(corrs.size(), numEachNode));
                curPointToSurfel.insert(curPointToSurfel.end(), newCorrs.cbegin(), newCorrs.cend());
            }
        }
        count += curPointToSurfel.size();
    }
    spdlog::info("total point to surfel count for LiDARs: {}", count);
    _viewer->AddPointToSurfel(associator->GetSurfelMap(), pointToSurfel, Viewer::VIEW_ASSOCIATION);

    return pointToSurfel;
}

std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>>
CalibSolver::DataAssociationForCameras() {
    if (!Configor::IsCameraIntegrated()) {
        return {};
    }

    std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>> corrs;
    for (const auto &[topic, sfmData] : _dataMagr->GetSfMData()) {
        spdlog::info("performing visual reprojection data association for camera '{}'...", topic);
        corrs[topic] =
            VisualReProjAssociator::Create(EnumCast::stringToEnum<CameraModelType>(
                                               Configor::DataStream::CameraTopics.at(topic).Type))
                ->Association(*sfmData, _parMagr->INTRI.Camera.at(topic));
        _viewer->AddVeta(sfmData, Viewer::VIEW_MAP);
        spdlog::info("visual reprojection sequences for '{}': {}", topic, corrs.at(topic).size());
    }
    return corrs;
}

std::vector<VisualPixelDynamic::Ptr> CalibSolver::CreateVisualPixelDynamicForRGBD(
    const std::list<RotOnlyVisualOdometer::FeatTrackingInfo> &trackInfoList) {
    std::vector<VisualPixelDynamic::Ptr> dynamics;
    for (const auto &trackInfo : trackInfoList) {
        for (const auto &[id, info] : trackInfo) {
            // for each tracked feature
            if (info.size() < 3) {
                continue;
            }
            // store in groups of three
            auto groupSize = info.size() - 2;
            auto iter0 = info.cbegin();
            for (int i = 0; i < static_cast<int>(groupSize); ++i, ++iter0) {
                auto iter1 = iter0, iter2 = iter0;
                std::advance(iter1, 1), std::advance(iter2, 2);
                std::array<std::pair<CameraFramePtr, Eigen::Vector2d>, 3> movement;
                movement.at(0) = {
                    iter0->first,
                    Eigen::Vector2d(iter0->second.undistorted.x, iter0->second.undistorted.y),
                };
                movement.at(1) = {
                    iter1->first,
                    Eigen::Vector2d(iter1->second.undistorted.x, iter1->second.undistorted.y),
                };
                movement.at(2) = {
                    iter2->first,
                    Eigen::Vector2d(iter2->second.undistorted.x, iter2->second.undistorted.y),
                };
                dynamics.push_back(VisualPixelDynamic::Create(movement));
            }
        }
    }
    return dynamics;
}

CalibSolver::BackUp::Ptr CalibSolver::BatchOptimization(
    OptOption::Option optOption,
    const std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> &ptsCorrs,
    const std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>> &visualCorrs) {
    auto GetOptString = [](OptOption::Option opt) -> std::string {
        std::stringstream stringStream;
        stringStream << magic_enum::enum_flags_name(opt);
        return stringStream.str();
    };
    spdlog::info("Optimization option: {}", GetOptString(optOption));

    auto estimator = Estimator::Create(_splines, _parMagr);
    auto visualGlobalScale = std::make_shared<double>(1.0);

    switch (GetScaleType()) {
        case TimeDeriv::LIN_ACCE_SPLINE: {
            for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
                this->AddAcceFactor<TimeDeriv::LIN_ACCE_SPLINE>(estimator, topic, optOption);
                this->AddGyroFactor(estimator, topic, optOption);
            }
        } break;
        case TimeDeriv::LIN_VEL_SPLINE: {
            for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
                this->AddRadarFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
            }
            for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
                this->AddAcceFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
                this->AddGyroFactor(estimator, topic, optOption);
            }
        } break;
        case TimeDeriv::LIN_POS_SPLINE: {
            for (const auto &[lidarTopic, corrSeqVec] : ptsCorrs) {
                this->AddPointToSurfelFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, lidarTopic,
                                                                        corrSeqVec, optOption);
            }
            for (const auto &[camTopic, corrSeqVec] : visualCorrs) {
                OptOption::Option visualOpt = optOption;
                if (IsOptionWith(OptOption::Option::OPT_RS_CAM_READOUT_TIME, visualOpt)) {
                    if (IsRSCamera(camTopic)) {
                        spdlog::info(
                            "camera '{}' is a rolling shutter (RS) camera, "
                            "use optimization option 'OPT_RS_CAM_READOUT_TIME'",
                            camTopic);
                    } else {
                        visualOpt ^= OptOption::Option::OPT_RS_CAM_READOUT_TIME;
                        spdlog::info(
                            "camera '{}' is a global shutter (GS) camera, "
                            "remove optimization option 'OPT_RS_CAM_READOUT_TIME'",
                            camTopic);
                    }
                }
                this->AddVisualReprojectionFactor<TimeDeriv::LIN_POS_SPLINE>(
                    estimator, camTopic, corrSeqVec, visualGlobalScale.get(), visualOpt);
            }
            for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
                this->AddRadarFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic, optOption);
            }
            for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
                this->AddAcceFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic, optOption);
                this->AddGyroFactor(estimator, topic, optOption);
            }
        } break;
    }

    // make this problem full rank
    estimator->SetRefIMUParamsConstant();
    // estimator->FixFirSO3ControlPoint();

    auto sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

    // align states to the gravity
    AlignStatesToGravity();

    // for better map consistency in visualization, we update the veta every time
    for (const auto &[topic, reprojCorrVec] : visualCorrs) {
        auto &veta = _dataMagr->GetSfMData(topic);
        auto &intri = _parMagr->INTRI.Camera.at(topic);
        // compute pose
        for (const auto &[viewId, view] : veta->views) {
            auto SE3_CurCmToW = CurCmToW(view->timestamp, topic);
            if (SE3_CurCmToW == std::nullopt) {
                throw Status(Status::CRITICAL,
                             "can not find pose from B-splines for camera '{}'!!!", topic);
            }
            veta->poses.at(view->poseId) =
                ns_veta::Posed(SE3_CurCmToW->so3(), SE3_CurCmToW->translation());
        }

        for (const auto &reprojCorrSeq : reprojCorrVec) {
            // recover point in camera frame
            Eigen::Vector2d pInCamPlane = intri->ImgToCam(reprojCorrSeq->firObv.x);
            double depth = *visualGlobalScale * 1.0 / *reprojCorrSeq->invDepthFir;
            Eigen::Vector3d pInCam(pInCamPlane(0) * depth, pInCamPlane(1) * depth, depth);
            // transform point to world frame (we do not consider the RS effect here, which only
            // affects the visualization)
            auto pose = veta->poses.at(veta->views.at(reprojCorrSeq->firObvViewId)->poseId);
            Eigen::Vector3d pInW = pose.Rotation() * pInCam + pose.Translation();
            veta->structure.at(reprojCorrSeq->lmId).X = pInW;
        }
    }

    // these quantities need to be backup for Hessian matrix finding in ceres
    auto backUp = std::make_shared<BackUp>();
    backUp->estimator = estimator;
    backUp->visualGlobalScale = visualGlobalScale;
    // the inverse depth of each corr sequence is stored here
    backUp->visualCorrs = visualCorrs;
    return backUp;
}
}  // namespace ns_ikalibr