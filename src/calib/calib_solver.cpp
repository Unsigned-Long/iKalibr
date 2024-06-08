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
#include "calib/estimator.h"
#include "core/rotation_estimator.h"
#include "core/lidar_odometer.h"
#include "core/scan_undistortion.h"
#include "core/pts_association.h"
#include "pangolin/display/display.h"
#include "core/rot_only_vo.h"
#include "veta/landmark.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "core/colmap_data_io.h"
#include "core/visual_reproj_association.h"
#include "opencv2/highgui.hpp"
#include "util/tqdm.h"

_3_

namespace ns_ikalibr {

    // ----------
    // ImagesInfo
    // ----------

    std::optional<std::string> ImagesInfo::GetImagePath(ns_veta::IndexT id) const {
        auto iter = images.find(id);
        if (iter == images.cend()) { return {}; }
        else { return root_path + '/' + iter->second; }
    }

    std::optional<std::string> ImagesInfo::GetImageFilename(ns_veta::IndexT id) const {
        auto iter = images.find(id);
        if (iter == images.cend()) { return {}; }
        else { return iter->second; }
    }

    std::map<ns_veta::IndexT, std::string> ImagesInfo::GetImagesIdxToName() const {
        return images;
    }

    std::map<std::string, ns_veta::IndexT> ImagesInfo::GetImagesNameToIdx() const {
        std::map<std::string, ns_veta::IndexT> imgsInvKV;
        for (const auto &[k, v]: images) { imgsInvKV.insert({v, k}); }
        return imgsInvKV;
    }

    // -----------
    // CalibSolver
    // -----------

    CalibSolver::CalibSolver(CalibDataManager::Ptr calibDataManager, CalibParamManager::Ptr calibParamManager)
            : _dataMagr(std::move(calibDataManager)), _parMagr(std::move(calibParamManager)),
              _ceresOption(Estimator::DefaultSolverOptions(
                      Configor::Preference::AvailableThreads(), true, Configor::Preference::UseCudaInSolving)
              ), _viewer(nullptr), _solveFinished(false) {
        // create so3 and linear scale splines given start and end times, knot distances
        _splines = CreateSplineBundle(
                _dataMagr->GetCalibStartTimestamp(), _dataMagr->GetCalibEndTimestamp(),
                Configor::Prior::KnotTimeDist::SO3Spline, Configor::Prior::KnotTimeDist::ScaleSpline
        );

        // create viewer
        _viewer = Viewer::Create(_parMagr, _splines);
        if (!Configor::IsCameraIntegrated() && !Configor::IsLiDARIntegrated()) {
            auto modelPath = ros::package::getPath("ikalibr") + "/model/ikalibr.obj";
            if (std::filesystem::exists(modelPath)) {
                _viewer->AddObjEntity(modelPath, Viewer::VIEW_MAP);
                _viewer->AddObjEntity(modelPath, Viewer::VIEW_ASSOCIATION);
            } else {
                spdlog::warn("can not load models from '{}'!", modelPath);
            }
        }

        // pass the 'CeresViewerCallBack' to ceres option so that update the viewer after every interation in ceres
        _ceresOption.callbacks.push_back(new CeresViewerCallBack(_viewer));
        _ceresOption.update_state_every_iteration = true;

        // output spatiotemporal parameters after each iteration if needed
        if (Configor::Preference::OutputParamInEachIter) {
            _ceresOption.callbacks.push_back(new CeresDebugCallBack(_parMagr));
        }
    }

    CalibSolver::Ptr CalibSolver::Create(const CalibDataManager::Ptr &calibDataManager,
                                         const CalibParamManager::Ptr &calibParamManager) {
        return std::make_shared<CalibSolver>(calibDataManager, calibParamManager);
    }

    void CalibSolver::Process() {
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
                        DataAssociationForCameras()
                );
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
                        DataAssociationForCameras()
                );
                // 'curGlobalMap' and 'curUndistFramesInMap' would be deconstructed here
            }

            // align states to the gravity
            AlignStatesToGravity();
            _viewer->UpdateSplineViewer();

            _parMagr->ShowParamStatus();

            if (Configor::Preference::OutputParamInEachIter) {
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
            _backup->lidarCorrs = DataAssociationForLiDARs(std::get<0>(final), std::get<1>(final), 100000);
        }
        if (Configor::IsRadarIntegrated() && GetScaleType() == TimeDeriv::LIN_POS_SPLINE) {
            spdlog::info("build final radar map...");
            // radar map would be added to the viewer in this function
            _backup->radarMap = BuildGlobalMapOfRadar();
        }
        if (Configor::IsCameraIntegrated()) {
            for (const auto &[topic, sfmData]: _dataMagr->GetSfMData()) {
                _viewer->AddVeta(sfmData, Viewer::VIEW_MAP);
            }
        }

        _solveFinished = true;

        spdlog::info("Solving is finished! Focus on the viewer and press [ctrl+'s'] to save the current scene!");
        spdlog::info("Focus on the viewer and press ['w', 's', 'a', 'd'] to zoom spline viewer!");
    }

    std::tuple<IKalibrPointCloud::Ptr, std::map<std::string, std::vector<LiDARFrame::Ptr>>>
    CalibSolver::Initialization() {
        const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
        // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are poor
        const double st = std::max(so3Spline.MinTime(), scaleSpline.MinTime()) + Configor::Prior::TimeOffsetPadding;
        const double et = std::min(so3Spline.MaxTime(), scaleSpline.MaxTime()) - Configor::Prior::TimeOffsetPadding;

        // ----------------------------------------------------------------------------------------------------
        // Step 1: initialize (i) so3 spline of the reference IMU, (ii) extrinsic rotations, (iii) time offsets
        // ----------------------------------------------------------------------------------------------------
        spdlog::info("fitting rotation b-spline...");

        // here we  recover the so3 spline first using only the angular velocities from the reference IMU,
        // then estimates other quantities by involving angular velocity measurements from other IMUs.
        // For better readability, we could also optimize them together (not current version)
        auto estimator = Estimator::Create(_splines, _parMagr);
        this->AddGyroFactor(estimator, Configor::DataStream::ReferIMU, OptOption::Option::OPT_SO3_SPLINE);
        auto sum = estimator->Solve(_ceresOption);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

        estimator = Estimator::Create(_splines, _parMagr);
        auto optOption = OptOption::Option::OPT_SO3_BiToBr;
        if (Configor::Prior::OptTemporalParams) { optOption |= OptOption::Option::OPT_TO_BiToBr; }
        for (const auto &[topic, _]: Configor::DataStream::IMUTopics) {
            this->AddGyroFactor(estimator, topic, optOption);
        }
        // make this problem full rank
        estimator->SetRefIMUParamsConstant();
        // estimator->FixFirSO3ControlPoint();

        sum = estimator->Solve(_ceresOption);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

        if (Configor::Preference::OutputParamInEachIter) { SaveStageCalibParam(_parMagr, "stage_1_rot_fit"); }

        // --------------------------------------------------------------------
        // Step 2.1: perform rotation-only visual odometer to recover rotations
        // --------------------------------------------------------------------
        std::shared_ptr<tqdm> bar = nullptr;
        std::map<std::string, RotOnlyVisualOdometer::Ptr> rotOnlyOdom;
        if (Configor::IsCameraIntegrated()) {
            // how many features to maintain in each image
            constexpr int featNumPerImg = 300;
            // the min distance between two features (to ensure features are distributed uniformly)
            constexpr int minDist = 25;
            for (const auto &[topic, frameVec]: _dataMagr->GetCameraMeasurements()) {
                spdlog::info("perform rotation-only visual odometer to recover extrinsic rotations for '{}'...", topic);

                // estimates rotations
                auto odometer = RotOnlyVisualOdometer::Create(featNumPerImg, minDist, _parMagr->INTRI.Camera.at(topic));
                // estimates extrinsic rotation between the camera and reference IMU using the estimated rotations
                auto rotEstimator = RotationEstimator::Create();

                bar = std::make_shared<tqdm>();
                for (int i = 0; i < static_cast<int>(frameVec.size()); ++i) {
                    bar->progress(i, static_cast<int>(frameVec.size()));

                    // if tracking current frame failed, the rotation-only odometer would re-initialize
                    if (!odometer->GrabFrame(frameVec.at(i))) {
                        spdlog::warn("tracking failed when grab the '{}' image frame!!! try to reinitialize", i);
                    }

                    // we do not want to try to recover the extrinsic rotation too frequent
                    if ((odometer->GetRotations().size() < 50) || (odometer->GetRotations().size() % 5 != 0)) {
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
                    throw Status(
                            Status::ERROR,
                            "initialize rotation 'SO3_CmToBr' failed, this may be related to insufficiently excited motion or bad images."
                    );
                } else {
                    spdlog::info(
                            "extrinsic rotation of '{}' is recovered using '{:06}' frames", topic,
                            odometer->GetRotations().size()
                    );
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
            // where only the extrinsic rotations are considered, but here we use a continuous-time rotation-only hand-eye alignment
            // where both extrinsic rotations and temporal parameters are considered.
            spdlog::info("perform rotation alignment to initialize time offset for each camera...");
            estimator = Estimator::Create(_splines, _parMagr);
            optOption = OptOption::Option::OPT_SO3_CmToBr;
            if (Configor::Prior::OptTemporalParams) { optOption |= OptOption::Option::OPT_TO_CmToBr; }

            // perform time offset estimation and extrinsic rotation refinement
            for (const auto &[topic, _]: Configor::DataStream::CameraTopics) {
                const auto &rotations = rotOnlyOdom.at(topic)->GetRotations();
                // this field should be zero
                double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(topic);
                double weight = Configor::DataStream::CameraTopics.at(topic).Weight;

                for (int i = 0; i < static_cast<int>(rotations.size()) - 1; ++i) {
                    const auto &sPose = rotations.at(i), ePose = rotations.at(i + 1);
                    // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are poor
                    if (sPose.first + TO_CmToBr < st || ePose.first + TO_CmToBr > et) { continue; }

                    estimator->AddHandEyeRotationAlignmentForCamera(
                            topic, sPose.first, ePose.first, sPose.second, ePose.second, optOption, weight
                    );
                }
            }

            sum = estimator->Solve(_ceresOption);
            spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
        }

        // ---------------------------------
        // Step 2.3: perform SfM for cameras
        // ---------------------------------
        if (ns_ikalibr::Configor::IsCameraIntegrated()) {
            spdlog::info("perform SfM for each camera...");
            int needSfMCount = 0;
            for (const auto &[topic, data]: _dataMagr->GetCameraMeasurements()) {
                // load data if SfM has been performed
                auto veta = TryLoadSfMData(
                        topic,
                        // for rs camera, as the rs effect is not considered in SfM, we relax the landmark selection condition
                        IsRSCamera(topic) ? 2.0 : 1.0,
                        Configor::DataStream::CameraTopics.at(topic).TrackLengthMin
                );
                if (veta != nullptr) {
                    spdlog::info("SfM data for camera '{}' is valid!", topic);
                    spdlog::info("down sample SfM data for camera '{}'", topic);
                    // keep too many landmarks and features in estimator is not always good
                    DownsampleVeta(veta, 10000, 10);
                    // we store SfM datas in '_dataMagr' as it is the 'calibration data manager'
                    _dataMagr->SetSfMData(topic, veta);
                    // just for visualization
                    _viewer->AddVeta(veta, Viewer::VIEW_MAP);
                    spdlog::info("SfM info for topic '{}' after filtering: view count: {}, landmark count: {}",
                                 topic, veta->views.size(), veta->structure.size());
                    continue;
                } else {
                    spdlog::warn("SfM data for camera '{}' is not valid! Perform SfM using colmap first!", topic);
                }

                // if the camera is integrated, output images for SfM using thirdparty
                // note that the rotation priors of each frame from the extrinsic rotation and so3 spline
                // are utilized in this process to accelerate the feature matching
                auto sfm = VisionOnlySfM::Create(topic, data, _parMagr, so3Spline, _viewer);
                // we use the thirdparty library, i.e., colmap, to solve the SfM problem, rather than ours
                // sfm->PreProcess();
                // sfm->StructureFromMotion();
                spdlog::info("store images of '{}' for SfM...", topic);
                // output image frames and their corresponding information, which would be re-load after SfM are performed
                StoreImagesForSfM(topic, sfm->FindCovisibility(0.1));
                ++needSfMCount;
            }
            if (needSfMCount != 0) {
                throw ns_ikalibr::Status(
                        Status::FINE,
                        "images have been output to '{}/images'. "
                        "Perform SfM for each camera using command lines in file '/sfm_ws/sfm-command-line.txt'!",
                        ns_ikalibr::Configor::DataStream::OutputPath
                );
            }
        }

        // -----------------------------------------------------------------------
        // Step 2.4: refine time offsets of Cameras by hand eye rotation alignment
        // -----------------------------------------------------------------------
        if (Configor::IsCameraIntegrated()) {
            // perhaps this is not too necessary. The difference from the last continuous-time
            // rotation-only hand-eye alignment is that we use the rotations recovered by colmap-derived ones,
            // rather than rotation-only odometer-derived ones.
            spdlog::info("perform rotation alignment to refine time offset for each camera using SfM results...");
            estimator = Estimator::Create(_splines, _parMagr);
            optOption = OptOption::Option::OPT_SO3_CmToBr;
            if (Configor::Prior::OptTemporalParams) { optOption |= OptOption::Option::OPT_TO_CmToBr; }

            for (const auto &[camTopic, veta]: _dataMagr->GetSfMData()) {
                double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;

                const auto &frames = _dataMagr->GetCameraMeasurements(camTopic);
                std::vector<ns_ctraj::Posed> constructedFrames;
                constructedFrames.reserve(frames.size());

                // find constructed frames by SfM
                for (const auto &frame: frames) {
                    auto viewIter = veta->views.find(frame->GetId());
                    if (viewIter == veta->views.cend()) { continue; }
                    auto poseIter = veta->poses.find(viewIter->second->poseId);
                    if (poseIter == veta->poses.cend()) { continue; }
                    const auto &pose = poseIter->second;
                    constructedFrames.emplace_back(pose.Rotation(), pose.Translation(), frame->GetTimestamp());
                }

                for (int i = 0; i < static_cast<int>(constructedFrames.size()) - 1; ++i) {
                    const auto &sPose = constructedFrames.at(i), ePose = constructedFrames.at(i + 1);
                    // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are poor
                    if (sPose.timeStamp + TO_CmToBr < st || ePose.timeStamp + TO_CmToBr > et) { continue; }

                    estimator->AddHandEyeRotationAlignmentForCamera(
                            camTopic, sPose.timeStamp, ePose.timeStamp, sPose.so3, ePose.so3, optOption, weight
                    );
                }
            }

            sum = estimator->Solve(_ceresOption);
            spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
        }

        // -----------------------------------------------------------------------
        // Step 3.1: initialize extrinsic rotations of LiDARs if they are involved
        // -----------------------------------------------------------------------
        if (Configor::IsLiDARIntegrated()) {
            spdlog::info("LiDARs are integrated, initializing extrinsic rotations of LiDARs...");
            for (const auto &[topic, data]: _dataMagr->GetLiDARMeasurements()) {
                spdlog::info("performing ndt odometer for '{}' for extrinsic rotation initialization...", topic);

                auto lidarOdometer = LiDAROdometer::Create(
                        static_cast<float>(Configor::Prior::NDTLiDAROdometer::Resolution),
                        Configor::Preference::AvailableThreads()
                );
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
                    throw Status(
                            Status::ERROR,
                            "initialize rotation 'SO3_LkToBr' failed, this may be related to the 'NDTResolution' of lidar odometer."
                    );
                } else {
                    spdlog::info(
                            "extrinsic rotation of '{}' is recovered using '{:06}' frames", topic,
                            lidarOdometer->GetOdomPoseVec().size()
                    );
                }
                _viewer->AddCloud(
                        lidarOdometer->GetMap(), Viewer::VIEW_MAP, ns_viewer::Entity::GetUniqueColour(), 2.0f
                );
                _viewer->UpdateSensorViewer();
            }
            _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
        }

        // --------------------------------------------------
        // Step 3.2: Undistort lidar scans and rerun odometer
        // --------------------------------------------------
        std::map<std::string, LiDAROdometer::Ptr> lidarOdometers;
        std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFramesInScan;
        if (Configor::IsLiDARIntegrated()) {

            auto undistHelper = ScanUndistortion::Create(_splines, _parMagr);

            for (const auto &[topic, data]: _dataMagr->GetLiDARMeasurements()) {
                spdlog::info("undistort scans for lidar '{}'...", topic);

                // undistort rotation only using 'UNDIST_SO3' in initialization
                undistFramesInScan[topic] = undistHelper->UndistortToScan(
                        data, topic, ScanUndistortion::Option::UNDIST_SO3
                );
                const auto &undistFrames = undistFramesInScan.at(topic);

                spdlog::info("rerun odometer for lidar '{}' using undistorted scans...", topic);
                lidarOdometers[topic] = LiDAROdometer::Create(
                        static_cast<float>(Configor::Prior::NDTLiDAROdometer::Resolution),
                        Configor::Preference::AvailableThreads()
                );
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

                        if (curUndistFrame == nullptr || lastUndistFrame == nullptr) { continue; }

                        auto curLtoRef = this->CurLkToW(curUndistFrame->GetTimestamp(), topic);
                        auto lastLtoRef = this->CurLkToW(lastUndistFrame->GetTimestamp(), topic);

                        // if query pose successfully
                        if (curLtoRef && lastLtoRef) {
                            // the translation has not been initialized
                            predCurToLast = ns_ctraj::Posed(
                                    (*lastLtoRef).so3().inverse() * (*curLtoRef).so3(), Eigen::Vector3d::Zero()
                            ).T();
                        } else {
                            predCurToLast = Eigen::Matrix4d::Identity();
                        }
                    }
                    lidarOdometers.at(topic)->FeedFrame(curUndistFrame, predCurToLast, i < 100);
                }
                bar->finish();

                _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
                _viewer->AddCloud(
                        lidarOdometers.at(topic)->GetMap(), Viewer::VIEW_MAP, ns_viewer::Entity::GetUniqueColour(), 2.0f
                );
            }
        }

        // --------------------------------------------------------------------------
        // Step 3.3: initialize time offsets of LiDARs by hand eye rotation alignment
        // --------------------------------------------------------------------------
        if (Configor::IsLiDARIntegrated()) {

            spdlog::info("performing  hand eye rotation alignment for LiDARs...");
            estimator = Estimator::Create(_splines, _parMagr);
            optOption = OptOption::Option::OPT_SO3_LkToBr;
            if (Configor::Prior::OptTemporalParams) { optOption |= OptOption::Option::OPT_TO_LkToBr; }

            for (const auto &[lidarTopic, odometer]: lidarOdometers) {
                const auto &poseSeq = odometer->GetOdomPoseVec();
                double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
                double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;

                for (int i = 0; i < static_cast<int>(poseSeq.size()) - 1; ++i) {
                    const auto &sPose = poseSeq.at(i), ePose = poseSeq.at(i + 1);
                    // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are poor
                    if (sPose.timeStamp + TO_LkToBr < st || ePose.timeStamp + TO_LkToBr > et) { continue; }

                    estimator->AddHandEyeRotationAlignmentForLiDAR(
                            lidarTopic, sPose.timeStamp, ePose.timeStamp, sPose.so3, ePose.so3, optOption, weight
                    );
                }
            }

            sum = estimator->Solve(_ceresOption);
            spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
        }

        // ---------------------------------------------------------------------------
        // Step 4: initialize other spatial parameters using inertial-sensor alignment
        // ---------------------------------------------------------------------------
        spdlog::info("performing inertial alignment to initialize other spatial parameters...");
        estimator = Estimator::Create(_splines, _parMagr);

        // we do not optimization the already initialized extrinsic rotations (IMUs', Cameras', and LiDARs') here
        optOption =
                // lidar extrinsic translations
                OptOption::Option::OPT_POS_LkInBr |
                // camera extrinsic translations and visual scale
                OptOption::Option::OPT_POS_CmInBr | OptOption::Option::OPT_VISUAL_GLOBAL_SCALE |
                // radar extrinsics
                OptOption::Option::OPT_SO3_RjToBr | OptOption::Option::OPT_POS_RjInBr |
                // imu extrinsic translations and gravity
                OptOption::Option::OPT_POS_BiInBr | OptOption::Option::OPT_GRAVITY;

        // lidar-inertial alignment
        std::map<std::string, std::vector<Eigen::Vector3d>> linVelSeqLk;
        for (const auto &[lidarTopic, odometer]: lidarOdometers) {

            const auto &poseSeq = odometer->GetOdomPoseVec();
            double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);

            // create linear velocity sequence
            linVelSeqLk[lidarTopic] = std::vector<Eigen::Vector3d>(poseSeq.size(), Eigen::Vector3d::Zero());
            auto &curLidarLinVelSeq = linVelSeqLk.at(lidarTopic);
            double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;

            for (const auto &[imuTopic, imuFrames]: _dataMagr->GetIMUMeasurements()) {
                spdlog::info("add lidar-inertial alignment factors for '{}' and '{}'...", lidarTopic, imuTopic);

                static constexpr int STEP = 5;
                for (int i = 0; i < static_cast<int>(poseSeq.size()) - STEP; ++i) {
                    const auto &sPose = poseSeq.at(i), ePose = poseSeq.at(i + STEP);
                    // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are poor
                    if (sPose.timeStamp + TO_LkToBr < st || ePose.timeStamp + TO_LkToBr > et) { continue; }

                    if (ePose.timeStamp - sPose.timeStamp < 1E-3 || ePose.timeStamp - sPose.timeStamp > 1.0) {
                        continue;
                    }

                    estimator->AddLiDARInertialAlignment(
                            imuFrames, lidarTopic, imuTopic, sPose, ePose, odometer->GetMapTime(),
                            &curLidarLinVelSeq.at(i), &curLidarLinVelSeq.at(i + STEP), optOption, weight
                    );
                }
            }
        }

        // visual-inertial alignment
        std::map<std::string, std::vector<Eigen::Vector3d>> linVelSeqCm;
        std::map<std::string, double> visualScaleSeq;
        std::map<std::string, std::vector<ns_ctraj::Posed>> sfmPoseSeq;
        for (const auto &[camTopic, veta]: _dataMagr->GetSfMData()) {
            double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);

            const auto &frames = _dataMagr->GetCameraMeasurements(camTopic);
            // first frame to world
            ns_veta::Posed FirCtoW;
            double firCTime = 0.0;
            for (const auto &frame: frames) {
                auto viewIter = veta->views.find(frame->GetId());
                if (viewIter == veta->views.cend()) { continue; }
                auto poseIter = veta->poses.find(viewIter->second->poseId);
                if (poseIter == veta->poses.cend()) { continue; }
                FirCtoW = poseIter->second, firCTime = frame->GetTimestamp();
                break;
            }
            // set first valid camera frame as world frame (transform from world to first frame)
            PerformTransformForVeta(veta, FirCtoW.Inverse(), 1.0);

            // load poses
            sfmPoseSeq[camTopic] = {};
            auto &constructedFrames = sfmPoseSeq.at(camTopic);
            constructedFrames.reserve(frames.size());
            for (const auto &frame: frames) {
                auto viewIter = veta->views.find(frame->GetId());
                if (viewIter == veta->views.cend()) { continue; }
                auto poseIter = veta->poses.find(viewIter->second->poseId);
                if (poseIter == veta->poses.cend()) { continue; }
                const auto &pose = poseIter->second;
                constructedFrames.emplace_back(pose.Rotation(), pose.Translation(), frame->GetTimestamp());
            }

            // create linear velocity sequence
            linVelSeqCm[camTopic] = std::vector<Eigen::Vector3d>(constructedFrames.size(), Eigen::Vector3d::Zero());
            auto &curCamLinVelSeq = linVelSeqCm.at(camTopic);
            double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;
            // create scale
            visualScaleSeq[camTopic] = 1.0;
            auto &scale = visualScaleSeq.at(camTopic);

            for (const auto &[imuTopic, imuFrames]: _dataMagr->GetIMUMeasurements()) {
                spdlog::info("add visual-inertial alignment factors for '{}' and '{}'...", camTopic, imuTopic);

                static constexpr int STEP = 5;
                for (int i = 0; i < static_cast<int>(constructedFrames.size()) - STEP; ++i) {
                    const auto &sPose = constructedFrames.at(i), ePose = constructedFrames.at(i + STEP);
                    // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are poor
                    if (sPose.timeStamp + TO_CmToBr < st || ePose.timeStamp + TO_CmToBr > et) { continue; }

                    if (ePose.timeStamp - sPose.timeStamp < 1E-3 || ePose.timeStamp - sPose.timeStamp > 1.0) {
                        continue;
                    }

                    estimator->AddVisualInertialAlignment(
                            imuFrames, camTopic, imuTopic, sPose, ePose, firCTime,
                            &curCamLinVelSeq.at(i), &curCamLinVelSeq.at(i + STEP), &scale, optOption, weight
                    );
                }
            }
        }

        // inertial alignment (only when more than or equal to 2 num IMUs are involved)
        constexpr double dt = 0.1;
        std::vector<Eigen::Vector3d> linVelSeqBr(std::floor((et - st) / dt), Eigen::Vector3d::Zero());
        if (Configor::DataStream::IMUTopics.size() >= 2) {
            for (const auto &[topic, frames]: _dataMagr->GetIMUMeasurements()) {
                spdlog::info("add inertial alignment factors for '{}'...", topic);

                for (int i = 0; i < static_cast<int>(linVelSeqBr.size()) - 1; ++i) {
                    int sIdx = i, eIdx = i + 1;
                    double sTimeByBr = sIdx * dt + st, eTimeByBr = eIdx * dt + st;
                    Eigen::Vector3d *sVel = &linVelSeqBr.at(sIdx), *eVel = &linVelSeqBr.at(eIdx);

                    estimator->AddInertialAlignment(
                            frames, topic, sTimeByBr, eTimeByBr, sVel, eVel, optOption,
                            Configor::DataStream::IMUTopics.at(topic).AcceWeight
                    );
                }
            }
        }

        // radar-inertial alignment
        for (const auto &[radarTopic, radarMes]: _dataMagr->GetRadarMeasurements()) {
            double weight = Configor::DataStream::RadarTopics.at(radarTopic).Weight;
            double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(radarTopic);

            for (const auto &[imuTopic, frames]: _dataMagr->GetIMUMeasurements()) {
                spdlog::info("add radar-inertial alignment factors for '{}' and '{}'...", radarTopic, imuTopic);

                for (int i = 0; i < static_cast<int>(radarMes.size()) - 1; ++i) {
                    const auto &sArray = radarMes.at(i), eArray = radarMes.at(i + 1);

                    //spdlog::info("sAry count: {}, eAry count: {}",
                    //             sArray->GetTargets().size(), eArray->GetTargets().size());

                    // to estimate the radar velocity by linear least-squares solver
                    // the minim targets number required is 3
                    if (sArray->GetTargets().size() < 10 || eArray->GetTargets().size() < 10) { continue; }

                    if (sArray->GetTimestamp() + TO_RjToBr < st || eArray->GetTimestamp() + TO_RjToBr > et) {
                        continue;
                    }

                    estimator->AddRadarInertialAlignment(
                            frames, imuTopic, radarTopic, sArray, eArray, optOption, weight
                    );
                }
            }
        }

        // fix spatiotemporal parameters of reference sensor
        // make this problem full rank
        estimator->SetRefIMUParamsConstant();

        sum = estimator->Solve(_ceresOption);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

        if (Configor::Preference::OutputParamInEachIter) { SaveStageCalibParam(_parMagr, "stage_2_align"); }

        for (const auto &[topic, lidarOdom]: lidarOdometers) {
            _viewer->AddCloud(lidarOdom->GetMap(), Viewer::VIEW_MAP, ns_viewer::Entity::GetUniqueColour(), 2.0f);
        }
        // recover scale fro veta
        for (const auto &[camTopic, veta]: _dataMagr->GetSfMData()) {
            spdlog::info("visual global scale for camera '{}': {:.3f}", camTopic, visualScaleSeq.at(camTopic));
            PerformTransformForVeta(veta, ns_veta::Posed(), visualScaleSeq.at(camTopic));
            _viewer->AddVeta(veta, Viewer::VIEW_MAP);
        }
        // perform scale for 'sfmPoseSeq', which would used for scale spline recovery
        for (auto &[camTopic, poseSeq]: sfmPoseSeq) {
            const double scale = visualScaleSeq.at(camTopic);
            for (auto &pose: poseSeq) { pose.t *= scale; }
        }

        // deconstruct data
        linVelSeqLk.clear(), linVelSeqBr.clear(), linVelSeqCm.clear(), visualScaleSeq.clear();

        // ----------------------------------------------------------------------------
        // Step 5: recover scale spline (linear acceleration, velocity, or translation)
        // ----------------------------------------------------------------------------
        spdlog::info("performing scale spline recovery...");

        estimator = Estimator::Create(_splines, _parMagr);
        optOption = OptOption::Option::OPT_SCALE_SPLINE;

        switch (GetScaleType()) {
            case TimeDeriv::LIN_ACCE_SPLINE: {
                // only multiple imus are involved
                this->AddAcceFactor<TimeDeriv::LIN_ACCE_SPLINE>(estimator, Configor::DataStream::ReferIMU, optOption);
            }
                break;
            case TimeDeriv::LIN_VEL_SPLINE: {
                // only multiple radars and imus are involved
                constexpr int VelDeriv = TimeDeriv::Deriv<TimeDeriv::LIN_VEL_SPLINE, TimeDeriv::LIN_VEL>();
                optOption |= OptOption::Option::OPT_SO3_SPLINE;
                for (const auto &[topic, data]: _dataMagr->GetRadarMeasurements()) {
                    const double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(topic);
                    const auto &SO3_RjToBr = _parMagr->EXTRI.SO3_RjToBr.at(topic);
                    const Eigen::Vector3d &POS_RjInBr = _parMagr->EXTRI.POS_RjInBr.at(topic);
                    const double weight = Configor::DataStream::RadarTopics.at(topic).Weight;

                    for (const auto &ary: data) {
                        double timeByBr = ary->GetTimestamp() + TO_RjToBr;
                        if (ary->GetTargets().size() < 10 || !so3Spline.TimeStampInRange(timeByBr)) { continue; }

                        auto SO3_BrToW = so3Spline.Evaluate(timeByBr);
                        auto ANG_VEL_BrToWInW = SO3_BrToW * so3Spline.VelocityBody(timeByBr);
                        Eigen::Vector3d LIN_VEL_BrToWInW =
                                SO3_BrToW * SO3_RjToBr * ary->RadarVelocityFromStaticTargetArray() -
                                Sophus::SO3d::hat(ANG_VEL_BrToWInW) * (SO3_BrToW * POS_RjInBr);

                        estimator->AddLinearScaleConstraint<VelDeriv>(timeByBr, LIN_VEL_BrToWInW, optOption, weight);
                    }
                }
                // add acceleration factor using inertial measurements only from reference IMU
                this->AddAcceFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, Configor::DataStream::ReferIMU, optOption);
                this->AddGyroFactor(estimator, Configor::DataStream::ReferIMU, optOption);
                // if optimize time offsets, we first recover the scale spline, then construct a ls problem to recover time offsets
                if (Configor::Prior::OptTemporalParams) {
                    // we don't want to output the solving information
                    estimator->Solve(Estimator::DefaultSolverOptions(
                            Configor::Preference::AvailableThreads(), false, Configor::Preference::UseCudaInSolving
                    ));

                    estimator = Estimator::Create(_splines, _parMagr);
                    optOption = OptOption::Option::OPT_TO_RjToBr;
                    for (const auto &[topic, _]: Configor::DataStream::RadarTopics) {
                        this->AddRadarFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
                    }
                }

            }
                break;
            case TimeDeriv::LIN_POS_SPLINE: {
                CalibSolver::SplineBundleType::Ptr roughSplines;
                // this value equals to zero
                constexpr int PosDeriv = TimeDeriv::Deriv<TimeDeriv::LIN_POS_SPLINE, TimeDeriv::LIN_POS>();
                optOption |= OptOption::Option::OPT_SO3_SPLINE;

                std::vector<double> headTimeVec, tailTimeVec;
                for (const auto &[lidarTopic, odometer]: lidarOdometers) {
                    const double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
                    headTimeVec.push_back(odometer->GetOdomPoseVec().front().timeStamp + TO_LkToBr);
                    tailTimeVec.push_back(odometer->GetOdomPoseVec().back().timeStamp + TO_LkToBr);
                }
                for (const auto &[camTopic, poseSeq]: sfmPoseSeq) {
                    const double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                    headTimeVec.push_back(poseSeq.front().timeStamp + TO_CmToBr);
                    tailTimeVec.push_back(poseSeq.back().timeStamp + TO_CmToBr);
                }

                const double minTime = *std::min_element(headTimeVec.begin(), headTimeVec.end()) - 1E-3;
                const double maxTime = *std::max_element(tailTimeVec.begin(), tailTimeVec.end()) + 1E-3;

                // create rough splines (the time distance is larger than that from configure as our poses are not dense)
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

                for (const auto &[lidarTopic, odometer]: lidarOdometers) {
                    double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);
                    double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;
                    auto SE3_BrToLk = _parMagr->EXTRI.SE3_LkToBr(lidarTopic).inverse();

                    // the translation part is not exactly rigorous as the translation spline is not recovered yet
                    // this assumes that the origin of map frame of {Lk} is the same as that of the world frame {br0}
                    auto SE3_Lk0ToBr0 = this->CurLkToW(odometer->GetMapTime(), lidarTopic);

                    if (SE3_Lk0ToBr0 == std::nullopt) {
                        throw Status(
                                Status::CRITICAL,
                                "map time of '{}' is out of time range of splines!", lidarTopic
                        );
                    }

                    for (const auto &SE3_LkToLk0: odometer->GetOdomPoseVec()) {
                        auto SE3_BrToBr0 = *SE3_Lk0ToBr0 * SE3_LkToLk0.se3() * SE3_BrToLk;
                        estimator->AddSO3Constraint(
                                SE3_LkToLk0.timeStamp + TO_LkToBr, SE3_BrToBr0.so3(), optOption, weight
                        );
                        estimator->AddLinearScaleConstraint<PosDeriv>(
                                SE3_LkToLk0.timeStamp + TO_LkToBr, SE3_BrToBr0.translation(), optOption, weight
                        );
                    }
                }

                for (const auto &[camTopic, poseSeq]: sfmPoseSeq) {
                    double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
                    double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;
                    auto SE3_BrToCm = _parMagr->EXTRI.SE3_CmToBr(camTopic).inverse();

                    // the translation part is not exactly rigorous as the translation spline is not recovered yet
                    // this assumes that the origin of map frame of {Cm} is the same as that of the world frame {br0}
                    auto SE3_Cm0ToBr0 = this->CurCmToW(poseSeq.front().timeStamp, camTopic);

                    if (SE3_Cm0ToBr0 == std::nullopt) {
                        throw Status(
                                Status::CRITICAL,
                                "map time of '{}' is out of time range of splines!", camTopic
                        );
                    }

                    for (const auto &SE3_CmToCm0: poseSeq) {
                        auto SE3_BrToBr0 = *SE3_Cm0ToBr0 * SE3_CmToCm0.se3() * SE3_BrToCm;
                        estimator->AddSO3Constraint(
                                SE3_CmToCm0.timeStamp + TO_CmToBr, SE3_BrToBr0.so3(), optOption, weight
                        );
                        estimator->AddLinearScaleConstraint<PosDeriv>(
                                SE3_CmToCm0.timeStamp + TO_CmToBr, SE3_BrToBr0.translation(), optOption, weight
                        );
                    }
                }

                // add tail factors (constraints) to maintain enough observability
                estimator->AddLinScaleTailConstraint(optOption, 1.0);
                estimator->AddSO3TailConstraint(optOption, 1.0);

                // we don't want to output the solving information
                estimator->Solve(Estimator::DefaultSolverOptions(
                        Configor::Preference::AvailableThreads(), false, Configor::Preference::UseCudaInSolving
                ));

                const auto &rSo3Spline = roughSplines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
                const auto &rScaleSpline = roughSplines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

                estimator = Estimator::Create(_splines, _parMagr);

                for (double t = minTime; t < maxTime;) {
                    estimator->AddSO3Constraint(t, rSo3Spline.Evaluate(t), optOption, 1.0);
                    estimator->AddLinearScaleConstraint<PosDeriv>(t, rScaleSpline.Evaluate(t), optOption, 1.0);
                    t += 0.01;
                }
            }
                break;
        }
        sum = estimator->Solve(_ceresOption);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

        if (GetScaleType() == TimeDeriv::LIN_POS_SPLINE &&
            Configor::IsRadarIntegrated() && Configor::Prior::OptTemporalParams) {
            // in this case, the time offsets of radars have not been recovered
            estimator = Estimator::Create(_splines, _parMagr);
            for (const auto &[topic, _]: Configor::DataStream::RadarTopics) {
                this->AddRadarFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic, OptOption::Option::OPT_TO_RjToBr);
            }
            sum = estimator->Solve(_ceresOption);
            spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
        }

        if (Configor::Preference::OutputParamInEachIter) { SaveStageCalibParam(_parMagr, "stage_3_scale_fit"); }

        //------------------------------------------------------
        // Step 6: align initialized states to gravity direction
        //------------------------------------------------------
        spdlog::info("aligning all states to gravity direction...");
        AlignStatesToGravity();
        _viewer->UpdateSplineViewer();

        // --------------------------------------------------------------
        // Step 7.1 : trans veta to world frame if Cameras are integrated
        // --------------------------------------------------------------
        for (const auto &[camTopic, poseSeq]: sfmPoseSeq) {
            auto SE3_Cm0ToBr0 = this->CurCmToW(poseSeq.front().timeStamp, camTopic);

            if (SE3_Cm0ToBr0 == std::nullopt) {
                throw Status(
                        Status::CRITICAL,
                        "map time of '{}' is out of time range of splines!", camTopic
                );
            }
            PerformTransformForVeta(
                    _dataMagr->GetSfMData(camTopic),
                    ns_veta::Posed(SE3_Cm0ToBr0->so3(), SE3_Cm0ToBr0->translation()), 1.0
            );
        }

        // ----------------------------------------------------------------
        // Step 7.2 : build map and undisto frames if LiDARs are integrated
        // ----------------------------------------------------------------
        spdlog::info("build global map and undisto lidar frames in world...");
        IKalibrPointCloud::Ptr globalMap(new IKalibrPointCloud);
        std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFramesInMap;
        for (const auto &[lidarTopic, odometer]: lidarOdometers) {
            auto SE3_Lk0ToBr0 = this->CurLkToW(odometer->GetMapTime(), lidarTopic);

            if (SE3_Lk0ToBr0 == std::nullopt) {
                throw Status(
                        Status::CRITICAL,
                        "map time of '{}' is out of time range of splines!", lidarTopic
                );
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
                    pcl::transformPointCloud(*undistoScan->GetScan(), *scanInBr0, SE3_ScanToBr0.matrix().cast<float>());
                    curUndistFramesInMap.at(i) = LiDARFrame::Create(undistoScan->GetTimestamp(), scanInBr0);

                    *globalMap += *scanInBr0;
                }
            }
        }
        return {globalMap, undistFramesInMap};
    }

    std::tuple<IKalibrPointCloud::Ptr, std::map<std::string, std::vector<LiDARFrame::Ptr>>>
    CalibSolver::BuildGlobalMapOfLiDAR() {
        if (!Configor::IsLiDARIntegrated()) { return {}; }

        // -------------------------------------------
        // undisto frames and build map in world frame
        // -------------------------------------------
        auto undistHelper = ScanUndistortion::Create(_splines, _parMagr);

        std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFrames;
        IKalibrPointCloud::Ptr mapCloud(new IKalibrPointCloud);
        for (const auto &[topic, data]: _dataMagr->GetLiDARMeasurements()) {
            spdlog::info("undistort scans for lidar '{}'...", topic);
            undistFrames[topic] = undistHelper->UndistortToRef(data, topic, ScanUndistortion::Option::ALL);

            spdlog::info("marge scans to map...");
            for (auto &frame: undistFrames.at(topic)) {
                if (frame == nullptr) { continue; }
                *mapCloud += *frame->GetScan();
            }
        }

        IKalibrPointCloud::Ptr newMapCloud(new IKalibrPointCloud);
        std::vector<int> index;
        pcl::removeNaNFromPointCloud(*mapCloud, *newMapCloud, index);

        return {newMapCloud, undistFrames};
    }

    IKalibrPointCloud::Ptr CalibSolver::BuildGlobalMapOfRadar() {
        if (!Configor::IsRadarIntegrated() || GetScaleType() != TimeDeriv::LIN_POS_SPLINE) { return {}; }

        IKalibrPointCloud::Ptr radarCloud(new IKalibrPointCloud);

        for (const auto &[topic, data]: _dataMagr->GetRadarMeasurements()) {
            const double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(topic);
            IKalibrPointCloud::Ptr curRadarCloud(new IKalibrPointCloud);
            curRadarCloud->reserve(data.size() * data.front()->GetTargets().size());
            for (const auto &ary: data) {
                for (const auto &frame: ary->GetTargets()) {
                    auto SE3_CurRjToW = CurRjToW(frame->GetTimestamp(), topic);
                    if (SE3_CurRjToW == std::nullopt) { continue; }
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
        _viewer->AddStarMarkCloud(radarCloud, Viewer::VIEW_MAP);

        return radarCloud;
    }

    std::map<std::string, std::vector<PointToSurfelCorr::Ptr>>
    CalibSolver::DataAssociationForLiDARs(const IKalibrPointCloud::Ptr &map,
                                          const std::map<std::string, std::vector<LiDARFrame::Ptr>> &undistFrames,
                                          int ptsCountInEachScan) {
        if (!Configor::IsLiDARIntegrated()) { return {}; }

        // ---------------------------------
        // Step 1: down sample the map cloud
        // ---------------------------------
        pcl::VoxelGrid<IKalibrPoint> filter;
        filter.setInputCloud(map);
        auto size = static_cast<float>(Configor::Prior::LiDARDataAssociate::MapDownSample);
        filter.setLeafSize(size, size, size);

        IKalibrPointCloud::Ptr mapDownSampled(new IKalibrPointCloud);
        filter.filter(*mapDownSampled);

        _viewer->AddAlignedCloud(mapDownSampled, Viewer::VIEW_MAP, -_parMagr->GRAVITY.cast<float>(), 2.0f);

        // ------------------------------------------------
        // Step 2: perform data association for each frames
        // ------------------------------------------------
        spdlog::info("perform point to surfel association...");
        auto associator = PointToSurfelAssociator::Create(
                // we use the dense map to create data associator for high-perform point-to-surfel search
                map, Configor::Prior::LiDARDataAssociate::MapResolution,
                Configor::Prior::LiDARDataAssociate::MapDepthLevels
        );
        auto condition = PointToSurfelCondition();
        _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
        _viewer->AddSurfelMap(associator->GetSurfelMap(), condition, Viewer::VIEW_ASSOCIATION);
        _viewer->AddCloud(mapDownSampled, Viewer::VIEW_ASSOCIATION, ns_viewer::Colour::Black().WithAlpha(0.2f), 2.0f);

        // deconstruction
        mapDownSampled.reset();

        std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> pointToSurfel;

        std::size_t count = 0;
        std::shared_ptr<tqdm> bar = nullptr;
        for (const auto &[topic, framesInMap]: undistFrames) {
            const auto &rawFrames = _dataMagr->GetLiDARMeasurements(topic);

            // for each scan, we keep 'ptsCountInEachScan' point to surfel corrs
            pointToSurfel[topic] = {};
            auto &curPointToSurfel = pointToSurfel.at(topic);
            bar = std::make_shared<tqdm>();
            for (int i = 0; i < static_cast<int>(framesInMap.size()); ++i) {
                bar->progress(i, static_cast<int>(framesInMap.size()));

                if (framesInMap.at(i) == nullptr || rawFrames.at(i) == nullptr) { continue; }

                auto ptsVec = associator->Association(
                        framesInMap.at(i)->GetScan(), rawFrames.at(i)->GetScan(), condition
                );

                curPointToSurfel.insert(curPointToSurfel.end(), ptsVec.cbegin(), ptsVec.cend());
            }
            bar->finish();

            // downsample
            int expectCount = ptsCountInEachScan * static_cast<int>(rawFrames.size());
            if (static_cast<int>(curPointToSurfel.size()) > expectCount) {
                std::map<ufo::map::Node, std::vector<PointToSurfelCorr::Ptr>> nodes;
                for (const auto &corr: curPointToSurfel) { nodes[corr->node].push_back(corr); }
                std::size_t numEachNode = expectCount / nodes.size() + 1;
                curPointToSurfel.clear();

                // uniform sampling
                std::default_random_engine engine(std::chrono::steady_clock::now().time_since_epoch().count());
                for (const auto &[node, corrs]: nodes) {
                    auto newCorrs = SamplingWoutReplace2(engine, corrs, std::min(corrs.size(), numEachNode));
                    curPointToSurfel.insert(curPointToSurfel.end(), newCorrs.cbegin(), newCorrs.cend());
                }
            }
            count += curPointToSurfel.size();
        }
        spdlog::info("total point to surfel count for LiDARs: {}", count);
        _viewer->AddPointToSurfel(associator->GetSurfelMap(), pointToSurfel, Viewer::VIEW_ASSOCIATION);

        return pointToSurfel;
    }

    std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>> CalibSolver::DataAssociationForCameras() {
        if (!Configor::IsCameraIntegrated()) { return {}; }

        std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>> corrs;
        for (const auto &[topic, sfmData]: _dataMagr->GetSfMData()) {
            spdlog::info("performing visual reprojection data association for camera '{}'...", topic);
            corrs[topic] = VisualReProjAssociator::Create()->Association(
                    *sfmData, _parMagr->INTRI.Camera.at(topic)
            );
            _viewer->AddVeta(sfmData, Viewer::VIEW_MAP);
            spdlog::info("visual reprojection sequences for '{}': {}", topic, corrs.at(topic).size());
        }
        return corrs;
    }

    CalibSolver::BackUp::Ptr CalibSolver::BatchOptimization(OptOption::Option optOption,
                                                            const std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> &ptsCorrs,
                                                            const std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>> &visualCorrs) {

        auto GetOptString = [](OptOption::Option opt) -> std::string {
            std::stringstream stringStream;
            stringStream << opt;
            return stringStream.str();
        };
        spdlog::info("Optimization option: {}", GetOptString(optOption));

        auto estimator = Estimator::Create(_splines, _parMagr);
        auto visualGlobalScale = std::make_shared<double>(1.0);

        switch (GetScaleType()) {
            case TimeDeriv::LIN_ACCE_SPLINE: {
                for (const auto &[topic, _]: Configor::DataStream::IMUTopics) {
                    this->AddAcceFactor<TimeDeriv::LIN_ACCE_SPLINE>(estimator, topic, optOption);
                    this->AddGyroFactor(estimator, topic, optOption);
                }
            }
                break;
            case TimeDeriv::LIN_VEL_SPLINE: {
                for (const auto &[topic, _]: Configor::DataStream::RadarTopics) {
                    this->AddRadarFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
                }
                for (const auto &[topic, _]: Configor::DataStream::IMUTopics) {
                    this->AddAcceFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
                    this->AddGyroFactor(estimator, topic, optOption);
                }
            }
                break;
            case TimeDeriv::LIN_POS_SPLINE: {
                for (const auto &[lidarTopic, corrSeqVec]: ptsCorrs) {
                    this->AddPointToSurfelFactor<TimeDeriv::LIN_POS_SPLINE>(
                            estimator, lidarTopic, corrSeqVec, optOption
                    );
                }
                for (const auto &[camTopic, corrSeqVec]: visualCorrs) {
                    OptOption::Option visualOpt = optOption;
                    if (IsOptionWith(OptOption::Option::OPT_RS_CAM_READOUT_TIME, visualOpt)) {
                        if (IsRSCamera(camTopic)) {
                            spdlog::info("camera '{}' is a rolling shutter (RS) camera, "
                                         "use optimization option 'OPT_RS_CAM_READOUT_TIME'", camTopic);
                        } else {
                            visualOpt ^= OptOption::Option::OPT_RS_CAM_READOUT_TIME;
                            spdlog::info("camera '{}' is a global shutter (GS) camera, "
                                         "remove optimization option 'OPT_RS_CAM_READOUT_TIME'", camTopic);
                        }
                    }
                    this->AddVisualReprojectionFactor<TimeDeriv::LIN_POS_SPLINE>(
                            estimator, camTopic, corrSeqVec, visualGlobalScale.get(), visualOpt
                    );
                }
                for (const auto &[topic, _]: Configor::DataStream::RadarTopics) {
                    this->AddRadarFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic, optOption);
                }
                for (const auto &[topic, _]: Configor::DataStream::IMUTopics) {
                    this->AddAcceFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic, optOption);
                    this->AddGyroFactor(estimator, topic, optOption);
                }
            }
                break;
        }

        // make this problem full rank
        estimator->SetRefIMUParamsConstant();
        // estimator->FixFirSO3ControlPoint();

        auto sum = estimator->Solve(_ceresOption);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

        // for better map consistency in visualization, we update the veta every time
        for (const auto &[topic, reprojCorrVec]: visualCorrs) {
            auto &veta = _dataMagr->GetSfMData(topic);
            auto &intri = _parMagr->INTRI.Camera.at(topic);
            // compute pose
            for (const auto &[viewId, view]: veta->views) {
                auto SE3_CurCmToW = CurCmToW(view->timestamp, topic);
                if (SE3_CurCmToW == std::nullopt) {
                    throw Status(Status::CRITICAL, "can not find pose from B-splines for camera '{}'!!!", topic);
                }
                veta->poses.at(view->poseId) = ns_veta::Posed(SE3_CurCmToW->so3(), SE3_CurCmToW->translation());
            }

            for (const auto &reprojCorrSeq: reprojCorrVec) {
                // recover point in camera frame
                Eigen::Vector2d pInCamPlane = intri->ImgToCam(reprojCorrSeq->firObv.x);
                double depth = *visualGlobalScale * 1.0 / *reprojCorrSeq->invDepthFir;
                Eigen::Vector3d pInCam(pInCamPlane(0) * depth, pInCamPlane(1) * depth, depth);
                // transform point to world frame (we do not consider the RS effect here, which only affects the visualization)
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

    CalibSolver::~CalibSolver() {
        // solving is not performed or not finished as an exception is thrown
        if (!_solveFinished) { pangolin::QuitAll(); }
        // solving is finished (when use 'pangolin::QuitAll()', the window not quit immediately)
        while (_viewer->IsActive()) { std::this_thread::sleep_for(std::chrono::milliseconds(100)); }
    }

    std::optional<Sophus::SE3d> CalibSolver::CurBrToW(double timeByBr) {
        if (GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
            throw Status(Status::CRITICAL, "'CurBrToW' error, scale spline is not translation spline!!!");
        }
        const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        const auto &posSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
        if (!so3Spline.TimeStampInRange(timeByBr) || !posSpline.TimeStampInRange(timeByBr)) {
            return {};
        } else {
            return Sophus::SE3d(so3Spline.Evaluate(timeByBr), posSpline.Evaluate(timeByBr));
        }
    }

    std::optional<Sophus::SE3d> CalibSolver::CurLkToW(double timeByLk, const std::string &topic) {
        if (GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
            throw Status(Status::CRITICAL, "'CurLkToW' error, scale spline is not translation spline!!!");
        }
        double timeByBr = timeByLk + _parMagr->TEMPORAL.TO_LkToBr.at(topic);
        const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        const auto &posSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

        if (!so3Spline.TimeStampInRange(timeByBr) || !posSpline.TimeStampInRange(timeByBr)) {
            return {};
        } else {
            Sophus::SE3d curBrToW(so3Spline.Evaluate(timeByBr), posSpline.Evaluate(timeByBr));
            return curBrToW * _parMagr->EXTRI.SE3_LkToBr(topic);
        }
    }

    std::optional<Sophus::SE3d> CalibSolver::CurCmToW(double timeByCm, const std::string &topic) {
        if (GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
            throw Status(Status::CRITICAL, "'CurCmToW' error, scale spline is not translation spline!!!");
        }
        double timeByBr = timeByCm + _parMagr->TEMPORAL.TO_CmToBr.at(topic);
        const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        const auto &posSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

        if (!so3Spline.TimeStampInRange(timeByBr) || !posSpline.TimeStampInRange(timeByBr)) {
            return {};
        } else {
            Sophus::SE3d curBrToW(so3Spline.Evaluate(timeByBr), posSpline.Evaluate(timeByBr));
            return curBrToW * _parMagr->EXTRI.SE3_CmToBr(topic);
        }
    }

    std::optional<Sophus::SE3d> CalibSolver::CurRjToW(double timeByRj, const std::string &topic) {
        if (GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
            throw Status(Status::CRITICAL, "'CurRjToW' error, scale spline is not translation spline!!!");
        }
        double timeByBr = timeByRj + _parMagr->TEMPORAL.TO_RjToBr.at(topic);
        const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        const auto &posSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

        if (!so3Spline.TimeStampInRange(timeByBr) || !posSpline.TimeStampInRange(timeByBr)) {
            return {};
        } else {
            Sophus::SE3d curBrToW(so3Spline.Evaluate(timeByBr), posSpline.Evaluate(timeByBr));
            return curBrToW * _parMagr->EXTRI.SE3_RjToBr(topic);
        }
    }

    TimeDeriv::ScaleSplineType CalibSolver::GetScaleType() {
        if (Configor::IsLiDARIntegrated() || Configor::IsCameraIntegrated()) {
            return TimeDeriv::ScaleSplineType::LIN_POS_SPLINE;
        } else if (Configor::IsRadarIntegrated()) {
            return TimeDeriv::ScaleSplineType::LIN_VEL_SPLINE;
        } else {
            return TimeDeriv::ScaleSplineType::LIN_ACCE_SPLINE;
        }
    }

    CalibSolver::SplineBundleType::Ptr
    CalibSolver::CreateSplineBundle(double st, double et, double so3Dt, double scaleDt) {
        // create splines
        auto so3SplineInfo = ns_ctraj::SplineInfo(
                Configor::Preference::SO3_SPLINE, ns_ctraj::SplineType::So3Spline, st, et, so3Dt
        );
        auto scaleSplineInfo = ns_ctraj::SplineInfo(
                Configor::Preference::SCALE_SPLINE, ns_ctraj::SplineType::RdSpline, st, et, scaleDt
        );
        spdlog::info(
                "create spline bundle: start time: '{:.5f}', end time: '{:.5f}', so3 dt : '{:.5f}', scale dt: '{:.5f}'",
                st, et, so3Dt, scaleDt
        );
        return SplineBundleType::Create({so3SplineInfo, scaleSplineInfo});
    }

    void CalibSolver::AlignStatesToGravity() {
        auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
        // current gravity, velocities, and rotations are expressed in the reference frame
        // align them to the world frame whose negative z axis is aligned with the gravity vector
        auto SO3_RefToW = ObtainAlignedWtoRef(so3Spline.Evaluate(so3Spline.MinTime()), _parMagr->GRAVITY).inverse();
        _parMagr->GRAVITY = SO3_RefToW * _parMagr->GRAVITY;
        for (int i = 0; i < static_cast<int>(so3Spline.GetKnots().size()); ++i) {
            so3Spline.GetKnot(i) = SO3_RefToW * so3Spline.GetKnot(i);
        }
        for (int i = 0; i < static_cast<int>(scaleSpline.GetKnots().size()); ++i) {
            // attention: for three kinds of scale splines, this holds
            scaleSpline.GetKnot(i) = SO3_RefToW * scaleSpline.GetKnot(i);
        }
    }

    void CalibSolver::StoreImagesForSfM(const std::string &topic, const std::set<IndexPair> &matchRes) {

        // -------------
        // output images
        // -------------
        auto path = ns_ikalibr::Configor::DataStream::CreateImageStoreFolder(topic);
        if (path == std::nullopt) {
            throw ns_ikalibr::Status(
                    Status::CRITICAL,
                    "can not create path for image storing for topic: '{}'!!!", topic
            );
        }
        ImagesInfo info(topic, *path, {});
        const auto &frames = _dataMagr->GetCameraMeasurements(topic);

        int size = static_cast<int>(frames.size());
        const auto &intri = _parMagr->INTRI.Camera.at(topic);

        auto bar = std::make_shared<tqdm>();
        for (int i = 0; i != size; ++i) {
            bar->progress(i, size);
            const auto &frame = frames.at(i);
            // generate the image name
            std::string filename = std::to_string(frame->GetId()) + ".jpg";
            info.images[frame->GetId()] = filename;

            cv::Mat undistImg = CalibParamManager::ParIntri::UndistortImage(intri, frame->GetImage());

            // save image
            cv::imwrite(*path + "/" + filename, undistImg);
        }
        bar->finish();

        // -------------------
        // colmap command line
        // -------------------
        auto ws = ns_ikalibr::Configor::DataStream::CreateSfMWorkspace(topic);
        if (ws == std::nullopt) {
            throw ns_ikalibr::Status(
                    Status::CRITICAL,
                    "can not create workspace for SfM for topic: '{}'!!!", topic
            );
        }
        const std::string database_path = *ws + "/database.db";
        const std::string image_path = *path;
        const std::string match_list_path = *ws + "/matches.txt";
        const std::string output_path = *ws;

        auto logger = spdlog::basic_logger_mt("sfm_cmd", *ws + "/sfm-command-line.txt", true);
        // feature extractor
        logger->info("command line for 'feature_extractor' in colmap for topic '{}':\n"
                     "colmap feature_extractor "
                     "--database_path {} "
                     "--image_path {} "
                     "--ImageReader.camera_model PINHOLE "
                     "--ImageReader.single_camera 1 "
                     "--ImageReader.camera_params {:.3f},{:.3f},{:.3f},{:.3f}\n",
                     topic, database_path, image_path, intri->FocalX(), intri->FocalY(),
                     intri->PrincipalPoint()(0), intri->PrincipalPoint()(1));

        // feature match
        std::ofstream matchPairFile(match_list_path, std::ios::out);
        for (const auto &[view1Id, view2Id]: matchRes) {
            matchPairFile << std::to_string(view1Id) + ".jpg ";
            matchPairFile << std::to_string(view2Id) + ".jpg" << std::endl;
        }
        matchPairFile.close();

        logger->info("command line for 'matches_importer' in colmap for topic '{}':\n"
                     "colmap matches_importer "
                     "--database_path {} "
                     "--match_list_path {} "
                     "--match_type pairs\n",
                     topic, database_path, match_list_path);

        logger->info("------------------------------------------------------------------------------");
        logger->info("-  SfM Reconstruction in COLMAP [colmap gui] (recommend) or [colmap mapper]  -");
        logger->info("------------------------------------------------------------------------------");
        logger->info("performing SfM using [colmap gui] is suggested, rather than the command line, "
                     "which is very strict in initialization (finding initial image pair) and would cost lots of time!!!");
        // reconstruction
        logger->info("command line for 'feature_extractor' in [colmap gui] for topic '{}':\n"
                     "colmap gui "
                     "--database_path {} "
                     "--image_path {}",
                     topic, database_path, image_path, output_path);
        logger->info("------------------------------------------------------------------------------");
        double init_max_error = IsRSCamera(topic) ? 2.0 : 1.0;
        // reconstruction
        logger->info("command line for 'feature_extractor' in [colmap mapper] for topic '{}':\n"
                     "colmap mapper "
                     "--database_path {} "
                     "--image_path {} "
                     "--output_path {} "
                     "--Mapper.init_min_tri_angle 25 "
                     "--Mapper.init_max_error {} "
                     "--Mapper.tri_min_angle 3 "
                     "--Mapper.ba_refine_focal_length 0 "
                     "--Mapper.ba_refine_principal_point 0",
                     topic, database_path, image_path, output_path, init_max_error);
        logger->info("------------------------------------------------------------------------------\n");

        // format convert
        logger->info("command line for 'model_converter' in colmap for topic '{}':\n"
                     "colmap model_converter "
                     "--input_path {} "
                     "--output_path {} "
                     "--output_type TXT\n",
                     topic, output_path + "/0", output_path);
        logger->flush();
        spdlog::drop("sfm_cmd");

        std::ofstream file(ns_ikalibr::Configor::DataStream::GetImageStoreInfoFile(topic));
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::DataIOFormat());
        SerializeByOutputArchiveVariant(ar, Configor::Preference::DataIOFormat(), cereal::make_nvp("info", info));
    }

    ns_veta::Veta::Ptr CalibSolver::TryLoadSfMData(const std::string &topic, double errorThd, std::size_t trackLenThd) {
        // info file
        const auto infoFilename = ns_ikalibr::Configor::DataStream::GetImageStoreInfoFile(topic);
        if (!std::filesystem::exists(infoFilename)) {
            spdlog::warn("the info file, i.e., '{}', dose not exists!!!", infoFilename);
            return nullptr;
        }

        const auto &sfmWsPath = ns_ikalibr::Configor::DataStream::CreateSfMWorkspace(topic);
        if (!sfmWsPath) {
            spdlog::warn("the sfm workspace for topic '{}' dose not exists!!!", topic);
            return nullptr;
        }

        const auto camerasFilename = *sfmWsPath + "/cameras.txt";
        if (!std::filesystem::exists(camerasFilename)) {
            spdlog::warn("the cameras file, i.e., '{}', dose not exists!!!", camerasFilename);
            return nullptr;
        }

        // images
        const auto imagesFilename = *sfmWsPath + "/images.txt";
        if (!std::filesystem::exists(imagesFilename)) {
            spdlog::warn("the images file, i.e., '{}', dose not exists!!!", imagesFilename);
            return nullptr;
        }

        // points
        const auto ptsFilename = *sfmWsPath + "/points3D.txt";
        if (!std::filesystem::exists(ptsFilename)) {
            spdlog::warn("the points 3D file, i.e., '{}', dose not exists!!!", ptsFilename);
            return nullptr;
        }

        // cameras

        // load info file
        ImagesInfo info("", "", {});
        {
            std::ifstream file(infoFilename);
            auto ar = GetInputArchiveVariant(file, Configor::Preference::DataIOFormat());
            SerializeByInputArchiveVariant(ar, Configor::Preference::DataIOFormat(), cereal::make_nvp("info", info));
        }

        // load cameras
        auto cameras = ColMapDataIO::ReadCamerasText(camerasFilename);

        // load images
        auto images = ColMapDataIO::ReadImagesText(imagesFilename);

        // load landmarks
        auto points3D = ColMapDataIO::ReadPoints3DText(ptsFilename);

        auto veta = ns_veta::Veta::Create();

        // cameras
        assert(cameras.size() == 1);
        const auto &camera = cameras.cbegin()->second;
        assert(camera.params_.size() == 4);
        const auto &intriIdx = camera.camera_id_;
        auto intri = std::make_shared<ns_veta::PinholeIntrinsic>(*_parMagr->INTRI.Camera.at(topic));
        veta->intrinsics.insert({intriIdx, intri});

        // from images to our views and poses
        std::map<ns_veta::IndexT, CameraFrame::Ptr> ourIdxToCamFrame;
        for (const auto &frame: _dataMagr->GetCameraMeasurements(topic)) {
            ourIdxToCamFrame.insert({frame->GetId(), frame});
        }

        const auto &nameToOurIdx = info.GetImagesNameToIdx();
        for (const auto &[IdFromColmap, image]: images) {
            const auto &viewId = nameToOurIdx.at(image.name_);
            const auto &poseId = viewId;

            auto frameIter = ourIdxToCamFrame.find(viewId);
            // this frame is not involved in solving
            if (frameIter == ourIdxToCamFrame.cend()) { continue; }

            // view
            auto view = ns_veta::View::Create(
                    // timestamp (aligned)
                    frameIter->second->GetTimestamp(),
                    // index
                    viewId, intriIdx, poseId,
                    // width, height
                    intri->imgWidth, intri->imgHeight
            );
            veta->views.insert({viewId, view});

            // pose
            auto T_WorldToImg = ns_veta::Posed(image.QuatWorldToImg().matrix(), image.tvec_);
            // we store pose from camera to world
            veta->poses.insert({poseId, T_WorldToImg.Inverse()});
        }

        for (const auto &frame: _dataMagr->GetCameraMeasurements(topic)) {
            if (veta->views.count(frame->GetId()) == 0) {
                spdlog::warn(
                        "frame indexed as '{}' of camera '{}' is involved in solving but not reconstructed in SfM!!!",
                        frame->GetId(), topic);
            }
        }

        // from point3D to our structure
        for (const auto &[pt3dId, pt3d]: points3D) {
            // filter bad landmarks
            if (pt3d.error_ > errorThd || pt3d.track_.size() < trackLenThd) { continue; }

            auto &lm = veta->structure[pt3dId];
            lm.X = pt3d.xyz_;
            lm.color = pt3d.color_;

            for (const auto &track: pt3d.track_) {
                const auto &img = images.at(track.image_id);
                auto pt2d = img.points2D_.at(track.point2D_idx);

                if (pt3dId != pt2d.point3D_id_) {
                    spdlog::warn("'point3D_id_' of point3D and 'point3D_id_' of feature connected are in conflict!!!");
                    continue;
                }

                const auto viewId = nameToOurIdx.at(img.name_);
                // this frame is not involved in solving
                if (veta->views.find(viewId) == veta->views.cend()) { continue; }

                lm.obs.insert({viewId, ns_veta::Observation(pt2d.xy_, track.point2D_idx)});
            }
            if (lm.obs.size() < trackLenThd) { veta->structure.erase(pt3dId); }
        }

        return veta;
    }

    void CalibSolver::PerformTransformForVeta(const ns_veta::Veta::Ptr &veta,
                                              const ns_veta::Posed &curToNew, double scale) {
        // pose
        for (auto &[id, pose]: veta->poses) {
            pose.Translation() *= scale;
            pose = curToNew * pose;
        }

        // structure
        for (auto &[id, lm]: veta->structure) {
            lm.X *= scale;
            lm.X = curToNew(lm.X);
        }
    }

    bool CalibSolver::IsRSCamera(const std::string &camTopic) {
        auto model = EnumCast::stringToEnum<CameraModelType>(Configor::DataStream::CameraTopics.at(camTopic).Type);
        return IsOptionWith(CameraModelType::RS, model);
    }

    void CalibSolver::DownsampleVeta(const ns_veta::Veta::Ptr &veta, std::size_t lmNumThd, std::size_t obvNumThd) {
        std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
        if (veta->structure.size() > lmNumThd) {
            std::vector<ns_veta::IndexT> lmIdVec;
            lmIdVec.reserve(veta->structure.size());
            for (const auto &[lmId, lm]: veta->structure) { lmIdVec.push_back(lmId); }
            auto lmIdVecToMove = SamplingWoutReplace2(engine, lmIdVec, veta->structure.size() - lmNumThd);
            for (const auto &id: lmIdVecToMove) { veta->structure.erase(id); }
        }
        for (auto &[lmId, lm]: veta->structure) {
            if (lm.obs.size() < obvNumThd) { continue; }

            std::vector<ns_veta::IndexT> obvIdVec;
            obvIdVec.reserve(lm.obs.size());
            for (const auto &[viewId, obv]: lm.obs) { obvIdVec.push_back(viewId); }

            auto obvIdVecToMove = SamplingWoutReplace2(engine, obvIdVec, lm.obs.size() - obvNumThd);
            for (const auto &id: obvIdVecToMove) { lm.obs.erase(id); }
        }
    }

    void CalibSolver::SaveStageCalibParam(const CalibParamManager::Ptr &par, const std::string &desc) {
        const static std::string paramDir = Configor::DataStream::OutputPath + "/iteration/stage";
        if (!std::filesystem::exists(paramDir) && !std::filesystem::create_directories(paramDir)) {
            spdlog::warn("create directory failed: '{}'", paramDir);
        } else {
            const std::string paramFilename = paramDir + "/" + desc + ns_ikalibr::Configor::GetFormatExtension();
            par->Save(paramFilename, ns_ikalibr::Configor::Preference::DataIOFormat());
        }
    }

    // ------------------
    // CeresDebugCallBack
    // ------------------

    CeresDebugCallBack::CeresDebugCallBack(CalibParamManager::Ptr calibParamManager)
            : _parMagr(std::move(calibParamManager)),
              _outputDir(Configor::DataStream::OutputPath + "/iteration/epoch"), _idx(0) {
        if (std::filesystem::exists(_outputDir)) { std::filesystem::remove_all(_outputDir); }

        if (!std::filesystem::create_directories(_outputDir)) {
            spdlog::warn("create directory failed: '{}'", _outputDir);
        } else {
            _iterInfoFile = std::ofstream(_outputDir + "/epoch_info.csv", std::ios::out);
            _iterInfoFile << "cost,gradient,tr_radius(1/lambda)" << std::endl;
        }
    }

    ceres::CallbackReturnType CeresDebugCallBack::operator()(const ceres::IterationSummary &summary) {
        if (std::filesystem::exists(_outputDir)) {
            // save param
            const std::string paramFilename = _outputDir + "/ikalibr_param_" + std::to_string(_idx) +
                                              ns_ikalibr::Configor::GetFormatExtension();
            _parMagr->Save(paramFilename, ns_ikalibr::Configor::Preference::DataIOFormat());

            // save iter info
            _iterInfoFile << _idx << ','
                          << summary.cost << ','
                          << summary.gradient_norm << ','
                          << summary.trust_region_radius << std::endl;

            ++_idx;
        }
        return ceres::SOLVER_CONTINUE;
    }

    CeresDebugCallBack::~CeresDebugCallBack() {
        _iterInfoFile.close();
    }

    // -------------------
    // CeresViewerCallBack
    // -------------------
    CeresViewerCallBack::CeresViewerCallBack(Viewer::Ptr viewer) : _viewer(std::move(viewer)) {}

    ceres::CallbackReturnType CeresViewerCallBack::operator()(const ceres::IterationSummary &summary) {
        _viewer->UpdateSensorViewer().UpdateSplineViewer();
        return ceres::CallbackReturnType::SOLVER_CONTINUE;
    }
}