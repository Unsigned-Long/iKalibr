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

#include "calib/calib_data_manager.h"
#include "calib/calib_param_manager.h"
#include "calib/estimator.h"
#include "core/rotation_estimator.h"
#include "core/vision_only_sfm.h"
#include "opencv2/highgui.hpp"
#include "solver/calib_solver.h"
#include "util/tqdm.h"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepPosCameraInertialAlign() const {
    if (!Configor::IsPosCameraIntegrated()) {
        return;
    }
    const auto& so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto& scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    /**
     * we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are
     * poor
     */
    const double st = std::max(so3Spline.MinTime(), scaleSpline.MinTime()) +  // the max as start
                      Configor::Prior::TimeOffsetPadding;
    const double et = std::min(so3Spline.MaxTime(), scaleSpline.MaxTime()) -  // the min as end
                      Configor::Prior::TimeOffsetPadding;

    /**
     * the feature tracking is first performed for rotation-only visual odometry.
     * the rotation-only visual odometry means we only estimate the time-varying rotations of the
     * camera, which would be utilized for extrinsic rotation recovery
     */
    std::shared_ptr<tqdm> bar;
    std::map<std::string, RotOnlyVisualOdometer::Ptr> rotOnlyOdom;
    // how many features to maintain in each image
    constexpr int featNumPerImg = 300;
    // the min distance between two features (to ensure features are distributed uniformly)
    constexpr int minDist = 25;
    for (const auto& [topic, _] : Configor::DataStream::PosCameraTopics()) {
        const auto& frameVec = _dataMagr->GetCameraMeasurements(topic);
        spdlog::info(
            "perform rotation-only visual odometer to recover extrinsic rotations for '{}'...",
            topic);

        // estimates rotations
        auto intri = _parMagr->INTRI.Camera.at(topic);
        auto tracker = LKFeatureTracking::Create(featNumPerImg, minDist, intri);
        auto odometer = RotOnlyVisualOdometer::Create(tracker, intri);

        // sensor-inertial rotation estimator (linear least-squares problem)
        const auto rotEstimator = RotationEstimator::Create();

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

    /**
     * in last step, we use a discrete-time rotation-only hand-eye alignment, where only the
     * extrinsic rotations are considered, but here we use a continuous-time rotation-only
     * hand-eye alignment where both extrinsic rotations and temporal parameters are considered.
     */
    if (Configor::Prior::OptTemporalParams) {
        spdlog::info("perform rotation alignment to initialize time offset for each camera...");
        auto estimator = Estimator::Create(_splines, _parMagr);
        auto optOption = OptOption::OPT_SO3_CmToBr | OptOption::OPT_TO_CmToBr;

        // perform time offset estimation and extrinsic rotation refinement
        for (const auto& [topic, _] : Configor::DataStream::PosCameraTopics()) {
            const auto& rotations = rotOnlyOdom.at(topic)->GetRotations();
            // this field should be zero here
            double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(topic);
            double weight = Configor::DataStream::CameraTopics.at(topic).Weight;

            for (int i = 0; i < static_cast<int>(rotations.size()) - 1; ++i) {
                const auto &sRot = rotations.at(i), eRot = rotations.at(i + 1);
                // we throw the head and tail data as the rotations from the fitted SO3 Spline in
                // that range are poor
                if (sRot.first + TO_CmToBr < st || eRot.first + TO_CmToBr > et) {
                    continue;
                }

                estimator->AddHandEyeRotationAlignmentForCamera(
                    topic,        // the ros topic
                    sRot.first,   // the time of start rotation stamped by the camera
                    eRot.first,   // the time of end rotation stamped by the camera
                    sRot.second,  // the start rotation
                    eRot.second,  // the end rotation
                    optOption,    // the optimization option
                    weight        // the weight
                );
            }
        }

        auto sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    /**
     * perform SfM for each camera using rotation priori from the rotation spline and extrinsic
     * rotation. The SfM is performed using extern ColMap or GloMap. if SfM results are detected, we
     * load them into this program, otherwise, we output image sequence for SfM.
     */
    spdlog::info("perform SfM for each camera...");
    int needSfMCount = 0;
    for (const auto& [topic, _] : Configor::DataStream::PosCameraTopics()) {
        const auto& data = _dataMagr->GetCameraMeasurements(topic);
        // load data if SfM has been performed
        auto veta = TryLoadSfMData(
            // ros topic
            topic,
            // for rs camera, as the rs effect is not considered in SfM,
            // we relax the landmark selection condition
            IsRSCamera(topic) ? 2.0 : 1.0,
            // the track length threshold
            Configor::DataStream::CameraTopics.at(topic).TrackLengthMin);
        if (veta != nullptr) {
            /**
             * the SfM result data is valid fro this camera, we store it in the data manager
             */
            spdlog::info("SfM data for camera '{}' is valid!", topic);
            spdlog::info("down sample SfM data for camera '{}'", topic);

            // keep too many landmarks and features in estimator is not always good
            DownsampleVeta(
                // the visual meta data
                veta,
                // how many landmarks are maintained
                10000,
                // the track length threshold
                Configor::DataStream::CameraTopics.at(topic).TrackLengthMin);

            // we store SfM datas in '_dataMagr' as it is the 'calibration data manager'
            _dataMagr->SetSfMData(topic, veta);

            // just for visualization
            _viewer->AddVeta(veta, Viewer::VIEW_MAP);

            spdlog::info(
                "SfM info for topic '{}' after filtering: view count: {}, landmark count: {}",
                topic, veta->views.size(), veta->structure.size());
            continue;
        }

        /**
         * the SfM result data is invalid fro this camera, we output images for SfM
         */
        spdlog::warn("SfM data for camera '{}' is not valid! Perform SfM using colmap first!",
                     topic);

        /**
         * output images for SfM using thirdparty note that the rotation priors of each frame
         * from the extrinsic rotation and so3 spline are utilized in this process to accelerate
         * the feature matching
         */
        auto sfm = VisionOnlySfM::Create(topic, data, _parMagr, so3Spline, _viewer);
        /**
         * we use the thirdparty library, i.e., colmap, to solve the SfM problem.
         * output image frames and their corresponding information, which would be re-load after
         * SfM are performed.
         */
        spdlog::info("store images of '{}' for SfM...", topic);
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

    /**
     * the recovered pose from SfM are highly accurate, so we refine the extrinsic rotation and
     * time offsets again here, using the continuous-time sensor-inertial extrinsic rotation
     * recovery
     */
    /**
     * perhaps this is not too necessary. The difference from the last continuous-time
     * rotation-only hand-eye alignment is that we use the rotations recovered by colmap-derived
     * ones, rather than rotation-only odometer-derived ones.
     */
    spdlog::info(
        "perform rotation alignment to refine time offset for each camera using SfM "
        "results...");
    auto estimator = Estimator::Create(_splines, _parMagr);
    auto optOption = OptOption::OPT_SO3_CmToBr;
    if (Configor::Prior::OptTemporalParams) {
        optOption |= OptOption::OPT_TO_CmToBr;
    }

    for (const auto& [camTopic, veta] : _dataMagr->GetSfMData()) {
        double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(camTopic);
        double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;

        const auto& frames = _dataMagr->GetCameraMeasurements(camTopic);
        std::vector<ns_ctraj::Posed> constructedFrames;
        constructedFrames.reserve(frames.size());

        // find constructed frames by SfM
        for (const auto& frame : frames) {
            auto viewIter = veta->views.find(frame->GetId());
            if (viewIter == veta->views.cend()) {
                // this frame is not constructed (grabbed)
                continue;
            }
            auto poseIter = veta->poses.find(viewIter->second->poseId);
            if (poseIter == veta->poses.cend()) {
                // this frame is not constructed (grabbed)
                continue;
            }
            const auto& pose = poseIter->second;
            constructedFrames.emplace_back(pose.Rotation(), pose.Translation(),
                                           frame->GetTimestamp());
        }

        for (int i = 0; i < static_cast<int>(constructedFrames.size()) - 1; ++i) {
            const auto& sPose = constructedFrames.at(i);
            const auto& ePose = constructedFrames.at(i + 1);
            // we throw the head and tail data as the rotations from the fitted SO3 Spline in
            // that range are poor
            if (sPose.timeStamp + TO_CmToBr < st || ePose.timeStamp + TO_CmToBr > et) {
                continue;
            }

            estimator->AddHandEyeRotationAlignmentForCamera(
                camTopic,         // the ros topic
                sPose.timeStamp,  // the time of start pose stamped by the camera
                ePose.timeStamp,  // the time of end pose stamped by the camera
                sPose.so3,        // the start rotation
                ePose.so3,        // the end rotation
                optOption,        // the optimization option
                weight            // the weight
            );
        }
    }

    auto sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
}
}  // namespace ns_ikalibr