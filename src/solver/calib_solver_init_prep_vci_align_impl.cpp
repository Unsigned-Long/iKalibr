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
#include "util/tqdm.h"
#include "spdlog/spdlog.h"
#include "calib/calib_data_manager.h"
#include "calib/calib_param_manager.h"
#include "viewer/viewer.h"
#include "opencv2/highgui.hpp"
#include "core/rotation_estimator.h"
#include "tiny-viewer/object/camera.h"
#include "calib/estimator.h"
#include "sensor/sensor_model.h"
#include "factor/data_correspondence.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepVelCameraInertialAlign() const {
    if (!Configor::IsVelCameraIntegrated()) {
        return;
    }
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
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
    // 'FeatTrackingInfo' is a list for each camera, as fail tracking leads to multiple pieces
    std::map<std::string, std::list<RotOnlyVisualOdometer::FeatTrackingInfo>> trackingInfo;
    // how many features to maintain in each image
    constexpr int featNumPerImg = 300;
    // the min distance between two features (to ensure features are distributed uniformly)
    constexpr int minDist = 25;
    for (const auto &[topic, _] : Configor::DataStream::VelCameraTopics()) {
        const auto &frameVec = _dataMagr->GetCameraMeasurements(topic);
        spdlog::info(
            "perform rotation-only visual odometer to recover extrinsic rotations for camera "
            "'{}'...",
            topic);

        // estimates rotations
        auto intri = _parMagr->INTRI.Camera.at(topic);
        auto tracker = LKFeatureTracking::Create(featNumPerImg, minDist, intri);
        auto odometer = RotOnlyVisualOdometer::Create(tracker, intri);

        // sensor-inertial rotation estimator (linear least-squares problem)
        auto rotEstimator = RotationEstimator::Create();

        auto bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(frameVec.size()); ++i) {
            bar->progress(i, static_cast<int>(frameVec.size()));

            /*
             * we try to compute the prior rotation to accelerate the feature tracking.
             * only the extrinsic rotation is recovered, we can compute such priori.
             */
            std::optional<Sophus::SO3d> SO3_LastToCur = std::nullopt;
            if (rotEstimator->SolveStatus()) {
                // such codes should be put out of the for loop, but fore better readability...
                const auto &SO3_CmToBr = _parMagr->EXTRI.SO3_CmToBr.at(topic);
                const auto &TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(topic);

                double lastTimeByBr = frameVec.at(i - 1)->GetTimestamp() + TO_CmToBr;
                double curTimeByBr = frameVec.at(i)->GetTimestamp() + TO_CmToBr;
                if (so3Spline.TimeStampInRange(lastTimeByBr) &&
                    so3Spline.TimeStampInRange(curTimeByBr)) {
                    // compute rotations at two timestamps
                    auto SO3_LastBrToW = so3Spline.Evaluate(lastTimeByBr);
                    auto SO3_CurBrToW = so3Spline.Evaluate(lastTimeByBr);
                    // compute the relative rotation of the reference imu
                    auto SO3_LastBrToCurBr = SO3_CurBrToW.inverse() * SO3_LastBrToW;
                    // assignment
                    SO3_LastToCur = SO3_CmToBr.inverse() * SO3_LastBrToCurBr * SO3_CmToBr;
                }
            }

            // if tracking current frame failed, the rotation-only odometer would re-initialize
            if (!odometer->GrabFrame(frameVec.at(i), SO3_LastToCur)) {
                spdlog::warn(
                    "tracking failed when grab the '{}' image frame!!! try to reinitialize", i);
                // save the tracking information
                trackingInfo[topic].push_back(odometer->GetLmTrackInfo());
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
                _parMagr->EXTRI.SO3_CmToBr.at(topic) = rotEstimator->GetSO3SensorToSpline();

                /**
                 * once the extrinsic rotation is recovered, if time offset is also required, we
                 * continue to recover it and refine extrineic rotation using
                 * continuous-time-based alignment
                 */
                if (Configor::Prior::OptTemporalParams) {
                    auto estimator = Estimator::Create(_splines, _parMagr);

                    auto optOption = OptOption::OPT_SO3_CmToBr | OptOption::OPT_TO_CmToBr;
                    double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(topic);
                    double weight = Configor::DataStream::CameraTopics.at(topic).Weight;

                    const auto &rotations = odometer->GetRotations();
                    for (int j = 0; j < static_cast<int>(rotations.size()) - 1; ++j) {
                        const auto &sRot = rotations.at(j), eRot = rotations.at(j + 1);
                        // we throw the head and tail data as the rotations from the fitted SO3
                        // Spline in that range are poor
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

                    // we don't want to output the solving information
                    auto optWithoutOutput = Estimator::DefaultSolverOptions(
                        Configor::Preference::AvailableThreads(),
                        false,  // do not output the solving information
                        Configor::Preference::UseCudaInSolving);

                    estimator->Solve(optWithoutOutput, _priori);
                }

                _viewer->UpdateSensorViewer();
            }
        }
        bar->finish();

        // add tracking info
        trackingInfo[topic].push_back(odometer->GetLmTrackInfo());

        /**
         * after all images are grabbed, if the extrinsic rotation is not recovered (use min
         * eigen value to check solve results), stop this program
         */
        if (!rotEstimator->SolveStatus()) {
            throw Status(Status::ERROR,
                         "initialize rotation 'SO3_CmToBr' failed, this may be related to "
                         "insufficiently excited motion or bad images.");
        }
    }
    _viewer->ClearViewer(Viewer::VIEW_MAP);
    cv::destroyAllWindows();

    /**
     * based on the tracking information, we construct the optical flow trace to estimate pixel
     * velocity of tracked featuers
     */
    for (const auto &[topic, trackInfoList] : trackingInfo) {
        const int trackThd = Configor::DataStream::CameraTopics.at(topic).TrackLengthMin;
        // store
        _dataMagr->SetVisualOpticalFlowTrace(
            // ros topic of this camera
            topic,
            // create the optical flow trace
            CreateOpticalFlowTrace(trackInfoList, trackThd));
    }

    spdlog::warn("this module is being developed!!!");
    std::cin.get();
}
}  // namespace ns_ikalibr