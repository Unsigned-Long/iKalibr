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
#include "core/optical_flow_trace.h"
#include "core/rotation_estimator.h"
#include "core/visual_velocity_sac.h"
#include "factor/data_correspondence.h"
#include "opencv2/highgui.hpp"
#include "solver/calib_solver.h"
#include "spdlog/spdlog.h"
#include "tiny-viewer/object/camera.h"
#include "util/tqdm.h"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void CalibSolver::InitPrepRGBDInertialAlign() const {
    if (!Configor::IsRGBDIntegrated()) {
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
    // 'FeatTrackingInfo' is a list for each rgbd camera, as fail tracking leads to multiple pieces
    std::map<std::string, std::list<RotOnlyVisualOdometer::FeatTrackingInfo>> RGBDTrackingInfo;
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
        auto intri = _parMagr->INTRI.RGBD.at(topic);
        auto tracker = LKFeatureTracking::Create(featNumPerImg, minDist, intri->intri);
        auto odometer = RotOnlyVisualOdometer::Create(tracker, intri->intri);

        // sensor-inertial rotation estimator (linear least-squares problem)
        auto rotEstimator = RotationEstimator::Create();

        auto bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(frameVec.size()); ++i) {
            bar->progress(i, static_cast<int>(frameVec.size()));

            if (i % 30 == 0) {
                /**
                 * we do not update the viewer too frequent, which would lead to heavy tasks
                 */
                _viewer->ClearViewer(Viewer::VIEW_MAP);
                // rgbd camera
                static auto rgbd = ns_viewer::CubeCamera::Create(
                    ns_viewer::Posef(), 0.04, ns_viewer::Colour(1.0f, 0.5f, 0.0f, 1.0f));
                _viewer->AddEntityLocal({rgbd}, Viewer::VIEW_MAP);
                // depth point could
                _viewer->AddRGBDFrame(frameVec.at(i), intri, Viewer::VIEW_MAP, true, 2.0f);
                // auto img = frameVec.at(i)->CreateColorDepthMap(intri, true);
                // cv::imshow("img", img);
                // cv::waitKey();
            }

            /*
             * we try to compute the prior rotation to accelerate the feature tracking.
             * only the extrinsic rotation is recovered, we can compute such priori.
             */
            std::optional<Sophus::SO3d> SO3_LastToCur = std::nullopt;
            if (rotEstimator->SolveStatus()) {
                // such codes should be put out of the for loop, but fore better readability...
                const auto SO3_DnToBr = _parMagr->EXTRI.SO3_DnToBr.at(topic);
                const auto TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(topic);
                double lastTimeByBr = frameVec.at(i - 1)->GetTimestamp() + TO_DnToBr;
                double curTimeByBr = frameVec.at(i)->GetTimestamp() + TO_DnToBr;
                if (so3Spline.TimeStampInRange(lastTimeByBr) &&
                    so3Spline.TimeStampInRange(curTimeByBr)) {
                    // compute rotations at two timestamps
                    auto SO3_LastBrToW = so3Spline.Evaluate(lastTimeByBr);
                    auto SO3_CurBrToW = so3Spline.Evaluate(lastTimeByBr);
                    // compute the relative rotation of the reference imu
                    auto SO3_LastBrToCurBr = SO3_CurBrToW.inverse() * SO3_LastBrToW;
                    // assignment
                    SO3_LastToCur = SO3_DnToBr.inverse() * SO3_LastBrToCurBr * SO3_DnToBr;
                }
            }

            // if tracking current frame failed, the rotation-only odometer would re-initialize
            if (!odometer->GrabFrame(frameVec.at(i), SO3_LastToCur)) {
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

                /**
                 * once the extrinsic rotation is recovered, if time offset is also required, we
                 * continue to recover it and refine extrineic rotation using
                 * continuous-time-based alignment
                 */
                if (Configor::Prior::OptTemporalParams) {
                    auto estimator = Estimator::Create(_splines, _parMagr);

                    auto optOption = OptOption::OPT_SO3_DnToBr | OptOption::OPT_TO_DnToBr;
                    double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(topic);
                    double weight = Configor::DataStream::RGBDTopics.at(topic).Weight;

                    const auto &rotations = odometer->GetRotations();
                    for (int j = 0; j < static_cast<int>(rotations.size()) - 1; ++j) {
                        const auto &sRot = rotations.at(j), eRot = rotations.at(j + 1);
                        // we throw the head and tail data as the rotations from the fitted SO3
                        // Spline in that range are poor
                        if (sRot.first + TO_DnToBr < st || eRot.first + TO_DnToBr > et) {
                            continue;
                        }

                        estimator->AddHandEyeRotationAlignmentForRGBD(
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
        RGBDTrackingInfo[topic].push_back(odometer->GetLmTrackInfo());

        /**
         * after all images are grabbed, if the extrinsic rotation is not recovered (use min
         * eigen value to check solve results), stop this program
         */
        if (!rotEstimator->SolveStatus()) {
            throw Status(Status::ERROR,
                         "initialize rotation 'SO3_DnToBr' failed, this may be related to "
                         "insufficiently excited motion or bad images.");
        }
    }
    _viewer->ClearViewer(Viewer::VIEW_MAP);
    cv::destroyAllWindows();

    /**
     * based on the tracking information, we construct the optical flow trace to estimate pixel
     * velocity of tracked featuers
     */
    for (const auto &[topic, trackInfoList] : RGBDTrackingInfo) {
        const int trackThd = Configor::DataStream::RGBDTopics.at(topic).TrackLengthMin;
        // store
        _dataMagr->SetVisualOpticalFlowTrace(
            // ros topic of this rgbd camera
            topic,
            // create the optical flow trace
            CreateOpticalFlowTrace(trackInfoList, trackThd));
    }

    /**
     * resort the optical flow correspondences based on the camera frame index
     */
    std::map<std::string, std::map<CameraFrame::Ptr, std::vector<OpticalFlowCorr::Ptr>>>
        opticalFlowInFrame;
    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        const auto &traceVec = _dataMagr->GetVisualOpticalFlowTrace(topic);
        const auto &intri = _parMagr->INTRI.RGBD.at(topic);
        // the rs exposure factor to compute the real visual timestamps
        const auto &rsExpFactor =
            CameraModel::RSCameraExposureFactor(EnumCast::stringToEnum<CameraModelType>(
                Configor::DataStream::RGBDTopics.at(topic).Type));
        auto &curOpticalFlowInFrame = opticalFlowInFrame[topic];

        for (const auto &trace : traceVec) {
            auto midCamFrame = trace->GetMidCameraFrame();
            // we use actual depth here (intrinsics is not nullptr)
            const auto &corr = trace->CreateOpticalFlowCorr(rsExpFactor, intri);
            if (corr->depth < 1E-3 /* 1 mm */) {
                continue;
            }
            // a valid depth
            curOpticalFlowInFrame[midCamFrame].emplace_back(corr);

#define VISUALIZE_OPTICAL_FLOW_TRACE 0
#if VISUALIZE_OPTICAL_FLOW_TRACE
            if (Eigen::Vector2d vel = corr->MidPointVel(readout);  // the pixel velocity
                vel.norm() > 2.0 * Configor::Prior::LossForOpticalFlowFactor) {
                // show the visual pixel trace image (features, mid-point pixel velocity)
                auto img = trace->CreateOpticalFlowMat(_parMagr->INTRI.RGBD.at(topic)->intri,
                                                       corr->MidPointVel(readout));
                Eigen::Vector2d mp = corr->MidPoint();
                const auto filename = fmt::format(
                    "{}/{}-{}-{}.png", Configor::DataStream::DebugPath, corr->frame->GetId(),
                    static_cast<int>(mp(0)), static_cast<int>(mp(1)));
                // save this image to disk
                cv::imwrite(filename, img);
                cv::imshow("img", img);
                cv::waitKey(0);
            }
#endif
        }
    }

    /**
     * for each camera frame of each topic, we estimate the linear velocity, and store them in
     * 'rgbdBodyFrameVels' [topic, camera frame, body-frame velocity]
     */
    auto &rgbdBodyFrameVels = _initAsset->rgbdBodyFrameVels;
    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        spdlog::info("estimate RGBD-derived linear velocities for '{}'...", topic);
        const auto &readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);
        const auto &rgbdIntri = _parMagr->INTRI.RGBD.at(topic);
        const double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(topic);
        const Sophus::SO3d &SO3_DnToBr = _parMagr->EXTRI.SO3_DnToBr.at(topic);
        for (const auto &[frame, ofVec] : opticalFlowInFrame.at(topic)) {
            const double timeByBr = frame->GetTimestamp() + TO_DnToBr;
            // at least two measurements are required, here we up the ante
            if (timeByBr < st || timeByBr > et || ofVec.size() < 5) {
                continue;
            }
            auto res = VisualVelocitySacProblem::VisualVelocityEstimationRANSAC(
                ofVec,             // pixel velocity sequence in this image
                readout,           // the readout time of the associated camera
                rgbdIntri->intri,  // the visual intrinsics
                timeByBr,          // the time stamped by the reference imu
                so3Spline,         // the rotation spline
                SO3_DnToBr         // the extrinsic rotation
            );
            if (res) {
                rgbdBodyFrameVels[topic].emplace_back(frame, *res);
#define VISUALIZE_RGBD_ONLY_VEL_EST 0
#if VISUALIZE_RGBD_ONLY_VEL_EST
                // feature, velocity, depth
                std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>> dynamics;
                dynamics.reserve(ofVec.size());
                for (const auto &ofCorr : ofVec) {
                    if (ofCorr->depth < 1E-3 /* 1mm */) {
                        continue;
                    }
                    dynamics.emplace_back(
                        ofCorr->MidPoint(),            // feature position
                        ofCorr->MidPointVel(readout),  // pixel velocity of the middle point
                        ofCorr->depth                  // the depth of the middle point
                    );
                }
                auto img = VisualVelocityEstimator::DrawVisualVelocityMat(
                    dynamics, rgbdIntri->intri, timeByBr, so3Spline, SO3_DnToBr, *res, frame, 0.25);
#endif
            }
        }
        /**
         * the obtained rgbd-frame velocities may not time-ordered, thus we sort these quantities
         * based on their time stamps
         */
        auto &curRGBDVels = rgbdBodyFrameVels.at(topic);
        std::sort(curRGBDVels.begin(), curRGBDVels.end(), [](const auto &p1, const auto &p2) {
            return p1.first->GetTimestamp() < p2.first->GetTimestamp();
        });
    }
}
}  // namespace ns_ikalibr