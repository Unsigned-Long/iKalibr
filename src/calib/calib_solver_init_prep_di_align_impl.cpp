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
#include "opencv2/highgui.hpp"
#include "core/rotation_estimator.h"
#include "util/tqdm.h"
#include "core/optical_flow_trace.h"
#include "core/visual_velocity_sac.h"
#include "factor/rgbd_velocity_factor.hpp"
#include "tiny-viewer/object/camera.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void CalibSolver::InitPrepRGBDInertialAlign() {
    if (!Configor::IsRGBDIntegrated()) {
        return;
    }
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    // we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are
    // poor
    const double st =
        std::max(so3Spline.MinTime(), scaleSpline.MinTime()) + Configor::Prior::TimeOffsetPadding;
    const double et =
        std::min(so3Spline.MaxTime(), scaleSpline.MaxTime()) - Configor::Prior::TimeOffsetPadding;

    // ---------------------------------------------------------------------------
    // perform rotation-only visual odometer to recover rotations and time offsets
    // ---------------------------------------------------------------------------
    // 'FeatTrackingInfo' is a list for each rgbd camera, as fail tracking leads to multiple pieces
    std::map<std::string, std::list<RotOnlyVisualOdometer::FeatTrackingInfo>> RGBDTrackingInfo;
    {
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

            auto bar = std::make_shared<tqdm>();
            auto intri = _parMagr->INTRI.RGBD.at(topic);
            for (int i = 0; i < static_cast<int>(frameVec.size()); ++i) {
                bar->progress(i, static_cast<int>(frameVec.size()));

                if (i % 30 == 0) {
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
                        auto estimator = Estimator::Create(_splines, _parMagr);

                        auto optOption =
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
        _viewer->ClearViewer(Viewer::VIEW_MAP);
        cv::destroyAllWindows();
    }

    // --------------------------------------------------
    // estimate rgbd-derived body-frame linear velocities
    // --------------------------------------------------
    for (const auto &[topic, trackInfoList] : RGBDTrackingInfo) {
        // store
        _dataMagr->SetRGBDOpticalFlowTrace(topic,
                                           CreateOpticalFlowTraceForRGBD(trackInfoList, topic));
    }
    // topic, camera frame, body-frame velocity
    auto &rgbdBodyFrameVels = _initAsset->rgbdBodyFrameVels;
    for (const auto &[topic, traceVec] : _dataMagr->GetRGBDOpticalFlowTrace()) {
        spdlog::info("estimate RGBD-derived linear velocities for '{}'...", topic);
        const auto &intri = _parMagr->INTRI.RGBD.at(topic);
        const auto &rsExposureFactor =
            CameraModel::RSCameraExposureFactor(EnumCast::stringToEnum<CameraModelType>(
                Configor::DataStream::RGBDTopics.at(topic).Type));

        // reorganize rgbd-traceVec, store them by frame index
        const auto &readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);
        std::map<CameraFrame::Ptr, std::vector<OpticalFlowCorr::Ptr>> dynamicsInFrame;
        for (const auto &trace : traceVec) {
            auto midCamFrame = trace->GetMidCameraFrame();
            // we use actual depth here (intrinsics is not nullptr)
            const auto &corr = trace->CreateOpticalFlowCorr(rsExposureFactor, intri);
            if (corr->depth < 1E-3 /* 1 mm */) {
                continue;
            }
            // a valid depth
            dynamicsInFrame[midCamFrame].emplace_back(corr);

            // if (Eigen::Vector2d vel = corr->MidPointVel(readout);
            //     vel.norm() > Configor::Prior::LossForRGBDFactor) {
            //     // show the visual pixel trace image (features, mid-point pixel velocity)
            //     auto img = trace->CreateOpticalFlowMat(_parMagr->INTRI.RGBD.at(topic)->intri,
            //                                               corr->MidPointVel(readout));
            //     Eigen::Vector2d mp = corr->MidPoint();
            //     const auto filename =
            //         fmt::format("{}/{}-{}-{}.png", Configor::DataStream::DebugPath,
            //                     corr->frame->GetId(), int(mp(0)), int(mp(1)));
            //     cv::imwrite(filename, img);
            //     cv::imshow("img", img);
            //     cv::waitKey(0);
            // }
        }

        // estimate rgbd-derived linear velocities for each frame
        const auto &rgbdIntri = _parMagr->INTRI.RGBD.at(topic);
        const double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(topic);
        const Sophus::SO3d &SO3_DnToBr = _parMagr->EXTRI.SO3_DnToBr.at(topic);
        for (const auto &[frame, corrVec] : dynamicsInFrame) {
            const double timeByBr = frame->GetTimestamp() + TO_DnToBr;
            // at least two measurements are required, here we up the ante
            if (timeByBr < st || timeByBr > et || corrVec.size() < 5) {
                continue;
            }
            auto res = VisualVelocitySacProblem::VisualVelocityEstimationRANSAC(
                corrVec, readout, rgbdIntri->intri, timeByBr, so3Spline, SO3_DnToBr);
            if (res) {
                rgbdBodyFrameVels[topic].emplace_back(frame, *res);
                // auto img = VisualVelocityEstimator::DrawVisualVelocityMat(
                //     curDynamics, rgbdIntri->intri, timeByBr, so3Spline, SO3_DnToBr, *res, frame,
                //     0.25);
            }
        }
        // sort timestamps
        auto &curRGBDVels = rgbdBodyFrameVels.at(topic);
        std::sort(curRGBDVels.begin(), curRGBDVels.end(), [](const auto &p1, const auto &p2) {
            return p1.first->GetTimestamp() < p2.first->GetTimestamp();
        });
    }
}
}  // namespace ns_ikalibr