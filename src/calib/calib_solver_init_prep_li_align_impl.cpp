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
#include "core/rotation_estimator.h"
#include "core/lidar_odometer.h"
#include "core/scan_undistortion.h"
#include "util/tqdm.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void CalibSolver::InitPrepLiDARInertialAlign() {
    if (!Configor::IsLiDARIntegrated()) {
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

    // -------------------------------------------------------------
    // initialize extrinsic rotations of LiDARs if they are involved
    // -------------------------------------------------------------
    {
        spdlog::info("LiDARs are integrated, initializing extrinsic rotations of LiDARs...");
        for (const auto &[topic, data] : _dataMagr->GetLiDARMeasurements()) {
            spdlog::info(
                "performing ndt odometer for '{}' for extrinsic rotation initialization...", topic);

            auto lidarOdometer = LiDAROdometer::Create(
                static_cast<float>(Configor::Prior::NDTLiDAROdometer::Resolution),
                Configor::Preference::AvailableThreads());
            auto rotEstimator = RotationEstimator::Create();
            auto bar = std::make_shared<tqdm>();
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
    // undistort lidar scans and rerun odometer
    // --------------------------------------------------
    auto &lidarOdometers = _initAsset->lidarOdometers;
    auto &undistFramesInScan = _initAsset->undistFramesInScan;
    {
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
            auto bar = std::make_shared<tqdm>();
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

    // ----------------------------------------------------------------
    // initialize time offsets of LiDARs by hand eye rotation alignment
    // ----------------------------------------------------------------
    {
        spdlog::info("performing  hand eye rotation alignment for LiDARs...");
        auto estimator = Estimator::Create(_splines, _parMagr);
        auto optOption = OptOption::Option::OPT_SO3_LkToBr;
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

        auto sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }
}
}  // namespace ns_ikalibr