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
#include "core/lidar_odometer.h"
#include "solver/calib_solver.h"
#include "util/utils_tpl.hpp"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void CalibSolver::InitSensorInertialAlign() const {
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

    spdlog::info("performing inertial alignment to initialize other spatial parameters...");

    /**
     * the gravity vector would be recovered in this stage, for better converage performance, we
     * assign the gravity roughly, f = a - g, g = a - f
     */
    Eigen::Vector3d firRefAcce =
        _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU).front()->GetAcce();
    // g = gDir * gNorm, where gDir = normalize(a - f), by assume the acceleration is zero
    _parMagr->GRAVITY = -firRefAcce.normalized() * Configor::Prior::GravityNorm;
    spdlog::info("rough assigned gravity in world frame: ['{:.3f}', '{:.3f}', '{:.3f}']",
                 _parMagr->GRAVITY(0), _parMagr->GRAVITY(1), _parMagr->GRAVITY(2));

    auto estimator = Estimator::Create(_splines, _parMagr);

    /**
     * we do not optimization the already initialized extrinsic rotations (IMUs', Cameras', and
     * LiDARs') here
     */
    auto optOption =
        // lidar extrinsic translations
        OptOption::OPT_POS_LkInBr |
        // camera extrinsic translations and visual scale
        OptOption::OPT_POS_CmInBr | OptOption::OPT_VISUAL_GLOBAL_SCALE |
        // radar extrinsics
        OptOption::OPT_SO3_RjToBr | OptOption::OPT_POS_RjInBr |
        // rgbd extrinsics
        OptOption::OPT_POS_DnInBr |
        // imu extrinsic translations and gravity
        OptOption::OPT_POS_BiInBr | OptOption::OPT_GRAVITY;

    // make sure 'ALIGN_STEP' * 'TIME INTERVAL OF TWO FRAMES' == 0.1
    static constexpr double DESIRED_TIME_INTERVAL = 0.5 /* 0.5 sed */;
    static constexpr double MIN_ALIGN_TIME = 1E-3 /* 0.001 sed */;
    static constexpr double MAX_ALIGN_TIME = 1.0 /* 1.0 sed */;

    // lidar-inertial alignment
    std::map<std::string, std::vector<Eigen::Vector3d>> linVelSeqLk;
    for (const auto &[lidarTopic, odometer] : _initAsset->lidarOdometers) {
        const auto &poseSeq = odometer->GetOdomPoseVec();
        double TO_LkToBr = _parMagr->TEMPORAL.TO_LkToBr.at(lidarTopic);

        // create linear velocity sequence
        linVelSeqLk[lidarTopic] =
            std::vector<Eigen::Vector3d>(poseSeq.size(), Eigen::Vector3d::Zero());
        auto &curLidarLinVelSeq = linVelSeqLk.at(lidarTopic);
        double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;

        const auto &refIMUFrames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);

        const int ALIGN_STEP =
            std::max(1, int(DESIRED_TIME_INTERVAL * _dataMagr->GetLiDARAvgFrequency(lidarTopic)));

        spdlog::info("add lidar-inertial alignment factors for '{}' and '{}', align step: {}",
                     lidarTopic, Configor::DataStream::ReferIMU);

        for (int i = 0; i < static_cast<int>(poseSeq.size()) - ALIGN_STEP; ++i) {
            const auto &sPose = poseSeq.at(i), ePose = poseSeq.at(i + ALIGN_STEP);
            if (sPose.timeStamp + TO_LkToBr < st || ePose.timeStamp + TO_LkToBr > et) {
                continue;
            }

            if (ePose.timeStamp - sPose.timeStamp < MIN_ALIGN_TIME ||
                ePose.timeStamp - sPose.timeStamp > MAX_ALIGN_TIME) {
                continue;
            }

            estimator->AddLiDARInertialAlignment(
                refIMUFrames,                           // the imu frames
                lidarTopic,                             // the ros topic of the lidar
                Configor::DataStream::ReferIMU,         // the ros topic of the imu
                sPose,                                  // the start pose
                ePose,                                  // the end pose
                odometer->GetMapTime(),                 // the map time
                &curLidarLinVelSeq.at(i),               // the start velocity (to be estimated)
                &curLidarLinVelSeq.at(i + ALIGN_STEP),  // the end velocity (to be estimated)
                optOption,                              // the optimize option
                weight);                                // the weigh
        }
    }

    // visual-inertial alignment
    std::map<std::string, std::vector<Eigen::Vector3d>> linVelSeqCm;
    std::map<std::string, double> visualScaleSeq;
    auto &sfmPoseSeq = _initAsset->sfmPoseSeq;
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

        const int ALIGN_STEP =
            std::max(1, int(DESIRED_TIME_INTERVAL * _dataMagr->GetCameraAvgFrequency(camTopic)));

        spdlog::info("add visual-inertial alignment factors for '{}' and '{}', align step: {}",
                     camTopic, Configor::DataStream::ReferIMU, ALIGN_STEP);

        for (int i = 0; i < static_cast<int>(constructedFrames.size()) - ALIGN_STEP; ++i) {
            const auto &sPose = constructedFrames.at(i);
            const auto &ePose = constructedFrames.at(i + ALIGN_STEP);
            // we throw the head and tail data as the rotations from the fitted SO3 Spline in
            // that range are poor
            if (sPose.timeStamp + TO_CmToBr < st || ePose.timeStamp + TO_CmToBr > et) {
                continue;
            }

            if (ePose.timeStamp - sPose.timeStamp < MIN_ALIGN_TIME ||
                ePose.timeStamp - sPose.timeStamp > MAX_ALIGN_TIME) {
                continue;
            }

            estimator->AddVisualInertialAlignment(
                imuFrames,                            // the imu frames
                camTopic,                             // the ros topic of the camera
                Configor::DataStream::ReferIMU,       // the ros topic of the imu
                sPose,                                // the start pose
                ePose,                                // the end pose
                firCTime,                             // the map time
                &curCamLinVelSeq.at(i),               // the start velocity (to be estimated)
                &curCamLinVelSeq.at(i + ALIGN_STEP),  // the end velocity (to be estimated)
                &scale,                               // the visual scale (to be estimated)
                optOption,                            // the optimize option
                weight);                              // the weigh
        }
    }

    // inertial alignment (only when more than or equal to 2 num IMUs are involved)
    constexpr double dt = DESIRED_TIME_INTERVAL;
    std::vector<Eigen::Vector3d> linVelSeqBr(std::floor((et - st) / dt), Eigen::Vector3d::Zero());
    if (Configor::DataStream::IMUTopics.size() >= 2) {
        for (const auto &[topic, frames] : _dataMagr->GetIMUMeasurements()) {
            spdlog::info("add inertial alignment factors for '{}'...", topic);

            for (int i = 0; i < static_cast<int>(linVelSeqBr.size()) - 1; ++i) {
                int sIdx = i, eIdx = i + 1;
                double sTimeByBr = sIdx * dt + st, eTimeByBr = eIdx * dt + st;
                Eigen::Vector3d *sVel = &linVelSeqBr.at(sIdx), *eVel = &linVelSeqBr.at(eIdx);

                estimator->AddInertialAlignment(
                    frames,     // imu frames
                    topic,      // the ros topic of this imu
                    sTimeByBr,  // the start time stamped by the reference imu
                    eTimeByBr,  // the end time stamped by the reference imu
                    sVel,       // the start velocity (to be estimated)
                    eVel,       // the end velocity (to be estimated)
                    optOption,  // the optimize option
                    Configor::DataStream::IMUTopics.at(topic).AcceWeight);
            }
        }
    }

    // radar-inertial alignment
    for (const auto &[radarTopic, radarMes] : _dataMagr->GetRadarMeasurements()) {
        double weight = Configor::DataStream::RadarTopics.at(radarTopic).Weight;
        double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(radarTopic);

        const auto &refIMUFrames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);

        const int ALIGN_STEP =
            std::max(1, int(DESIRED_TIME_INTERVAL * _dataMagr->GetRadarAvgFrequency(radarTopic)));

        spdlog::info("add radar-inertial alignment factors for '{}' and '{}', align step: {}",
                     radarTopic, Configor::DataStream::ReferIMU, ALIGN_STEP);

        for (int i = 0; i < static_cast<int>(radarMes.size()) - ALIGN_STEP; ++i) {
            const auto &sArray = radarMes.at(i), eArray = radarMes.at(i + ALIGN_STEP);

            /**
             * to estimate the radar velocity by linear least-squares solver the minim targets
             * number required is 3
             */
            if (sArray->GetTargets().size() < 10 || eArray->GetTargets().size() < 10) {
                continue;
            }

            if (sArray->GetTimestamp() + TO_RjToBr < st ||
                eArray->GetTimestamp() + TO_RjToBr > et) {
                continue;
            }

            if (eArray->GetTimestamp() - sArray->GetTimestamp() < MIN_ALIGN_TIME ||
                eArray->GetTimestamp() - sArray->GetTimestamp() > MAX_ALIGN_TIME) {
                continue;
            }

            /**
             * here, we use 'AddRadarInertialRotRoughAlignment', rather than
             * 'AddRadarInertialAlignment'. the difference is that:
             * 'AddRadarInertialAlignment': align both the extrinsic rotation and translation
             * 'AddRadarInertialRotRoughAlignment': only align the extrinsic rotation (roughly)
             *
             * for better converage performance, we use 'AddRadarInertialRotRoughAlignment' here,
             * after the extrinsic rotation is recovered, we refine it and estimate the translation
             */
            estimator->AddRadarInertialRotRoughAlignment(
                refIMUFrames,                    // imu frames
                Configor::DataStream::ReferIMU,  // the ros topic of this imu
                radarTopic,                      // the ros topic of this radar
                sArray,                          // the start target array
                eArray,                          // the end target array
                optOption,                       // the optimization option
                weight);                         // the weight
        }
    }

    // rgbd-inertial alignment
    for (const auto &[rgbdTopic, bodyFrameVels] : _initAsset->rgbdBodyFrameVels) {
        double weight = Configor::DataStream::RGBDTopics.at(rgbdTopic).Weight;
        double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(rgbdTopic);

        const auto &frames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);

        const int ALIGN_STEP =
            std::max(1, int(DESIRED_TIME_INTERVAL * _dataMagr->GetRGBDAvgFrequency(rgbdTopic)));

        spdlog::info("add rgbd-inertial alignment factors for '{}' and '{}', align step: {}",
                     rgbdTopic, Configor::DataStream::ReferIMU, ALIGN_STEP);

        for (int i = 0; i < static_cast<int>(bodyFrameVels.size()) - ALIGN_STEP; ++i) {
            const auto &[sFrame, sVel] = bodyFrameVels.at(i);
            const auto &[eFrame, eVel] = bodyFrameVels.at(i + ALIGN_STEP);

            if (sFrame->GetTimestamp() + TO_DnToBr < st ||
                eFrame->GetTimestamp() + TO_DnToBr > et) {
                continue;
            }

            if (eFrame->GetTimestamp() - sFrame->GetTimestamp() < MIN_ALIGN_TIME ||
                eFrame->GetTimestamp() - sFrame->GetTimestamp() > MAX_ALIGN_TIME) {
                continue;
            }

            // we don't consider the readout time here (i.e., the camera frame time is treated as
            // the velocity of the midpoint)
            estimator->AddRGBDInertialAlignment(
                frames,                            // imu frames
                Configor::DataStream::ReferIMU,    // the ros topic of this imu
                rgbdTopic,                         // the ros topic of this rgbd camera
                bodyFrameVels.at(i),               // the start velocity
                bodyFrameVels.at(i + ALIGN_STEP),  // the end velocity
                optOption,                         // the optimization option
                weight);                           // the weight
        }
    }

    // (vel) camera-inertial alignment
    std::map<std::string, std::vector<double>> velCamLinVelScales;
    for (const auto &[topic, velDirs] : _initAsset->velCamBodyFrameVelDirs) {
        double weight = Configor::DataStream::CameraTopics.at(topic).Weight;
        double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(topic);
        const auto &frames = _dataMagr->GetIMUMeasurements(Configor::DataStream::ReferIMU);

        const int ALIGN_STEP =
            std::max(1, int(DESIRED_TIME_INTERVAL * _dataMagr->GetCameraAvgFrequency(topic)));

        spdlog::info(
            "add vel-visual-inertial velocity alignment factors for '{}' and '{}', align step: {}",
            topic, Configor::DataStream::ReferIMU, ALIGN_STEP);

        auto &curVelScales = velCamLinVelScales[topic];
        curVelScales = std::vector<double>(velDirs.size(), -1.0);
        for (int i = 0; i < static_cast<int>(velDirs.size()) - ALIGN_STEP; ++i) {
            const auto &[sFrame, sVel] = velDirs.at(i);
            const auto &[eFrame, eVel] = velDirs.at(i + ALIGN_STEP);

            if (sFrame->GetTimestamp() + TO_CmToBr < st ||
                eFrame->GetTimestamp() + TO_CmToBr > et) {
                continue;
            }

            if (eFrame->GetTimestamp() - sFrame->GetTimestamp() < MIN_ALIGN_TIME ||
                eFrame->GetTimestamp() - sFrame->GetTimestamp() > MAX_ALIGN_TIME) {
                continue;
            }

            // create scales for valid range
            curVelScales.at(i) = 1.0;
            curVelScales.at(i + ALIGN_STEP) = 1.0;
            estimator->AddVelVisualInertialAlignment(
                frames,                            // imu frames
                Configor::DataStream::ReferIMU,    // ros topic of the imu
                topic,                             // ros topic of the camera
                velDirs.at(i),                     // the direction of the start velocity
                &curVelScales.at(i),               // the scale of start velocity to be estimated
                velDirs.at(i + ALIGN_STEP),        // the direction of the end velocity
                &curVelScales.at(i + ALIGN_STEP),  // the scale of end velocity to be estimated
                optOption,                         // the optimization option
                weight);                           // the weight
        }
    }

    // fix spatiotemporal parameters of reference sensor
    // make this problem full rank
    estimator->SetRefIMUParamsConstant();

    auto sum = estimator->Solve(_ceresOption, this->_priori);
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

                /**
                 * to estimate the radar velocity by linear least-squares solver the minim targets
                 * number required is 3
                 */
                if (sArray->GetTargets().size() < 10 || eArray->GetTargets().size() < 10) {
                    continue;
                }

                if (sArray->GetTimestamp() + TO_RjToBr < st ||
                    eArray->GetTimestamp() + TO_RjToBr > et) {
                    continue;
                }

                estimator->AddRadarInertialAlignment(
                    frames,                          // imu frames
                    Configor::DataStream::ReferIMU,  // the ros topic of this imu
                    radarTopic,                      // the ros topic of this radar
                    sArray,                          // the start target array
                    eArray,                          // the end target array
                    optOption,                       // the optimization option
                    weight);                         // the weight
            }
        }
        estimator->SetRefIMUParamsConstant();
        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }

    if (Configor::IsVelCameraIntegrated()) {
        for (const auto &[topic, velDirs] : _initAsset->velCamBodyFrameVelDirs) {
            const auto &curVelCamLinVelScales = velCamLinVelScales.at(topic);
            auto &curVels = _initAsset->velCamBodyFrameVels[topic];
            for (int i = 0; i < static_cast<int>(curVelCamLinVelScales.size()); ++i) {
                if (const double velScale = curVelCamLinVelScales.at(i); velScale > 1E-3) {
                    const auto &[camFrame, velDir] = velDirs.at(i);
                    curVels.emplace_back(camFrame, velScale * velDir);
                }
            }
        }
    }

    for (const auto &[topic, lidarOdom] : _initAsset->lidarOdometers) {
        _viewer->AddCloud(lidarOdom->GetMap(), Viewer::VIEW_MAP,
                          ns_viewer::Entity::GetUniqueColour(), 2.0f);
    }
    /**
     * based on the estimated visual scales from the sensor-inertial alignment, we update the vetas
     */
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

    // linVelSeqLk.clear(), linVelSeqBr.clear(), linVelSeqCm.clear(), visualScaleSeq.clear();
}
}  // namespace ns_ikalibr