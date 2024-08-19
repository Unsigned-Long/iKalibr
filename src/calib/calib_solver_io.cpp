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
#include "calib/calib_solver_io.h"
#include "viewer/visual_lidar_covisibility.h"
#include "viewer/visual_ang_vel_drawer.h"
#include "viewer/visual_lin_vel_drawer.h"
#include "viewer/visual_gravity.h"
#include "viewer/visual_colorized_cloud_map.h"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "cereal/types/utility.hpp"
#include "cereal/types/list.hpp"
#include "util/tqdm.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

CalibSolverIO::CalibSolverIO(CalibSolver::Ptr solver)
    : _solver(std::move(solver)) {
    if (!_solver->_solveFinished) {
        spdlog::warn("calibration has not been performed!!! Do not try to save anything now!!!");
    }
}

CalibSolverIO::Ptr CalibSolverIO::Create(const CalibSolver::Ptr &solver) {
    return std::make_shared<CalibSolverIO>(solver);
}

void CalibSolverIO::SaveBSplines(int hz) {
    std::string saveDir = Configor::DataStream::OutputPath + "/splines";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving splines to dir: '{}'...", saveDir);
    } else {
        return;
    }
    {
        // sampled linear scale points
        const auto &so3Spline = _solver->_splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
        const auto &scaleSpline =
            _solver->_splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
        const double st = std::max(so3Spline.MinTime(), scaleSpline.MinTime());
        const double et = std::min(so3Spline.MaxTime(), scaleSpline.MaxTime());
        const double dt = 1.0 / hz;
        // pose container
        Eigen::aligned_vector<ns_ctraj::Posed> poseSeq;
        poseSeq.reserve(static_cast<std::size_t>((et - st) / dt));

        for (double t = st; t < et;) {
            if (_solver->_splines->TimeInRange(t, so3Spline) &&
                _solver->_splines->TimeInRange(t, scaleSpline)) {
                Sophus::SO3d curSO3_BrToBr0 = so3Spline.Evaluate(t);
                Eigen::Vector3d linAcce_BrToBr0InBr0 = scaleSpline.Evaluate(t);
                poseSeq.emplace_back(curSO3_BrToBr0, linAcce_BrToBr0InBr0, t);
            }
            t += dt;
        }
        auto filename = saveDir + "/samples" + Configor::GetFormatExtension();
        SavePoseSequence(poseSeq, filename, Configor::Preference::OutputDataFormat);
    }
    {
        // control points
        auto filename = saveDir + "/knots" + Configor::GetFormatExtension();
        std::ofstream file(filename);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(ar, Configor::Preference::OutputDataFormat,
                                        cereal::make_nvp("splines", *_solver->_splines));
    }
    spdlog::info("saving splines finished!");
}

void CalibSolverIO::SaveHessianMatrix() {
    std::string saveDir = Configor::DataStream::OutputPath + "/hessian";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving hessian matrix to dir: '{}'...", saveDir);
    } else {
        return;
    }

    std::vector<double *> parAddress;
    auto estimator = _solver->_backup->estimator;
    parAddress.reserve(estimator->NumParameterBlocks());
    std::vector<std::pair<std::string, int>> parOrderSize;
    parOrderSize.reserve(estimator->NumParameterBlocks());

    auto InvolveParameter = [&parAddress, &parOrderSize, &estimator](double *address,
                                                                     const std::string &name) {
        if (!estimator->HasParameterBlock(address)) {
            return;
        }
        parAddress.push_back(address);
        parOrderSize.emplace_back(name, estimator->ParameterBlockTangentSize(address));
    };

    // spatiotemporal parameters (extrinsics, time offsets, rs readout time)
    for (auto &[topic, item] : _solver->_parMagr->EXTRI.SO3_BiToBr) {
        InvolveParameter(item.data(), "SO3_BiToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->EXTRI.POS_BiInBr) {
        InvolveParameter(item.data(), "POS_BiInBr-" + topic);
    }

    for (auto &[topic, item] : _solver->_parMagr->EXTRI.SO3_RjToBr) {
        InvolveParameter(item.data(), "SO3_RjToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->EXTRI.POS_RjInBr) {
        InvolveParameter(item.data(), "POS_RjInBr-" + topic);
    }

    for (auto &[topic, item] : _solver->_parMagr->EXTRI.SO3_LkToBr) {
        InvolveParameter(item.data(), "SO3_LkToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->EXTRI.POS_LkInBr) {
        InvolveParameter(item.data(), "POS_LkInBr-" + topic);
    }

    for (auto &[topic, item] : _solver->_parMagr->EXTRI.SO3_CmToBr) {
        InvolveParameter(item.data(), "SO3_CmToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->EXTRI.POS_CmInBr) {
        InvolveParameter(item.data(), "POS_CmInBr-" + topic);
    }

    for (auto &[topic, item] : _solver->_parMagr->EXTRI.SO3_DnToBr) {
        InvolveParameter(item.data(), "SO3_DnToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->EXTRI.POS_DnInBr) {
        InvolveParameter(item.data(), "POS_DnInBr-" + topic);
    }

    for (auto &[topic, item] : _solver->_parMagr->TEMPORAL.TO_BiToBr) {
        InvolveParameter(&item, "TO_BiToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->TEMPORAL.TO_RjToBr) {
        InvolveParameter(&item, "TO_RjToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->TEMPORAL.TO_LkToBr) {
        InvolveParameter(&item, "TO_LkToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->TEMPORAL.TO_CmToBr) {
        InvolveParameter(&item, "TO_CmToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->TEMPORAL.TO_DnToBr) {
        InvolveParameter(&item, "TO_DnToBr-" + topic);
    }
    for (auto &[topic, item] : _solver->_parMagr->TEMPORAL.RS_READOUT) {
        InvolveParameter(&item, "RS_READOUT-" + topic);
    }

    // intrinsics
    for (auto &[topic, intri] : _solver->_parMagr->INTRI.IMU) {
        InvolveParameter(intri->ACCE.BIAS.data(), "ACCE.BIAS-" + topic);
        InvolveParameter(intri->GYRO.BIAS.data(), "GYRO.BIAS-" + topic);
    }

    for (auto &[topic, intri] : _solver->_parMagr->INTRI.Camera) {
        InvolveParameter(intri->FXAddress(), "FX" + topic);
        InvolveParameter(intri->FYAddress(), "FY" + topic);
        InvolveParameter(intri->CXAddress(), "CX" + topic);
        InvolveParameter(intri->CYAddress(), "CY" + topic);
    }

    for (auto &[topic, intri] : _solver->_parMagr->INTRI.RGBD) {
        InvolveParameter(intri->intri->FXAddress(), "FX" + topic);
        InvolveParameter(intri->intri->FYAddress(), "FY" + topic);
        InvolveParameter(intri->intri->CXAddress(), "CX" + topic);
        InvolveParameter(intri->intri->CYAddress(), "CY" + topic);
        InvolveParameter(&intri->alpha, "ALPHA" + topic);
        InvolveParameter(&intri->beta, "BETA" + topic);
    }

    // gravity
    InvolveParameter(_solver->_parMagr->GRAVITY.data(), "GRAVITY");

    // control points
    auto &so3Spline = _solver->_splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    auto &scaleSpline = _solver->_splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    for (int i = 0; i < std::min(static_cast<int>(so3Spline.GetKnots().size()), 15); ++i) {
        InvolveParameter(so3Spline.GetKnot(i).data(), "SO3_KNOT-" + std::to_string(i));
    }
    for (int i = 0; i < std::min(static_cast<int>(scaleSpline.GetKnots().size()), 15); ++i) {
        InvolveParameter(scaleSpline.GetKnot(i).data(), "SCALE_KNOT-" + std::to_string(i));
    }

    auto hessianMat =
        estimator->GetHessianMatrix(parAddress, Configor::Preference::AvailableThreads());
    auto filename = saveDir + "/hessian" + Configor::GetFormatExtension();
    std::ofstream file(filename);
    auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
    SerializeByOutputArchiveVariant(
        ar, Configor::Preference::OutputDataFormat, cereal::make_nvp("row", hessianMat.rows()),
        cereal::make_nvp("col", hessianMat.cols()), cereal::make_nvp("hessian", hessianMat),
        cereal::make_nvp("par_order_size", parOrderSize));
    spdlog::info("saving hessian matrix finished!");
}

void CalibSolverIO::VerifyVisualLiDARConsistency() {
    if (!Configor::IsLiDARIntegrated()) {
        return;
    }
    if (!Configor::IsCameraIntegrated() || !Configor::IsRGBDIntegrated()) {
        return;
    }

    std::string saveDir = Configor::DataStream::OutputPath + "/consistency";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving consistency images to dir: '{}'...", saveDir);
    } else {
        return;
    }

    auto covisibility = VisualLiDARCovisibility::Create(_solver->_backup->lidarMap);
    // for cameras
    std::shared_ptr<tqdm> bar;
    for (const auto &[topic, data] : _solver->_dataMagr->GetCameraMeasurements()) {
        spdlog::info("verify consistency between LiDAR and camera '{}'...", topic);

        auto subSaveDir = saveDir + '/' + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        const auto &intri = _solver->_parMagr->INTRI.Camera.at(topic);
        std::vector<std::pair<ns_veta::IndexT, Sophus::SE3d>> poseVec;
        poseVec.reserve(data.size());
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(data.size()); ++i) {
            bar->progress(i, static_cast<int>(data.size()));

            const auto &frame = data.at(i);
            auto pose = _solver->CurCmToW(frame->GetTimestamp(), topic);
            if (pose == std::nullopt) {
                continue;
            }

            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            auto filenameDepth = subSaveDir + '/' + std::to_string(frame->GetId()) + "-d.png";
            auto filenameColor = subSaveDir + '/' + std::to_string(frame->GetId()) + "-c.jpg";

            // undistorted gray image
            cv::Mat undistImgColor, res;
            undistImgColor =
                CalibParamManager::ParIntri::UndistortImage(intri, frame->GetColorImage());

            // depth image
            auto [depthImg, colorImg] = covisibility->CreateCovisibility(*pose, intri);
            cv::imwrite(filename, undistImgColor);
            cv::imwrite(filenameDepth, depthImg);
            cv::imwrite(filenameColor, colorImg);
            poseVec.emplace_back(frame->GetId(), *pose);

            // connect
            cv::hconcat(undistImgColor, colorImg, res);
            cv::imshow("Covisibility Image", res);
            cv::waitKey(1);
        }
        bar->finish();
        // save pose vector
        auto filename = subSaveDir + "/pose" + ns_ikalibr::Configor::GetFormatExtension();
        std::ofstream file(filename, std::ios::out);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(
            ar, Configor::Preference::OutputDataFormat,
            cereal::make_nvp(
                "note",
                std::string("The pose denotes transformation from camera frame to world frame!")),
            cereal::make_nvp("pose_seq", poseVec));
        spdlog::info("verify consistency for camera '{}' finished", topic);
    }
    cv::destroyAllWindows();

    // for rgbds
    for (const auto &[topic, data] : _solver->_dataMagr->GetRGBDMeasurements()) {
        spdlog::info("verify consistency between LiDAR and rgbd '{}'...", topic);

        auto subSaveDir = saveDir + '/' + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        const auto &intri = _solver->_parMagr->INTRI.RGBD.at(topic);
        std::vector<std::pair<ns_veta::IndexT, Sophus::SE3d>> poseVec;
        poseVec.reserve(data.size());
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(data.size()); ++i) {
            bar->progress(i, static_cast<int>(data.size()));

            const auto &frame = data.at(i);
            auto pose = _solver->CurDnToW(frame->GetTimestamp(), topic);
            if (pose == std::nullopt) {
                continue;
            }

            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            auto filenameDepth = subSaveDir + '/' + std::to_string(frame->GetId()) + "-d.png";
            auto filenameColor = subSaveDir + '/' + std::to_string(frame->GetId()) + "-c.jpg";
            auto filenameRaw = subSaveDir + '/' + std::to_string(frame->GetId()) + "-r.jpg";

            // undistorted gray image
            cv::Mat undistImgColor, res;
            undistImgColor =
                CalibParamManager::ParIntri::UndistortImage(intri->intri, frame->GetColorImage());

            // depth image
            auto [depthImg, colorImg] = covisibility->CreateCovisibility(*pose, intri->intri);
            cv::imwrite(filename, undistImgColor);
            cv::imwrite(filenameDepth, depthImg);
            cv::imwrite(filenameColor, colorImg);
            cv::imwrite(filenameRaw, frame->CreateColorDepthMap(intri, false));
            poseVec.emplace_back(frame->GetId(), *pose);

            // connect
            cv::hconcat(undistImgColor, colorImg, res);
            cv::imshow("Covisibility Image", res);
            cv::waitKey(1);
        }
        bar->finish();
        // save pose vector
        auto filename = subSaveDir + "/pose" + ns_ikalibr::Configor::GetFormatExtension();
        std::ofstream file(filename, std::ios::out);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(
            ar, Configor::Preference::OutputDataFormat,
            cereal::make_nvp(
                "note",
                std::string("The pose denotes transformation from rgbd frame to world frame!")),
            cereal::make_nvp("pose_seq", poseVec));
        spdlog::info("verify consistency for rgbd '{}' finished", topic);
    }
    cv::destroyAllWindows();

    spdlog::info("verify consistency finished!");
}

void CalibSolverIO::SaveVisualKinematics() {
    if (!Configor::IsCameraIntegrated() && !Configor::IsRGBDIntegrated()) {
        return;
    }

    std::string saveDir = Configor::DataStream::OutputPath + "/kinematics";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving visual kinematics to dir: '{}'...", saveDir);
    } else {
        return;
    }

    std::shared_ptr<tqdm> bar;
    // gravity
    for (const auto &[topic, data] : _solver->_dataMagr->GetCameraMeasurements()) {
        spdlog::info("create visual images with gravity for camera '{}'...", topic);

        auto subSaveDir = saveDir + "/gravity/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        auto gravityDrawer = VisualGravityDrawer::Create(
            topic, _solver->_dataMagr->GetSfMData(topic), _solver->_splines, _solver->_parMagr);
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(data.size()); ++i) {
            bar->progress(i, static_cast<int>(data.size()));
            const auto &frame = data.at(i);
            cv::Mat res = gravityDrawer->CreateGravityImg(frame);
            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            cv::imwrite(filename, res);
            cv::imshow("Visual Gravity", res);
            cv::waitKey(1);
        }
        bar->finish();
    }
    cv::destroyAllWindows();

    // gravity for rgbds
    for (const auto &[topic, data] : _solver->_backup->rgbdCorrs) {
        spdlog::info("create visual images with gravity for rgbd '{}'...", topic);

        auto subSaveDir = saveDir + "/gravity/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        auto gravityDrawer =
            RGBDVisualGravityDrawer::Create(topic, data, _solver->_splines, _solver->_parMagr);
        const auto &frames = _solver->_dataMagr->GetRGBDMeasurements(topic);
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(frames.size()); ++i) {
            bar->progress(i, static_cast<int>(frames.size()));
            const auto &frame = frames.at(i);
            cv::Mat res = gravityDrawer->CreateGravityImg(frame);
            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            cv::imwrite(filename, res);
            cv::imshow("Visual Gravity", res);
            cv::waitKey(1);
        }
        bar->finish();
    }
    cv::destroyAllWindows();

    // linear velocities
    for (const auto &[topic, data] : _solver->_dataMagr->GetCameraMeasurements()) {
        spdlog::info("create visual linear velocity images for camera '{}'...", topic);

        auto subSaveDir = saveDir + "/lin_vel/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        auto linVelDrawer = VisualLinVelDrawer::Create(topic, _solver->_dataMagr->GetSfMData(topic),
                                                       _solver->_splines, _solver->_parMagr);

        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(data.size()); ++i) {
            bar->progress(i, static_cast<int>(data.size()));
            const auto &frame = data.at(i);
            cv::Mat res = linVelDrawer->CreateLinVelImg(frame);
            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            cv::imwrite(filename, res);
            cv::imshow("Visual Linear Velocity", res);
            cv::waitKey(1);
        }
        bar->finish();
    }
    cv::destroyAllWindows();

    // linear velocities for rgbds
    for (const auto &[topic, data] : _solver->_backup->rgbdCorrs) {
        spdlog::info("create visual linear velocity images for rgbd '{}'...", topic);

        auto subSaveDir = saveDir + "/lin_vel/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        auto linVelDrawer =
            RGBDVisualLinVelDrawer::Create(topic, data, _solver->_splines, _solver->_parMagr);

        const TimeDeriv::ScaleSplineType &scaleSplineType = ns_ikalibr::CalibSolver::GetScaleType();
        const auto &frames = _solver->_dataMagr->GetRGBDMeasurements(topic);
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(frames.size()); ++i) {
            bar->progress(i, static_cast<int>(frames.size()));
            const auto &frame = frames.at(i);
            cv::Mat res = linVelDrawer->CreateLinVelImg(frame, scaleSplineType);
            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            cv::imwrite(filename, res);
            cv::imshow("Visual Linear Velocity", res);
            cv::waitKey(1);
        }
        bar->finish();
    }
    cv::destroyAllWindows();

    // angular velocities
    for (const auto &[topic, data] : _solver->_dataMagr->GetCameraMeasurements()) {
        spdlog::info("create visual angular velocity images for camera '{}'...", topic);

        auto subSaveDir = saveDir + "/ang_vel/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        auto angVelDrawer = VisualAngVelDrawer::Create(topic, _solver->_dataMagr->GetSfMData(topic),
                                                       _solver->_splines, _solver->_parMagr);

        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(data.size()); ++i) {
            bar->progress(i, static_cast<int>(data.size()));
            const auto &frame = data.at(i);
            cv::Mat res = angVelDrawer->CreateAngVelImg(frame);
            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            cv::imwrite(filename, res);
            cv::imshow("Visual Angular Velocity", res);
            cv::waitKey(1);
        }
        bar->finish();
    }

    // angular velocities for rgbds
    for (const auto &[topic, data] : _solver->_backup->rgbdCorrs) {
        spdlog::info("create visual angular velocity images for rgbd '{}'...", topic);

        auto subSaveDir = saveDir + "/ang_vel/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        auto angVelDrawer =
            RGBDVisualAngVelDrawer::Create(topic, data, _solver->_splines, _solver->_parMagr);

        const auto &frames = _solver->_dataMagr->GetRGBDMeasurements(topic);
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(frames.size()); ++i) {
            bar->progress(i, static_cast<int>(frames.size()));
            const auto &frame = frames.at(i);
            cv::Mat res = angVelDrawer->CreateAngVelImg(frame);
            auto filename = subSaveDir + '/' + std::to_string(frame->GetId()) + ".jpg";
            cv::imwrite(filename, res);
            cv::imshow("Visual Angular Velocity", res);
            cv::waitKey(1);
        }
        bar->finish();
    }

    cv::destroyAllWindows();
    spdlog::info("saving visual kinematics finished!");
}

void CalibSolverIO::SaveVisualColorizedMap() {
    if (!Configor::IsCameraIntegrated() || !Configor::IsLiDARIntegrated()) {
        return;
    }

    std::string saveDir = Configor::DataStream::OutputPath + "/maps";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving visual colorized map to dir: '{}'...", saveDir);
    } else {
        return;
    }

    for (const auto &[topic, frames] : _solver->_dataMagr->GetCameraMeasurements()) {
        spdlog::info("create colorized map using camera '{}'...", topic);

        auto subSaveDir = saveDir + "/camera/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        auto shader =
            ColorizedCloudMap::Create(topic, frames, _solver->_dataMagr->GetSfMData(topic),
                                      _solver->_splines, _solver->_parMagr);
        auto colorizedMap = shader->Colorize(_solver->_backup->lidarMap);
        auto filename = subSaveDir + "/colorized_map.pcd";
        if (pcl::io::savePCDFile(filename, *colorizedMap, true) == -1) {
            spdlog::warn("save colorized map as : '{}' failed!", filename);
        }
    }
    spdlog::info("saving visual colorized map finished!");
}

void CalibSolverIO::SaveAlignedInertialMes() {
    auto &scaleSpline = _solver->_splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    auto &so3Spline = _solver->_splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);

    // folder
    std::string saveDir = Configor::DataStream::OutputPath + "/residuals/lin_acce";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving aligned inertial measurements to dir: '{}'...", saveDir);
    } else {
        return;
    }

    // align B-spline-derived measurements to each IMU
    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        auto subSaveDir = saveDir + "/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        const IMUIntrinsics::Ptr &intri = _solver->_parMagr->INTRI.IMU.at(topic);
        const auto &SO3_BiToBr = _solver->_parMagr->EXTRI.SO3_BiToBr.at(topic);
        const Eigen::Vector3d &POS_BiInBr = _solver->_parMagr->EXTRI.POS_BiInBr.at(topic);
        const double &timeOffset = _solver->_parMagr->TEMPORAL.TO_BiToBr.at(topic);

        // inertial measurements
        std::list<IMUFrame> rawMes, estMes, diff;

        for (const auto &item : _solver->_dataMagr->GetIMUMeasurements(topic)) {
            double t = item->GetTimestamp() + timeOffset;

            if (!_solver->_splines->TimeInRange(t, scaleSpline)) {
                continue;
            }
            if (!_solver->_splines->TimeInRange(t, so3Spline)) {
                continue;
            }

            rawMes.push_back(*item);

            const auto &SO3_curBrToW = so3Spline.Evaluate(t);
            const Eigen::Vector3d &angVelInW = SO3_curBrToW * so3Spline.VelocityBody(t);
            const Eigen::Vector3d &angAcceInW = SO3_curBrToW * so3Spline.AccelerationBody(t);
            const Eigen::Matrix3d &angVelMat = Sophus::SO3d::hat(angVelInW);
            const Eigen::Matrix3d &angAcceMat = Sophus::SO3d::hat(angAcceInW);
            Eigen::Vector3d linAcce;
            switch (ns_ikalibr::CalibSolver::GetScaleType()) {
                case TimeDeriv::LIN_ACCE_SPLINE: {
                    constexpr int derive =
                        TimeDeriv::Deriv<TimeDeriv::LIN_ACCE_SPLINE, TimeDeriv::LIN_ACCE>();
                    linAcce = scaleSpline.Evaluate<derive>(t);
                } break;
                case TimeDeriv::LIN_VEL_SPLINE: {
                    constexpr int derive =
                        TimeDeriv::Deriv<TimeDeriv::LIN_VEL_SPLINE, TimeDeriv::LIN_ACCE>();
                    linAcce = scaleSpline.Evaluate<derive>(t);
                } break;
                case TimeDeriv::LIN_POS_SPLINE: {
                    constexpr int derive =
                        TimeDeriv::Deriv<TimeDeriv::LIN_POS_SPLINE, TimeDeriv::LIN_ACCE>();
                    linAcce = scaleSpline.Evaluate<derive>(t);
                } break;
            }

            const auto &est = IMUIntrinsics::KinematicsToInertialMes(
                item->GetTimestamp(),
                linAcce + (angAcceMat + angVelMat * angVelMat) * SO3_curBrToW.matrix() * POS_BiInBr,
                SO3_curBrToW * so3Spline.VelocityBody(t), SO3_curBrToW * SO3_BiToBr,
                _solver->_parMagr->GRAVITY);

            estMes.push_back(*intri->InvolveIntri(est));

            diff.emplace_back(item->GetTimestamp(),
                              rawMes.back().GetGyro() - estMes.back().GetGyro(),
                              rawMes.back().GetAcce() - estMes.back().GetAcce());
        }
        std::ofstream file(subSaveDir + "/inertial_mes" + Configor::GetFormatExtension(),
                           std::ios::out);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(
            ar, Configor::Preference::OutputDataFormat, cereal::make_nvp("raw_inertial", rawMes),
            cereal::make_nvp("est_inertial", estMes), cereal::make_nvp("inertial_diff", diff));
    }

    const IMUIntrinsics::Ptr &refIntri =
        _solver->_parMagr->INTRI.IMU.at(Configor::DataStream::ReferIMU);
    // align inertial measurements of all IMUs to the reference IMU
    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        auto subSaveDir = saveDir + "/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        const IMUIntrinsics::Ptr &intri = _solver->_parMagr->INTRI.IMU.at(topic);
        const auto &SO3_BiToBr = _solver->_parMagr->EXTRI.SO3_BiToBr.at(topic);
        const Eigen::Vector3d &POS_BiInBr = _solver->_parMagr->EXTRI.POS_BiInBr.at(topic);
        const double &timeOffset = _solver->_parMagr->TEMPORAL.TO_BiToBr.at(topic);

        std::list<IMUFrame> estMes;

        for (const auto &item : _solver->_dataMagr->GetIMUMeasurements(topic)) {
            double t = item->GetTimestamp() + timeOffset;

            if (!_solver->_splines->TimeInRange(t, scaleSpline)) {
                continue;
            }
            if (!_solver->_splines->TimeInRange(t, so3Spline)) {
                continue;
            }

            auto mesInIdeal = intri->RemoveIntri(item);

            const auto &SO3_curBrToW = so3Spline.Evaluate(t);
            const Eigen::Vector3d &angVelInW = SO3_curBrToW * so3Spline.VelocityBody(t);
            const Eigen::Vector3d &angAcceInW = SO3_curBrToW * so3Spline.AccelerationBody(t);
            const Eigen::Matrix3d &angVelMat = Sophus::SO3d::hat(angVelInW);
            const Eigen::Matrix3d &angAcceMat = Sophus::SO3d::hat(angAcceInW);

            auto res = IMUIntrinsics::InertialMesToKinematics(
                t, mesInIdeal, SO3_curBrToW * SO3_BiToBr, _solver->_parMagr->GRAVITY);

            auto alignedMes = IMUIntrinsics::KinematicsToInertialMes(
                std::get<0>(res),
                std::get<1>(res) -
                    (angAcceMat + angVelMat * angVelMat) * SO3_curBrToW.matrix() * POS_BiInBr,
                std::get<2>(res), SO3_curBrToW, _solver->_parMagr->GRAVITY);

            estMes.push_back(*refIntri->InvolveIntri(alignedMes));
        }
        std::ofstream file(subSaveDir + "/aligned_mes_to_ref" + Configor::GetFormatExtension(),
                           std::ios::out);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(ar, Configor::Preference::OutputDataFormat,
                                        cereal::make_nvp("aligned_inertial", estMes));
    }
    spdlog::info("saving aligned inertial measurements finished!");
}

void CalibSolverIO::SaveVisualReprojectionError() {
    if (!Configor::IsCameraIntegrated()) {
        return;
    }

    // folder
    std::string saveDir = Configor::DataStream::OutputPath + "/residuals/reproj_error";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving visual reprojection errors to dir: '{}'...", saveDir);
    } else {
        return;
    }

    using VisualProjFactor =
        VisualReProjFactor<Configor::Prior::SplineOrder,
                           TimeDeriv::Deriv<TimeDeriv::LIN_POS_SPLINE, TimeDeriv::LIN_POS>()>;

    for (const auto &[topic, corrsVec] : _solver->_backup->visualCorrs) {
        const double TO_CmToBr = _solver->_parMagr->TEMPORAL.TO_CmToBr.at(topic);
        const double READOUT_TIME = _solver->_parMagr->TEMPORAL.RS_READOUT.at(topic);

        const auto &intri = _solver->_parMagr->INTRI.Camera.at(topic);
        const double FX = intri->FocalX(), FX_INV = 1.0 / FX;
        const double FY = intri->FocalY(), FY_INV = 1.0 / FY;
        const double CX = intri->PrincipalPoint()(0);
        const double CY = intri->PrincipalPoint()(1);
        const double GLOBAL_SCALE = 1.0;

        const auto SE3_CmToBr = _solver->_parMagr->EXTRI.SE3_CmToBr(topic);

        auto subSaveDir = saveDir + "/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }
        std::list<Eigen::Vector2d> reprojErrors;

        for (const auto &corrs : corrsVec) {
            const double INV_DEPTH = *corrs->invDepthFir;
            const double DEPTH = 1.0 / INV_DEPTH;

            for (const auto &_corr : corrs->corrs) {
                // calculate the so3 and lin scale offset for i-feat
                double timeIByBr = _corr.ti + TO_CmToBr + _corr.li * READOUT_TIME;
                auto SE3_BrToBr0_I = _solver->CurBrToW(timeIByBr);
                if (SE3_BrToBr0_I == std::nullopt) {
                    continue;
                }

                // calculate the so3 and lin scale offset for j-feat
                auto timeJByBr = _corr.tj + TO_CmToBr + _corr.lj * READOUT_TIME;
                auto SE3_BrToBr0_J = _solver->CurBrToW(timeJByBr);
                if (SE3_BrToBr0_J == std::nullopt) {
                    continue;
                }

                Sophus::SE3d SE3_BrIToBrJ = SE3_BrToBr0_J->inverse() * *SE3_BrToBr0_I;
                Sophus::SE3d SE3_CmIToCmJ = SE3_CmToBr.inverse() * SE3_BrIToBrJ * SE3_CmToBr;

                Eigen::Vector3d PI;
                VisualProjFactor::TransformImgToCam<double>(&FX_INV, &FY_INV, &CX, &CY, _corr.fi,
                                                            &PI);
                PI *= DEPTH * GLOBAL_SCALE;

                Eigen::Vector3d PJ = SE3_CmIToCmJ * PI;
                PJ /= PJ(2);
                Eigen::Vector2d fjPred;
                VisualProjFactor::TransformCamToImg<double>(&FX, &FY, &CX, &CY, PJ, &fjPred);

                Eigen::Vector2d residuals = fjPred - _corr.fj;
                reprojErrors.push_back(residuals);
            }
        }

        std::ofstream file(subSaveDir + "/residuals" + Configor::GetFormatExtension(),
                           std::ios::out);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(ar, Configor::Preference::OutputDataFormat,
                                        cereal::make_nvp("reproj_errors", reprojErrors));
    }
    spdlog::info("saving visual reprojection errors finished!");
}

bool CalibSolverIO::SavePoseSequence(const Eigen::aligned_vector<ns_ctraj::Posed> &poseSeq,
                                     const std::string &filename,
                                     CerealArchiveType::Enum archiveType) {
    std::ofstream file(filename);
    auto ar = GetOutputArchiveVariant(file, archiveType);
    SerializeByOutputArchiveVariant(ar, archiveType, cereal::make_nvp("pose_seq", poseSeq));
    return true;
}

bool CalibSolverIO::TryCreatePath(const std::string &path) {
    if (!std::filesystem::exists(path) && !std::filesystem::create_directories(path)) {
        spdlog::warn("create directory failed: '{}'", path);
        return false;
    } else {
        return true;
    }
}

void CalibSolverIO::SaveLiDARMaps() {
    if (ns_ikalibr::CalibSolver::GetScaleType() != TimeDeriv::ScaleSplineType::LIN_POS_SPLINE ||
        !Configor::IsLiDARIntegrated()) {
        return;
    }

    std::string saveDir = Configor::DataStream::OutputPath + "/maps";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving maps to dir: '{}'...", saveDir);
    } else {
        return;
    }

    auto subSaveDir = saveDir + "/lidar";
    if (!TryCreatePath(subSaveDir)) {
        spdlog::warn("create sub directory to save lidar map failed: '{}'", subSaveDir);
    } else {
        std::map<ufo::map::Node, std::vector<PointToSurfelCorr::Ptr>> nodes;
        std::size_t count = 0;
        for (const auto &[topic, corrVec] : _solver->_backup->lidarCorrs) {
            for (const auto &corr : corrVec) {
                nodes[corr->node].push_back(corr), ++count;
            }
        }

        ColorPointCloud::Ptr surfelCloud(new ColorPointCloud);
        surfelCloud->reserve(count);
        for (const auto &[node, corrVec] : nodes) {
            auto color = ns_viewer::Entity::GetUniqueColour();
            for (const auto &corr : corrVec) {
                ColorPoint p;
                p.x = static_cast<float>(corr->pInMap(0));
                p.y = static_cast<float>(corr->pInMap(1));
                p.z = static_cast<float>(corr->pInMap(2));
                p.r = static_cast<std::uint8_t>(color.r * 255.0f);
                p.g = static_cast<std::uint8_t>(color.g * 255.0f);
                p.b = static_cast<std::uint8_t>(color.b * 255.0f);
                p.a = static_cast<std::uint8_t>(color.a * 255.0f);
                surfelCloud->push_back(p);
            }
        }
        auto filename = subSaveDir + "/surfel_map.pcd";
        spdlog::info("save global lidar surfel map...");
        if (pcl::io::savePCDFile(filename, *surfelCloud, true) == -1) {
            spdlog::warn("save surfel lidar map as : '{}' failed!", filename);
        } else {
            spdlog::info("save global lidar surfel map as '{}'", filename);
        }

        // deconstruct
        surfelCloud.reset();

        filename = subSaveDir + "/gravity_aligned_map.pcd";

        // create colorized map by aligning to the gravity
        auto colorizedMap =
            ns_viewer::AlignedCloud<IKalibrPoint>::Create(
                _solver->_backup->lidarMap, -_solver->_parMagr->GRAVITY.cast<float>(), 2.0f)
                ->GetCloud()
                .GetCloud();
        spdlog::info("save global lidar map...");
        if (pcl::io::savePCDFile(filename, *colorizedMap, true) == -1) {
            spdlog::warn("save aligned lidar map as : '{}' failed!", filename);
        } else {
            spdlog::info("save global lidar map as '{}'", filename);
        }
    }
}

void CalibSolverIO::SaveVisualMaps() {
    if (ns_ikalibr::CalibSolver::GetScaleType() != TimeDeriv::ScaleSplineType::LIN_POS_SPLINE) {
        return;
    }
    if (!Configor::IsCameraIntegrated() || !Configor::IsRGBDIntegrated()) {
        return;
    }

    std::string saveDir = Configor::DataStream::OutputPath + "/maps";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving maps to dir: '{}'...", saveDir);
    } else {
        return;
    }

    for (const auto &[topic, sfmData] : _solver->_dataMagr->GetSfMData()) {
        spdlog::info("save veta for camera '{}'...", topic);

        auto subSaveDir = saveDir + "/camera/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory to save veta for '{}' failed: '{}'", topic,
                         subSaveDir);
            continue;
        }

        auto filename = subSaveDir + "/veta.bin";
        if (!ns_veta::Save(*sfmData, filename, ns_veta::Veta::ALL)) {
            spdlog::warn("create sub directory to save veta for '{}' failed: '{}'", topic, saveDir);
        } else {
            spdlog::info("save veta for camera '{}' as '{}'", topic, filename);
        }
    }

    if (Configor::IsRGBDIntegrated()) {
        auto subSaveDir = saveDir + "/rgbd";
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory to save rgbd point cloud map failed: '{}'",
                         subSaveDir);
        } else {
            spdlog::info("save rgbd point cloud map for rgbd...");
            auto map = _solver->BuildGlobalMapOfRGBD();
            auto filename = subSaveDir + "/rgbd_map.pcd";
            if (pcl::io::savePCDFile(filename, *map, true) == -1) {
                spdlog::warn("save rgbd map as : '{}' failed!", filename);
            } else {
                spdlog::info("save rgbd map as '{}'", filename);
            }
        }
    }
}

void CalibSolverIO::SaveRadarMaps() {
    if (ns_ikalibr::CalibSolver::GetScaleType() != TimeDeriv::ScaleSplineType::LIN_POS_SPLINE ||
        !Configor::IsRadarIntegrated()) {
        return;
    }

    std::string saveDir = Configor::DataStream::OutputPath + "/maps";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving maps to dir: '{}'...", saveDir);
    } else {
        return;
    }

    auto subSaveDir = saveDir + "/radar";
    if (!TryCreatePath(subSaveDir)) {
        spdlog::warn("create sub directory to save radar target cloud map failed: '{}'",
                     subSaveDir);
    } else {
        spdlog::info("save target cloud map for radar...");
        auto radarMap = _solver->_backup->radarMap;
        auto filename = subSaveDir + "/target_map.pcd";
        if (pcl::io::savePCDFile(filename, *radarMap, true) == -1) {
            spdlog::warn("save radar map as : '{}' failed!", filename);
        } else {
            spdlog::info("save radar map as '{}'", filename);
        }
    }
}

void CalibSolverIO::SaveRadarDopplerError() {
    if (!Configor::IsRadarIntegrated()) {
        return;
    }

    // folder
    std::string saveDir = Configor::DataStream::OutputPath + "/residuals/doppler_error";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving radar radar doppler errors to dir: '{}'...", saveDir);
    } else {
        return;
    }

    const auto &so3Spline = _solver->_splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _solver->_splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    auto scaleType = ns_ikalibr::CalibSolver::GetScaleType();

    for (const auto &[topic, data] : _solver->_dataMagr->GetRadarMeasurements()) {
        auto subSaveDir = saveDir + "/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }
        std::list<double> dopplerErrors;

        const double TO_RjToBr = _solver->_parMagr->TEMPORAL.TO_RjToBr.at(topic);
        const auto SE3_RjToBr = _solver->_parMagr->EXTRI.SE3_RjToBr(topic);

        for (const auto &ary : data) {
            for (const auto &tar : ary->GetTargets()) {
                double timeByBr = tar->GetTimestamp() + TO_RjToBr;
                if (!so3Spline.TimeStampInRange(timeByBr) ||
                    !scaleSpline.TimeStampInRange(timeByBr)) {
                    continue;
                }
                // query
                Sophus::SO3d SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
                Eigen::Vector3d ANG_VEL_BrToBr0InBr = so3Spline.VelocityBody(timeByBr);
                Eigen::Vector3d ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * ANG_VEL_BrToBr0InBr;

                Eigen::Vector3d LIN_VEL_BrInBr0;
                if (scaleType == TimeDeriv::LIN_VEL_SPLINE) {
                    LIN_VEL_BrInBr0 = scaleSpline.Evaluate<0>(timeByBr);
                } else if (scaleType == TimeDeriv::LIN_POS_SPLINE) {
                    LIN_VEL_BrInBr0 = scaleSpline.Evaluate<1>(timeByBr);
                } else {
                    throw Status(Status::ERROR,
                                 "unknown scale spline type when compute the radar residuals");
                }

                Eigen::Vector3d tarInRj = tar->GetTargetXYZ();
                Eigen::Vector1d v1 = -tarInRj.transpose() * SE3_RjToBr.so3().matrix().transpose() *
                                     SO3_BrToBr0.matrix().transpose() *
                                     (-Sophus::SO3d::hat(SO3_BrToBr0 * SE3_RjToBr.translation()) *
                                          ANG_VEL_BrToBr0InBr0 +
                                      LIN_VEL_BrInBr0);

                double v2 = tar->GetRadialVelocity();
                double error = tar->GetInvRange() * v1(0) - v2;
                dopplerErrors.push_back(error);
            }
        }

        std::ofstream file(subSaveDir + "/residuals" + Configor::GetFormatExtension(),
                           std::ios::out);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(ar, Configor::Preference::OutputDataFormat,
                                        cereal::make_nvp("doppler_errors", dopplerErrors));
    }
    spdlog::info("saving radar doppler errors finished!");
}

void CalibSolverIO::SaveRGBDVelocityError() {
    if (!Configor::IsRGBDIntegrated()) {
        return;
    }

    // folder
    std::string saveDir = Configor::DataStream::OutputPath + "/residuals/rgbd_vel_error";
    if (TryCreatePath(saveDir)) {
        spdlog::info("saving radar rgbd velocity errors to dir: '{}'...", saveDir);
    } else {
        return;
    }

    const auto &so3Spline = _solver->_splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _solver->_splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    auto scaleType = ns_ikalibr::CalibSolver::GetScaleType();

    for (const auto &[topic, corrVec] : _solver->_backup->rgbdCorrs) {
        auto subSaveDir = saveDir + "/" + topic;
        if (!TryCreatePath(subSaveDir)) {
            spdlog::warn("create sub directory for '{}' failed: '{}'", topic, subSaveDir);
            continue;
        }

        const auto &parMagr = _solver->_parMagr;
        const double readout = parMagr->TEMPORAL.RS_READOUT.at(topic);
        const double TO_DnToBr = parMagr->TEMPORAL.TO_DnToBr.at(topic);
        auto SE3_DnToBr = parMagr->EXTRI.SE3_DnToBr(topic);
        const auto &intri = parMagr->INTRI.RGBD.at(topic);

        std::list<Eigen::Vector2d> velErrors;

        for (const auto &corr : corrVec) {
            double timeByBr = corr->MidPointTime(readout) + TO_DnToBr;
            if (!so3Spline.TimeStampInRange(timeByBr) || !scaleSpline.TimeStampInRange(timeByBr)) {
                continue;
            }

            auto SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
            Sophus::SO3d SO3_BrToDn = SE3_DnToBr.so3().inverse();

            Eigen::Vector3d ANG_VEL_BrToBr0InBr = so3Spline.VelocityBody(timeByBr);
            Eigen::Vector3d ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * ANG_VEL_BrToBr0InBr;
            Eigen::Vector3d ANG_VEL_DnToBr0InDn = SO3_BrToDn * ANG_VEL_BrToBr0InBr;

            Eigen::Vector3d LIN_VEL_BrToBr0InBr0;
            switch (scaleType) {
                case TimeDeriv::LIN_ACCE_SPLINE:
                    // this would not happen
                    continue;
                case TimeDeriv::LIN_VEL_SPLINE:
                    LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<0>(timeByBr);
                    break;
                case TimeDeriv::LIN_POS_SPLINE:
                    LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<1>(timeByBr);
                    break;
            }

            Eigen::Vector3d LIN_VEL_DnToBr0InBr0 =
                LIN_VEL_BrToBr0InBr0 -
                Sophus::SO3d::hat(SO3_BrToBr0 * SE3_DnToBr.translation()) * ANG_VEL_BrToBr0InBr0;
            Eigen::Vector3d LIN_VEL_DnToBr0InDn =
                SO3_BrToDn * SO3_BrToBr0.inverse() * LIN_VEL_DnToBr0InBr0;

            const double FX = intri->intri->FocalX(), FY = intri->intri->FocalY();
            const double CX = intri->intri->PrincipalPoint()(0),
                         CY = intri->intri->PrincipalPoint()(1);

            // map depth using alpha and beta
            const double depth = intri->ActualDepth(corr->depth);

            Eigen::Matrix<double, 2, 3> subAMat, subBMat;
            RGBDVelocityFactor<0, 0>::SubMats<double>(&FX, &FY, &CX, &CY, corr->MidPoint(),
                                                      &subAMat, &subBMat);
            Eigen::Vector2d pred =
                1.0 / depth * subAMat * LIN_VEL_DnToBr0InDn + subBMat * ANG_VEL_DnToBr0InDn;

            Eigen::Vector2d mes = corr->MidPointVel(readout);

            Eigen::Vector2d residuals = pred - mes;
            velErrors.push_back(residuals);
        }

        std::ofstream file(subSaveDir + "/residuals" + Configor::GetFormatExtension(),
                           std::ios::out);
        auto ar = GetOutputArchiveVariant(file, Configor::Preference::OutputDataFormat);
        SerializeByOutputArchiveVariant(ar, Configor::Preference::OutputDataFormat,
                                        cereal::make_nvp("rgbd_vel_errors", velErrors));
    }

    spdlog::info("saving rgbd velocity errors finished!");
}
}  // namespace ns_ikalibr