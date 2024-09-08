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

#include "nofree/data_collect_demo.h"
#include "tiny-viewer/core/multi_viewer.h"
#include "tiny-viewer/object/imu.h"
#include "tiny-viewer/entity/line.h"
#include "ros/package.h"
#include "filesystem"
#include "spdlog/spdlog.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void DataCollectionMotionDemo::Run(const std::array<std::vector<ns_viewer::Posef>, 2> &poseSeq) {
    auto [poses, views] = poseSeq;
    auto config =
        ns_viewer::MultiViewerConfigor({"GLOBAL", "FOLLOW"}, "iKalibr Data Collection Demo");
    config.grid.at("GLOBAL").showGrid = false;
    config.grid.at("FOLLOW").showGrid = false;
    config.camera.at("GLOBAL").initPos = {4.0f, 4.0f, 2.0f};
    ns_viewer::MultiViewer viewer(config);
    viewer.RunInMultiThread();
    viewer.WaitForActive(-1.0);
    auto modelPath = ros::package::getPath("ikalibr") + "/model/ikalibr.obj";
    if (std::filesystem::exists(modelPath)) {
        viewer.AddObjEntity(modelPath, "GLOBAL");
        viewer.AddObjEntity(modelPath, "FOLLOW");
    } else {
        spdlog::warn("can not load models from '{}'!", modelPath);
    }

    std::vector<ns_viewer::Entity::Ptr> poseEntities(poses.size() - 1);
    for (int i = 0; i < static_cast<int>(poses.size()) - 1; ++i) {
        poseEntities.at(i) = ns_viewer::Line::Create(
            poses.at(i).translation, poses.at(i + 1).translation, ns_viewer::Colour::Black());
    }
    viewer.AddEntity(poseEntities, "GLOBAL");
    viewer.AddEntity(poseEntities, "FOLLOW");

    auto lastTime = std::chrono::steady_clock::now();
    int curIdx = 0;
    std::map<std::string, std::vector<std::size_t>> lastEntities;

    while (viewer.IsActive()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(int(Dt * 1E3 * 1E-1)));
        std::chrono::duration<double, std::milli> curWaitTime =
            std::chrono::steady_clock::now() - lastTime;
        if (curWaitTime.count() > Dt * 1E3) {
            // reset time
            lastTime = std::chrono::steady_clock::now();

            // update viewer
            viewer.RemoveEntity(lastEntities["GLOBAL"], "GLOBAL");
            viewer.RemoveEntity(lastEntities["FOLLOW"], "FOLLOW");

            const int minIdx = std::max(curIdx - KeepNum, 0);
            const int maxIdx = curIdx, num = maxIdx - minIdx + 1;
            std::vector<ns_viewer::Entity::Ptr> entities(num);
            for (int i = minIdx; i <= curIdx; ++i) {
                entities.at(i - minIdx) = ns_viewer::Coordinate::Create(
                    poses.at(i), CoordSizeMin + (float)(i - minIdx) / (float)num * CoordSizeRange);
            }
            entities.push_back(ns_viewer::IMU::Create(poses.at(curIdx)));
            lastEntities["GLOBAL"] = viewer.AddEntity(entities, "GLOBAL");

            // view
            lastEntities["FOLLOW"] = {
                viewer.AddEntity(ns_viewer::IMU::Create(poses.at(curIdx)), "FOLLOW")};
            viewer.SetCamView(views.at(curIdx), "FOLLOW");

            ++curIdx;
            if (curIdx > int(poses.size() - 1)) {
                curIdx = 0;
            }
        }
    }
}

std::array<std::vector<ns_viewer::Posef>, 2> DataCollectionMotionDemo::GeneratePoseSeq() {
    std::vector<ns_viewer::Posef> poses(PoseSize), views(PoseSize);
    for (int i = 0; i < PoseSize; ++i) {
        const float cTheta = std::cos(DTheta * (float)i);
        const float sTheta = std::sin(DTheta * (float)i);

        Eigen::Vector3f CircleInOrigin(Radius * cTheta, Radius * sTheta, 0.0f);
        Eigen::Matrix3f CircleToOrigin;
        CircleToOrigin.col(0) = Eigen::Vector3f(cTheta, sTheta, 0.0f);
        CircleToOrigin.col(1) = Eigen::Vector3f(-sTheta, cTheta, 0.0f);
        CircleToOrigin.col(2) = Eigen::Vector3f(0.0f, 0.0f, 1.0f);

        Eigen::Vector3f EightInCircle;
        const float alpha = std::fmod((float)i * Dt, EightTime) / EightTime * (2.0 * M_PI);
        EightInCircle(0) = EightWidth * std::cos(alpha);
        EightInCircle(1) = 0.0f;
        EightInCircle(2) = EightHeight * std::sin(alpha) * std::cos(alpha);

        auto &pose = poses.at(i);
        pose.timeStamp = Dt * i;
        pose.translation = CircleToOrigin * EightInCircle + CircleInOrigin;
        pose.rotation.col(1) = CircleToOrigin.col(1);
        Eigen::Vector3f midHeader =
            CircleInOrigin + Eigen::Vector3f(0.0f, 0.0f, 5.0f * EightHeight);
        pose.rotation.col(2) = (midHeader - pose.translation).normalized();
        pose.rotation.col(0) = pose.rotation.col(1).cross(pose.rotation.col(2));

        // views
        auto &view = views.at(i);
        view.timeStamp = Dt * i;

        view.rotation.col(0) = pose.rotation.col(0);
        view.rotation.col(1) = -pose.rotation.col(2);
        view.rotation.col(2) = pose.rotation.col(1);
        view.translation = pose.translation - 0.2f * view.rotation.col(2);
    }
    return {poses, views};
}

std::array<std::vector<ns_viewer::Posef>, 2> DataCollectionMotionDemo::GeneratePoseSeqLivox() {
    std::vector<ns_viewer::Posef> poses(PoseSize), views(PoseSize);
    for (int i = 0; i < PoseSize; ++i) {
        auto &pose = poses.at(i);
        pose.timeStamp = Dt * i;

        const float sTheta = std::sin((M_PI * 6 / (PoseSize * 0.17f)) * (float)i);
        const float cTheta = std::cos((M_PI * 6 / (PoseSize * 0.17f)) * (float)i);

        const float theta = sTheta * M_PI / 6;

        if (i < PoseSize * 0.17f * 1.0f) {
            pose.translation = Eigen::Vector3f ::Zero();
            pose.rotation.col(0) = Eigen::Vector3f::UnitX();
            pose.rotation.col(1) =
                Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitX()) * Eigen::Vector3f::UnitY();
            pose.rotation.col(2) = pose.rotation.col(0).cross(pose.rotation.col(1));
        } else if (i < PoseSize * 0.17f * 2.0f) {
            pose.translation = Eigen::Vector3f ::Zero();
            pose.rotation.col(1) = Eigen::Vector3f::UnitY();
            pose.rotation.col(0) =
                Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitY()) * Eigen::Vector3f::UnitX();
            pose.rotation.col(2) = pose.rotation.col(0).cross(pose.rotation.col(1));
        } else if (i < PoseSize * 0.17f * 3.0f) {
            pose.translation = Eigen::Vector3f ::Zero();
            pose.rotation.col(2) = Eigen::Vector3f::UnitZ();
            pose.rotation.col(0) =
                Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitZ()) * Eigen::Vector3f::UnitX();
            pose.rotation.col(1) = pose.rotation.col(2).cross(pose.rotation.col(0));
        } else if (i < PoseSize * 0.17f * 4.0f) {
            pose.rotation = Eigen::Matrix3f ::Identity();
            pose.translation = Eigen::Vector3f(0.0f, cTheta, sTheta);
        } else if (i < PoseSize * 0.17f * 5.0f) {
            pose.rotation = Eigen::Matrix3f ::Identity();
            pose.translation = Eigen::Vector3f(cTheta, 0.0f, sTheta);
        } else {
            pose.rotation = Eigen::Matrix3f ::Identity();
            pose.translation = Eigen::Vector3f(cTheta, sTheta, 0.0f);
        }

        // views
        auto &view = views.at(i);
        view.timeStamp = Dt * i;

        view.rotation.col(0) = pose.rotation.col(0);
        view.rotation.col(1) = -pose.rotation.col(2);
        view.rotation.col(2) = pose.rotation.col(1);
        view.translation = pose.translation - 0.2f * view.rotation.col(2);
    }
    return {poses, views};
}
}  // namespace ns_ikalibr
