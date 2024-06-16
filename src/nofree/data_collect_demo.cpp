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

_3_

namespace ns_ikalibr {
    void DataCollectionMotionDemo::Run() {
        auto [poses, cams] = GeneratePoseSeq();
        ns_viewer::MultiViewer viewer(
                ns_viewer::MultiViewerConfigor({"GLOBAL", "FOLLOW"}, "iKalibr Data Collection Demo")
        );
        viewer.RunInMultiThread();
        viewer.WaitForActive(-1.0);

        std::vector<ns_viewer::Entity::Ptr> poseEntities(poses.size() - 1);
        for (int i = 0; i < static_cast<int>(poses.size()) - 1; ++i) {
            poseEntities.at(i) = ns_viewer::Line::Create(
                    poses.at(i).translation, poses.at(i + 1).translation, ns_viewer::Colour::Black()
            );
        }
        viewer.AddEntity(poseEntities, "GLOBAL");
        viewer.AddEntity(poseEntities, "FOLLOW");

        auto lastTime = std::chrono::steady_clock::now();
        int curIdx = 0;
        std::map<std::string, std::vector<std::size_t>> lastEntities;

        while (viewer.IsActive()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(int(Dt * 1E3 * 1E-1)));
            std::chrono::duration<double, std::milli> curWaitTime = std::chrono::steady_clock::now() - lastTime;
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
                            poses.at(i), CoordSizeMin + (float) (i - minIdx) / (float) num * CoordSizeRange
                    );
                }
                lastEntities["FOLLOW"] = viewer.AddEntity(entities, "FOLLOW");
                entities.push_back(ns_viewer::CubeCamera::Create(cams.at(curIdx)));
                lastEntities["GLOBAL"] = viewer.AddEntity(entities, "GLOBAL");

                // camera
                viewer.SetCamView(cams.at(curIdx), "FOLLOW");

                ++curIdx;
                if (curIdx > int(poses.size() - 1)) { curIdx = 0; }
            }
        }
    }

    std::array<std::vector<ns_viewer::Posef>, 2> DataCollectionMotionDemo::GeneratePoseSeq() {
        std::vector<ns_viewer::Posef> poses(PoseSize), cams(PoseSize);
        for (int i = 0; i < PoseSize; ++i) {
            const float cTheta = std::cos(DTheta * (float) i);
            const float sTheta = std::sin(DTheta * (float) i);

            Eigen::Vector3f CircleInOrigin(Radius * cTheta, Radius * sTheta, 0.0f);
            Eigen::Matrix3f CircleToOrigin;
            CircleToOrigin.col(0) = Eigen::Vector3f(cTheta, sTheta, 0.0f);
            CircleToOrigin.col(1) = Eigen::Vector3f(-sTheta, cTheta, 0.0f);
            CircleToOrigin.col(2) = Eigen::Vector3f(0.0f, 0.0f, 1.0f);

            Eigen::Vector3f EightInCircle;
            const float alpha = std::fmod((float) i * Dt, EightTime) / EightTime * (2.0 * M_PI);
            EightInCircle(0) = EightWidth * std::cos(alpha);
            EightInCircle(1) = 0.0f;
            EightInCircle(2) = EightHeight * std::sin(alpha) * std::cos(alpha);

            auto &pose = poses.at(i);
            pose.timeStamp = Dt * i;
            pose.translation = CircleToOrigin * EightInCircle + CircleInOrigin;
            pose.rotation.col(0) = pose.translation.normalized();
            pose.rotation.col(1) = Eigen::Vector3f(
                    -pose.translation(1), pose.translation(0), 0.0f
            ).normalized();
            pose.rotation.col(2) = pose.rotation.col(0).cross(pose.rotation.col(1));

            // cameras
            auto &cam = cams.at(i);
            cam.timeStamp = Dt * i;
            cam.translation = CircleInOrigin
                              // z-axis bias
                              + Eigen::Vector3f(0.0, 0.0, CameraZBias)
                              // x-axis and y-axis bias
                              + Eigen::Vector3f(sTheta, -cTheta, 0.0f) * CameraYBias;

            cam.rotation.col(2) = (pose.translation - cam.translation).normalized();
            cam.rotation.col(0) = Eigen::Vector3f(
                    cam.rotation.col(2)(1), -cam.rotation.col(2)(0), 0.0f
            ).normalized();
            cam.rotation.col(1) = cam.rotation.col(2).cross(cam.rotation.col(0));
        }
        return {poses, cams};
    }
}

