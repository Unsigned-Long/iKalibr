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

#include "ros/ros.h"
#include "spdlog/spdlog.h"
#include "util/status.hpp"
#include "util/utils.h"
#include "pcl/filters/voxel_grid.h"
#include "ufo/map/predicate/predicates.h"
#include "tiny-viewer/core/multi_viewer.h"
#include "util/cloud_define.hpp"
#include "spdlog/fmt/bundled/color.h"
#include "util/utils_tpl.hpp"
#include "tiny-viewer/object/aligned_cloud.hpp"
#include "tiny-viewer/core/pose.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

template <class PointType>
typename pcl::PointCloud<PointType>::Ptr CloudFilter(
    const typename pcl::PointCloud<PointType>::Ptr &cloud, float minZ, float maxZ, float leafSize) {
    typename pcl::PointCloud<PointType>::Ptr downSampled(new pcl::PointCloud<PointType>);
    if (leafSize > 0.0) {
        pcl::VoxelGrid<PointType> filter;
        filter.setInputCloud(cloud);
        filter.setLeafSize(leafSize, leafSize, leafSize);
        filter.filter(*downSampled);
    } else {
        downSampled = cloud;
    }

    typename pcl::PointCloud<PointType>::Ptr newCloud(new pcl::PointCloud<PointType>);
    newCloud->reserve(downSampled->size());

    for (const auto &p : downSampled->points) {
        if (p.z < minZ || p.z > maxZ) {
            continue;
        }
        newCloud->push_back(p);
    }
    return newCloud;
}

ns_viewer::Posef GetCameraPose(const float radius, const float height, const float rotRate) {
    static float radAngle = 0.0;

    float ca = std::cos(radAngle), sa = std::sin(radAngle);
    Eigen::Vector3f pos(radius * ca, radius * sa, height);

    Eigen::Vector3f zAxis = -pos.normalized();
    Eigen::Vector3f xAxis(-sa, ca, 0.0);
    Eigen::Vector3f yAxis = zAxis.cross(xAxis);

    Eigen::Matrix3f rotMat;
    rotMat.col(0) = xAxis;
    rotMat.col(1) = yAxis;
    rotMat.col(2) = zAxis;

    radAngle += rotRate;
    if (radAngle > 2 * M_PI) {
        radAngle -= 2 * M_PI;
    }

    return ns_viewer::Posef(rotMat, pos);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_lidar_map_viewer");
    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // load parameters
        auto preDir = ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_lidar_map_viewer/pre_dir");
        spdlog::info("pre directory: '{}'", preDir);

        auto mode = ns_ikalibr::GetParamFromROS<int>("/ikalibr_lidar_map_viewer/mode");
        if (mode != 0 && mode != 1 && mode != 2) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR, "wrong mode is set!");
        }

        auto outputDir =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_lidar_map_viewer/output_dir");

        auto leafSize = (float)ns_ikalibr::GetParamFromROS<double>(
            "/ikalibr_lidar_map_viewer/down_sample_leaf_size");
        auto minZ = (float)ns_ikalibr::GetParamFromROS<double>("/ikalibr_lidar_map_viewer/min_z");
        auto maxZ = (float)ns_ikalibr::GetParamFromROS<double>("/ikalibr_lidar_map_viewer/max_z");
        spdlog::info("min z value: '{:.3f}', max z value: '{:.3f}', leaf size: '{:.3f}'", minZ,
                     maxZ, leafSize);

        auto camRadius =
            (float)ns_ikalibr::GetParamFromROS<double>("/ikalibr_lidar_map_viewer/cam_radius");
        auto camHeight =
            (float)ns_ikalibr::GetParamFromROS<double>("/ikalibr_lidar_map_viewer/cam_height");
        auto rotRate =
            (float)ns_ikalibr::GetParamFromROS<double>("/ikalibr_lidar_map_viewer/rot_rate");
        spdlog::info("camera radius: '{:.3f}', camera height: '{:.3f}', rot rate: '{:.3f}'",
                     camRadius, camHeight, rotRate);

        auto winScaleStr =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_lidar_map_viewer/win_scale");
        auto wh = ns_ikalibr::SplitString(winScaleStr, ':');
        if (wh.size() != 2) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR, "wrong scale of window is set!");
        }
        double scale = std::stod(wh[0]) / std::stod(wh[1]);
        spdlog::info("scale of window: '{:.3f}'", scale);

        auto winGridSize =
            ns_ikalibr::GetParamFromROS<int>("/ikalibr_lidar_map_viewer/win_grid_size");
        if (winGridSize < 10) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "wrong size of window grid is set!");
        }
        spdlog::info("size of window grid: '{}'", winGridSize);

        auto wsDirVecSrc = ns_ikalibr::GetParamFromROS<std::vector<std::string>>(
            "/ikalibr_lidar_map_viewer/ws_dir_vec");
        std::vector<std::pair<std::string, int>> wsDirVec;
        for (const auto &dir : wsDirVecSrc) {
            auto ws = preDir + '/' + dir;
            auto alignedMap = ws + "/maps/lidar_aligned_map.pcd";
            auto surfelMap = ws + "/maps/lidar_surfel_map.pcd";
            if (mode == 0) {
                wsDirVec.emplace_back(alignedMap, 0);
            } else if (mode == 1) {
                wsDirVec.emplace_back(surfelMap, 1);
            } else {
                // mode == 2
                wsDirVec.emplace_back(alignedMap, 0);
                wsDirVec.emplace_back(surfelMap, 1);
            }
        }

        int row = std::ceil(std::sqrt(static_cast<double>(wsDirVec.size()) / scale));
        int col = std::ceil(static_cast<double>(wsDirVec.size()) / row);
        spdlog::info("rows: '{}', cols: '{}'", row, col);

        // create viewer names
        std::vector<std::string> viewerNames(row * col);
        for (int i = 0; i < static_cast<int>(viewerNames.size()); ++i) {
            if (i < static_cast<int>(wsDirVec.size())) {
                viewerNames.at(i) = wsDirVec.at(i).first;
            } else {
                viewerNames.at(i) = std::to_string(i);
            }
        }

        ns_viewer::MultiViewerConfigor configor(viewerNames, "lidar-map-viewer");
        configor.window.width = col * winGridSize;
        configor.window.height = row * winGridSize;
        configor.output.dataOutputPath = outputDir;

        for (const auto &name : viewerNames) {
            auto &cam = configor.camera.at(name);
            cam.height = winGridSize;
            cam.width = winGridSize;
            cam.cx = winGridSize * 0.5;
            cam.cy = winGridSize * 0.5;

            configor.grid.at(name).showGrid = false;
        }

        ns_viewer::MultiViewer viewer(configor);
        viewer.RunInMultiThread();
        // set init camera poses
        for (const auto &name : viewerNames) {
            viewer.SetCamView(GetCameraPose(camRadius, camHeight, rotRate), name);
        }

        PosPointCloud::Ptr alignedCloudCopy = nullptr;
        for (const auto &[filename, type] : wsDirVec) {
            spdlog::info("load pcd from '{}'...", filename);

            if (type == 0) {
                PosPointCloud::Ptr alignedCloud(new PosPointCloud);
                if (pcl::io::loadPCDFile(filename, *alignedCloud) == -1) {
                    spdlog::warn("load lidar aligned map from '{}' failed!", filename);
                    continue;
                }
                alignedCloud = CloudFilter<PosPoint>(alignedCloud, minZ, maxZ, leafSize);
                viewer.AddEntity(
                    ns_viewer::AlignedCloud<PosPoint>::Create(alignedCloud, {0, 0, -1}), filename);
                // aligned cloud is stored for surfel map visualization
                if (mode == 2) {
                    alignedCloudCopy = alignedCloud;
                }
            } else {
                // type == 1

                // load surfel map
                ColorPointCloud::Ptr surfelCloud(new ColorPointCloud);
                if (pcl::io::loadPCDFile(filename, *surfelCloud) == -1) {
                    spdlog::warn("load lidar surfel map from '{}' failed!", filename);
                    continue;
                }

                surfelCloud = CloudFilter<ColorPoint>(surfelCloud, minZ, maxZ, leafSize * 0.5f);
                viewer.AddEntity(ns_viewer::Cloud<ColorPoint>::Create(surfelCloud), filename);
                surfelCloud.reset();

                // load aligned map
                PosPointCloud::Ptr alignedCloud(new PosPointCloud);

                // if aligned is loaded, use it
                if (alignedCloudCopy != nullptr) {
                    alignedCloud = alignedCloudCopy;
                    alignedCloudCopy = nullptr;
                } else {
                    std::string alignedFilename = filename, search = "lidar_surfel_map",
                                replace = "lidar_aligned_map";
                    size_t pos = alignedFilename.find(search);
                    if (pos != std::string::npos) {
                        alignedFilename.replace(pos, search.length(), replace);
                    }

                    if (pcl::io::loadPCDFile(alignedFilename, *alignedCloud) == -1) {
                        spdlog::warn("load aligned lidar map from '{}' failed!", alignedFilename);
                        continue;
                    }
                    alignedCloud = CloudFilter<PosPoint>(alignedCloud, minZ, maxZ, leafSize);
                }

                pcl::PointCloud<pcl::PointXYZ>::Ptr alignedCloudRaw(
                    new pcl::PointCloud<pcl::PointXYZ>);
                pcl::copyPointCloud(*alignedCloud, *alignedCloudRaw);
                alignedCloud.reset();

                auto color = ns_viewer::Colour::Black().WithAlpha(0.2f);
                viewer.AddEntity(ns_viewer::Cloud<pcl::PointXYZ>::Create(alignedCloudRaw,
                                                                         DefaultPointSize, color),
                                 filename);
                alignedCloudRaw.reset();
            }
        }
        alignedCloudCopy.reset();

        if (rotRate > 0.0) {
            ros::start();
            ros::Rate r(25);
            while (ros::ok() && viewer.IsActive()) {
                for (const auto &[filename, type] : wsDirVec) {
                    viewer.SetCamView(GetCameraPose(camRadius, camHeight, rotRate), filename);
                }
                r.sleep();
            }
        }

    } catch (const ns_ikalibr::IKalibrStatus &status) {
        // if error happened, print it
        static const auto FStyle = fmt::emphasis::italic | fmt::fg(fmt::color::green);
        static const auto WECStyle = fmt::emphasis::italic | fmt::fg(fmt::color::red);
        switch (status.flag) {
            case ns_ikalibr::Status::FINE:
                // this case usually won't happen
                spdlog::info(fmt::format(FStyle, "{}", status.what));
                break;
            case ns_ikalibr::Status::WARNING:
                spdlog::warn(fmt::format(WECStyle, "{}", status.what));
                break;
            case ns_ikalibr::Status::ERROR:
                spdlog::error(fmt::format(WECStyle, "{}", status.what));
                break;
            case ns_ikalibr::Status::CRITICAL:
                spdlog::critical(fmt::format(WECStyle, "{}", status.what));
                break;
        }
    } catch (const std::exception &e) {
        // an unknown exception not thrown by this program
        static const auto WECStyle = fmt::emphasis::italic | fmt::fg(fmt::color::red);
        spdlog::critical(fmt::format(WECStyle, "unknown error happened: '{}'", e.what()));
    }

    ros::shutdown();
    return 0;
}