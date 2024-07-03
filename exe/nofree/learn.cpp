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
#include "calib/calib_param_manager.h"
#include "util/status.hpp"
#include "nofree/ufomap_learner.hpp"
#include "nofree/logo_svg.h"
#include "nofree/data_collect_demo.h"
#include "spdlog/fmt/bundled/color.h"
#include "nofree/extract_intri.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_learn");
    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // std::vector<std::string> wsVec = {
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-14-52-11",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-14-55-10",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-15-22-04",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-15-30-05",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-16-34-01",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-16-42-32",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-16-44-05",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-16-49-39",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-16-55-29",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-17-10-20",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-17-12-05",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-17-17-17",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/iKalibr-data-2024-06-25-17-33-29",
        //     "/home/csl/ros_ws/iKalibr/src/ikalibr/data/real-world/"
        //     "iKalibr-data-2024-06-25-17-37-20"};
        // ns_ikalibr::ExtractIntriParam::ToDisk(wsVec, false, true);

        ns_ikalibr::DataCollectionMotionDemo::Run(
            ns_ikalibr::DataCollectionMotionDemo::GeneratePoseSeq());
        ns_ikalibr::DataCollectionMotionDemo::Run(
            ns_ikalibr::DataCollectionMotionDemo::GeneratePoseSeqLivox());

        // ns_ikalibr::LoGoMaker::SaveToFile("/home/csl/ros_ws/iKalibr/src/ikalibr/docs/img/logo.svg");

        // ns_ikalibr::UFOMapLearner::Learn();

        // show spatiotemporal relationship in viewer
        // ns_viewer::ViewerConfigor config;
        // config.camera.initPos = {0.4f, 0.4f, 0.3f};
        // ns_viewer::Viewer viewer(config);
        // viewer.RemoveEntity();
        // ns_ikalibr::CalibParamManager::Load(
        //         "/home/csl/ros_ws/iKalibr/src/ikalibr_sim/bag/sp_param_gt.yaml"
        // )->VisualizationSensors(viewer);
        // viewer.RunInSingleThread();

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