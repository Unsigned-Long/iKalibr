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

#include "nofree/bag_merge.h"
#include "util/cereal_archive_helper.hpp"
#include "spdlog/fmt/bundled/color.h"
#include "spdlog/spdlog.h"
#include "util/status.hpp"
#include "util/utils_tpl.hpp"
#include "filesystem"
#include "rosbag/bag.h"
#include "rosbag/view.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_bag_topic_downsample");

    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // load parameters
        auto iBagPath = ns_ikalibr::GetParamFromROS<std::string>(
            "/ikalibr_bag_topic_downsample/input_bag_path");
        if (!std::filesystem::exists(iBagPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR, "the bag path not exists!!! '{}'",
                                     iBagPath);
        } else {
            spdlog::info("the path of rosbag: '{}'", iBagPath);
        }

        auto topic = ns_ikalibr::GetParamFromROS<std::string>(
            "/ikalibr_bag_topic_downsample/topic_to_downsample");
        spdlog::info("ros topic to down sampled: '{}'", topic);

        auto oBagPath = ns_ikalibr::GetParamFromROS<std::string>(
            "/ikalibr_bag_topic_downsample/output_bag_path");
        if (std::filesystem::exists(oBagPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::WARNING,
                                     "the output bag exists: '{}', please delete it first!",
                                     oBagPath);
        } else {
            spdlog::info("the path of output rosbag: '{}'", oBagPath);
        }

        auto frequency =
            ns_ikalibr::GetParamFromROS<double>("/ikalibr_bag_topic_downsample/desired_frequency");
        if (frequency < 1E-3) {
            spdlog::error("invalid desired frequency: '{:.3f}'", frequency);
        } else {
            spdlog::info("the desired frequency of topic '{}': '{:.3f}'", topic, frequency);
        }

        // open rosbag
        auto srcBag = std::make_unique<rosbag::Bag>();
        srcBag->open(iBagPath, rosbag::BagMode::Read);
        auto view = rosbag::View();
        view.addQuery(*srcBag, rosbag::TopicQuery(topic));
        spdlog::info("load message of topic '{}' and sort them by timestamps...", topic);
        std::list<rosbag::View::iterator> msgs;
        for (auto iter = view.begin(); iter != view.end(); ++iter) {
            msgs.push_back(iter);
        }
        struct {
            bool operator()(const rosbag::View::iterator &i1,
                            const rosbag::View::iterator &i2) const {
                return i1->getTime().toSec() < i2->getTime().toSec();
            }
        } customLess;
        // sort messages according to the time stamps
        msgs.sort(customLess);

        auto dstBag = std::make_unique<rosbag::Bag>();
        dstBag->open(oBagPath, rosbag::BagMode::Write);
        double lastTimeSed = msgs.front()->getTime().toSec();
        const double DeltaTime = 1.0 / frequency;
        spdlog::info("write message of topic '{}'...", topic);
        for (const auto &iter : msgs) {
            if (iter->getTime().toSec() - lastTimeSed > DeltaTime) {
                dstBag->write(topic, iter->getTime(), *iter, iter->getConnectionHeader());
                lastTimeSed = iter->getTime().toSec();
            }
        }
        spdlog::info("write message of topic '{}' finished!", topic);
        dstBag->close();
        srcBag->close();

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