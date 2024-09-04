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
#include "filesystem"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_bag_merge");

    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // load settings
        auto configPath =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_bag_merge/config_path");
        spdlog::info("loading configure from json file '{}'...", configPath);

        auto mergeConfigor = ns_ikalibr::MergeConfigor::LoadConfigure(
            configPath, ns_ikalibr::CerealArchiveType::Enum::YAML);
        if (!mergeConfigor) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::CRITICAL, "load configure file failed!");
        } else {
            mergeConfigor->PrintMainFields();
        }

        if (std::filesystem::exists(mergeConfigor->_outputBagPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::CRITICAL,
                                     "the output bag exists! delete it then rerun as I do not "
                                     "known whether it's valuable!");
        }

        auto [begTime, endTime] = ns_ikalibr::BagMerger::Create(mergeConfigor)->Process();
        if (mergeConfigor->_bagBeginTime < 0.0 && mergeConfigor->_bagEndTime < 0.0) {
            spdlog::info("do not perform rosbag filter, as configured.");
        } else {
            spdlog::info("perform rosbag filter, as configured.");
            spdlog::info("origin time span: '{:.5f}' to '{:.5f}'", begTime.toSec(),
                         endTime.toSec());

            if (mergeConfigor->_bagEndTime > 0.0) {
                auto tarEndTime = begTime + ros::Duration(mergeConfigor->_bagEndTime);
                if (tarEndTime > endTime) {
                    spdlog::warn(
                        "desired end time '{:.5f}' is out of the bag's data range, set end time to "
                        "'{:.5f}'.",
                        tarEndTime.toSec(), endTime.toSec());
                } else {
                    endTime = tarEndTime;
                }
            }
            if (mergeConfigor->_bagBeginTime > 0.0) {
                auto tarBegTime = begTime + ros::Duration(mergeConfigor->_bagBeginTime);
                if (tarBegTime >= endTime) {
                    spdlog::warn(
                        "desired begin time '{:.5f}' is out of the bag's data range, set begin "
                        "time to '{:.5f}'.",
                        tarBegTime.toSec(), begTime.toSec());
                } else {
                    begTime = tarBegTime;
                }
            }

            spdlog::info("configured time span: '{:.5f}' to '{:.5f}'", begTime.toSec(),
                         endTime.toSec());
            auto newPath = std::filesystem::path(mergeConfigor->_outputBagPath);
            newPath.replace_filename("filtered-" + newPath.filename().string());

            auto res = std::system(fmt::format("rosbag filter {} {} "
                                               "\"t.to_sec() > {:.5f} and t.to_sec() < {:.5f}\"",
                                               mergeConfigor->_outputBagPath, newPath.string(),
                                               begTime.toSec(), endTime.toSec())
                                       .c_str());
            if (res != 0) {
                spdlog::warn("filter rosbag failed!!!");
            } else {
                spdlog::warn("filter rosbag as '{}'", newPath.string());
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