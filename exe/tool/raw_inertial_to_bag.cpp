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
#include "util/utils.h"
#include "util/utils_tpl.hpp"
#include "util/status.hpp"
#include "cereal/types/utility.hpp"
#include "rosbag/bag.h"
#include "sensor_msgs/Imu.h"
#include "filesystem"
#include "spdlog/fmt/bundled/color.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_raw_inertial_to_bag");
    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // load parameters
        auto rawInertialPath = ns_ikalibr::GetParamFromROS<std::string>(
            "/ikalibr_raw_inertial_to_bag/raw_inertial_path");
        spdlog::info("the path of raw inertial measurements: '{}'", rawInertialPath);
        if (!std::filesystem::exists(rawInertialPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "the raw inertial path not exists!!! '{}'", rawInertialPath);
        }

        auto imuTopic =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_raw_inertial_to_bag/imu_topic");
        spdlog::info("the topic of imu: '{}'", imuTopic);

        auto bagPath =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_raw_inertial_to_bag/bag_path");
        spdlog::info("the path of rosbag to output: '{}'", bagPath);
        auto parPath = std::filesystem::path(bagPath).parent_path();
        if (!std::filesystem::exists(parPath) && !std::filesystem::create_directories(parPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "the folder to output rosbag not exists!!! '{}'",
                                     parPath.string());
        }

        auto scaleFactorsStr =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_raw_inertial_to_bag/scale_factors");
        auto scaleFactorsVec = ns_ikalibr::SplitString(scaleFactorsStr, ';');
        if (scaleFactorsVec.size() != 3) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "the element count in 'scale_factors' should be three!!!");
        }
        double stampToSedScale = std::stod(scaleFactorsVec.at(0));
        spdlog::info("raw stamp to second-unit timestamp scale: '{}'", stampToSedScale);
        if (stampToSedScale <= 0.0) {
            throw ns_ikalibr::Status(
                ns_ikalibr::Status::ERROR,
                "raw stamp to second-unit timestamp scale should be positive!!! '{}'");
        }
        double gyroScale = std::stod(scaleFactorsVec.at(1));
        spdlog::info("gyroscope measurements scale: '{}'", gyroScale);
        double acceScale = std::stod(scaleFactorsVec.at(2));
        spdlog::info("accelerator measurements scale: '{}'", acceScale);

        auto headLineCount =
            ns_ikalibr::GetParamFromROS<int>("/ikalibr_raw_inertial_to_bag/head_line_count");
        spdlog::info("head line count: '{}'", headLineCount);

        auto splitCharacter = ns_ikalibr::GetParamFromROS<std::string>(
            "/ikalibr_raw_inertial_to_bag/split_character");
        if (splitCharacter.size() != 3) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "the split character should be only one!!!");
        }
        char splitor = splitCharacter.at(1);

        auto inertialOrder =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_raw_inertial_to_bag/inertial_order");
        auto orderStrVec = ns_ikalibr::SplitString(inertialOrder, ';');
        if (orderStrVec.size() != 7) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "the element count in 'inertial_order' should be seven!!!");
        }
        spdlog::info("inertial order: '{}'", inertialOrder);
        constexpr static const char *TIME = "t";
        constexpr static const char *AX = "ax";
        constexpr static const char *AY = "ay";
        constexpr static const char *AZ = "az";
        constexpr static const char *GX = "gx";
        constexpr static const char *GY = "gy";
        constexpr static const char *GZ = "gz";
        std::set<std::string> itemSet = {TIME, AX, AY, AZ, GX, GY, GZ};
        for (const auto &item : orderStrVec) {
            if (itemSet.count(item) == 0) {
                throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR, "the item '{}' is invalid!!!",
                                         item);
            }
        }

        std::map<std::string, int> index;
        for (int idx = 0; idx < static_cast<int>(orderStrVec.size()); ++idx) {
            const auto &item = orderStrVec.at(idx);
            index.insert({item, idx});
            std::cout << item << ": " << idx << " ";
        }
        std::cout << std::endl;
        if (index.size() != 7) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "there are same items in 'inertial_order'seven!!!");
        }

        auto dstBag = std::make_unique<rosbag::Bag>();
        dstBag->open(bagPath, rosbag::BagMode::Write);

        std::ifstream file(rawInertialPath, std::ios::in);
        std::string line;
        for (int i = 0; i < headLineCount; ++i) {
            std::getline(file, line);
        }
        while (std::getline(file, line)) {
            auto elems = ns_ikalibr::SplitString(line, splitor, true);
            if (elems.size() != 7) {
                spdlog::warn("wrong decoded line: {}", line);
                continue;
            }
            sensor_msgs::Imu imu;
            imu.header.stamp = ros::Time(std::stod(elems.at(index.at(TIME))) * stampToSedScale);
            imu.linear_acceleration.x = std::stod(elems.at(index.at(AX))) * acceScale;
            imu.linear_acceleration.y = std::stod(elems.at(index.at(AY))) * acceScale;
            imu.linear_acceleration.z = std::stod(elems.at(index.at(AZ))) * acceScale;
            imu.angular_velocity.x = std::stod(elems.at(index.at(GX))) * gyroScale;
            imu.angular_velocity.y = std::stod(elems.at(index.at(GY))) * gyroScale;
            imu.angular_velocity.z = std::stod(elems.at(index.at(GZ))) * gyroScale;
            // std::cout << imu.header.stamp << std::endl;
            // std::cout << imu.linear_acceleration << std::endl;
            // std::cout << imu.angular_velocity << std::endl;
            dstBag->write(imuTopic, imu.header.stamp, imu);
            // std::cin.get();
        }
        file.close();
        dstBag->close();
        spdlog::info("raw inertial in '{}' have been writen to rosbag as '{}'!", rawInertialPath,
                     bagPath);

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