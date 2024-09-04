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

#include "util/utils.h"
#include "spdlog/spdlog.h"
#include "util/status.hpp"
#include "util/utils_tpl.hpp"
#include "cereal/types/utility.hpp"
#include "opencv2/imgcodecs.hpp"
#include "sensor_msgs/Image.h"
#include "cv_bridge/cv_bridge.h"
#include "rosbag/bag.h"
#include "filesystem"
#include "spdlog/fmt/bundled/color.h"
#include "regex"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_imgs_to_bag");
    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // load parameters
        auto imgPath = ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_imgs_to_bag/imgs_path");
        spdlog::info("the path of images: '{}'", imgPath);
        if (!std::filesystem::exists(imgPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "the raw image path not exists!!! '{}'", imgPath);
        }
        auto nameStampRegexStr =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_imgs_to_bag/valid_name_stamp_regex");
        spdlog::info("the regex pattern of image names and stamps: '{}'", nameStampRegexStr);

        auto imgsTopic =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_imgs_to_bag/imgs_topic");
        spdlog::info("the topic of images: '{}'", imgsTopic);

        auto bagPath = ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_imgs_to_bag/bag_path");
        spdlog::info("the path of rosbag to output: '{}'", bagPath);
        auto parPath = std::filesystem::path(bagPath).parent_path();
        if (!std::filesystem::exists(parPath) && !std::filesystem::create_directories(parPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::ERROR,
                                     "the folder to output rosbag not exists!!! '{}'",
                                     parPath.string());
        }

        auto encoding = ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_imgs_to_bag/encoding");
        spdlog::info("the encoding of images: '{}'", encoding);

        auto downsampleNum =
            ns_ikalibr::GetParamFromROS<int>("/ikalibr_imgs_to_bag/downsample_num");
        spdlog::info("downsample images: grab an image per {} images to rosbag", downsampleNum);

        auto filenames = ns_ikalibr::FilesInDir(imgPath);

        // erase filenames that are not images
        if (!nameStampRegexStr.empty()) {
            std::regex fileNameRegex(nameStampRegexStr);
            auto iter = std::remove_if(
                filenames.begin(), filenames.end(), [&fileNameRegex](const std::string &filename) {
                    std::string sName = std::filesystem::path(filename).filename().string();
                    return !std::regex_match(sName, fileNameRegex);
                });
            filenames.erase(iter, filenames.cend());
        }

        std::sort(filenames.begin(), filenames.end());

        if (filenames.empty()) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::CRITICAL,
                                     "there is no any matched image in this directory '{}'.",
                                     imgPath);
        }

        spdlog::info("there are '{}' images in this directory.", filenames.size());
        spdlog::info("the first image is '{}'", filenames.front());
        spdlog::info("the  last image is '{}'", filenames.back());

        auto nameAsStamp = ns_ikalibr::GetParamFromROS<bool>("/ikalibr_imgs_to_bag/name_as_stamp");
        spdlog::info("if use name as time stamp: '{}'", nameAsStamp);

        std::vector<double> timestamps(filenames.size());

        if (nameAsStamp) {
            auto nameToStampScale =
                ns_ikalibr::GetParamFromROS<double>("/ikalibr_imgs_to_bag/name_to_stamp_scale");
            spdlog::info("when using name as time stamp, the scale is: '{}'", nameToStampScale);

            for (int i = 0; i < static_cast<int>(filenames.size()); ++i) {
                std::regex fileNameRegex;

                if (!nameStampRegexStr.empty()) {
                    fileNameRegex = std::regex(nameStampRegexStr);
                } else {
                    fileNameRegex = std::regex(R"((\d+)\..+)");
                }

                std::string sName = std::filesystem::path(filenames.at(i)).filename().string();
                std::smatch smatch;

                if (std::regex_search(sName, smatch, fileNameRegex) && smatch.size() >= 2) {
                    timestamps.at(i) = std::stod(smatch[1]) * nameToStampScale;
                } else {
                    spdlog::warn("filename '{}' was not matched using regex expression!!!", sName);
                }
            }

        } else {
            ros::Time::init();
            auto imgFrequency =
                ns_ikalibr::GetParamFromROS<int>("/ikalibr_imgs_to_bag/img_frequency");
            spdlog::info("if not use name as time stamp, the frequency is: '{}'", imgFrequency);
            const double dt = 1.0 / imgFrequency;

            for (int i = 0; i < static_cast<int>(filenames.size()); ++i) {
                timestamps.at(i) = i * dt;
            }
        }

        auto imgFrequency = ns_ikalibr::GetParamFromROS<int>("/ikalibr_imgs_to_bag/img_frequency");
        spdlog::info("if not use name as time stamp, the frequency is: '{}'", imgFrequency);

        auto dstBag = std::make_unique<rosbag::Bag>();
        dstBag->open(bagPath, rosbag::BagMode::Write);

        for (int i = 0; i < static_cast<int>(filenames.size()); ++i) {
            // downsample
            if (i % downsampleNum != 0) {
                continue;
            }

            auto filename = filenames.at(i);
            auto img = cv::imread(filename, cv::IMREAD_UNCHANGED);
            if (img.empty()) {
                spdlog::warn("invalid image: '{}'!!!", filename);
                continue;
            }

            auto time = timestamps.at(i);
            spdlog::info("filename: '{}', time: '{:.3f}'",
                         std::filesystem::path(filename).filename().string(), time);

            cv_bridge::CvImage cvImage;
            cvImage.image = img;
            cvImage.header.stamp = ros::Time(time);
            cvImage.encoding = encoding;

            sensor_msgs::Image sensorImage;
            cvImage.toImageMsg(sensorImage);
            dstBag->write(imgsTopic, cvImage.header.stamp, sensorImage);
        }

        dstBag->close();
        spdlog::info("images in '{}' have been writen to rosbag as '{}'!", imgPath, bagPath);

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