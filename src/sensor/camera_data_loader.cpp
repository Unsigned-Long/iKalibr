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

#include "sensor/camera_data_loader.h"
#include "sensor_msgs/Image.h"
#include "util/status.hpp"
#include "spdlog/fmt/fmt.h"

_3_

namespace ns_ikalibr {

    CameraDataLoader::CameraDataLoader(CameraModelType model) : _model(model) {}

    CameraDataLoader::Ptr CameraDataLoader::GetLoader(const std::string &modelStr) {
        // try extract radar model
        CameraModelType model;
        try {
            model = EnumCast::stringToEnum<CameraModelType>(modelStr);
        } catch (...) {
            throw Status(
                    Status::WARNING,
                    "Unsupported camera Type: '{}'. "
                    "Currently supported camera types are: \n"
                    "1. SENSOR_IMAGE_GS: https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                    "2. SENSOR_IMAGE_RS_FIRST: first-row exposure, https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                    "3. SENSOR_IMAGE_RS_MID: middle-row exposure, https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                    "4. SENSOR_IMAGE_RS_LAST: last-row exposure, https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                    "...\n"
                    "If you need to use other camera types, "
                    "please 'Issues' us on the profile of the github repository.",
                    modelStr
            );
        }
        CameraDataLoader::Ptr imuDataLoader;
        switch (model) {
            case CameraModelType::SENSOR_IMAGE_GS:
            case CameraModelType::SENSOR_IMAGE_RS_FIRST:
            case CameraModelType::SENSOR_IMAGE_RS_MID:
            case CameraModelType::SENSOR_IMAGE_RS_LAST:
                imuDataLoader = SensorImageLoader::Create(model);
                break;
            default:
                throw Status(
                        Status::WARNING,
                        "Unsupported camera Type: '{}'. "
                        "Currently supported camera types are: \n"
                        "1. SENSOR_IMAGE_GS: https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                        "2. SENSOR_IMAGE_RS_FIRST: first-row exposure, https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                        "3. SENSOR_IMAGE_RS_MID: middle-row exposure, https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                        "4. SENSOR_IMAGE_RS_LAST: last-row exposure, https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
                        "...\n"
                        "If you need to use other camera types, "
                        "please 'Issues' us on the profile of the github repository.",
                        modelStr
                );
        }
        return imuDataLoader;
    }

    CameraModelType CameraDataLoader::GetCameraModel() const {
        return _model;
    }

    // -----------------
    // SensorImageLoader
    // -----------------
    SensorImageLoader::SensorImageLoader(CameraModelType model) : CameraDataLoader(model) {}

    SensorImageLoader::Ptr SensorImageLoader::Create(CameraModelType model) {
        return std::make_shared<SensorImageLoader>(model);
    }

    CameraFrame::Ptr SensorImageLoader::UnpackFrame(const rosbag::MessageInstance &msgInstance) {
        // imu data item
        sensor_msgs::ImageConstPtr msg = msgInstance.instantiate<sensor_msgs::Image>();

        CheckMessage<sensor_msgs::Image>(msg);

        cv::Mat img, gImg, cImg;
        cv_bridge::toCvCopy(msg, msg->encoding)->image.copyTo(img);
        if (msg->encoding == sensor_msgs::image_encodings::BGR8 ||
            msg->encoding == sensor_msgs::image_encodings::TYPE_8UC3) {
            // color image
            cImg = img;
            cv::cvtColor(img, gImg, cv::COLOR_BGR2GRAY);
        } else if (msg->encoding == sensor_msgs::image_encodings::MONO8 ||
                   msg->encoding == sensor_msgs::image_encodings::TYPE_8UC1) {
            // grey image
            gImg = img;
            cv::cvtColor(img, cImg, cv::COLOR_GRAY2BGR);
        } else {
            throw Status(
                    Status::CRITICAL,
                    "Unsupported sensor_msgs::image_encodings type: '{}'."
                    "Only 'sensor_msgs::image_encodings::BGR8|TYPE_8UC3' and 'sensor_msgs::image_encodings::MONO8|TYPE_8UC1' "
                    "are supported currently! To decode more type images, please contact us or add code "
                    "by yourself at: line-'{}' in file-'{}'",
                    msg->encoding, __LINE__, __FILE__
            );
        }

        return CameraFrame::Create(msg->header.stamp.toSec(), gImg, cImg);
    }

}