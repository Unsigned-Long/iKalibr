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
#include "sensor_msgs/CompressedImage.h"
#include "cv_bridge/cv_bridge.h"
#include "util/status.hpp"
#include "spdlog/fmt/fmt.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

CameraDataLoader::CameraDataLoader(CameraModelType model)
    : _model(model) {}

CameraDataLoader::Ptr CameraDataLoader::GetLoader(const std::string &modelStr) {
    // try extract radar model
    CameraModelType model;
    try {
        model = EnumCast::stringToEnum<CameraModelType>(modelStr);
    } catch (...) {
        throw Status(Status::ERROR, CameraModel::UnsupportedCameraModelMsg(modelStr));
    }
    CameraDataLoader::Ptr dataLoader;
    switch (model) {
        case CameraModelType::SENSOR_IMAGE_GS:
        case CameraModelType::SENSOR_IMAGE_RS_FIRST:
        case CameraModelType::SENSOR_IMAGE_RS_MID:
        case CameraModelType::SENSOR_IMAGE_RS_LAST:
            dataLoader = SensorImageLoader::Create(model);
            break;
        case CameraModelType::SENSOR_IMAGE_COMP_GS:
        case CameraModelType::SENSOR_IMAGE_COMP_RS_FIRST:
        case CameraModelType::SENSOR_IMAGE_COMP_RS_MID:
        case CameraModelType::SENSOR_IMAGE_COMP_RS_LAST:
            dataLoader = SensorImageCompLoader::Create(model);
            break;
        default:
            throw Status(Status::ERROR, CameraModel::UnsupportedCameraModelMsg(modelStr));
    }
    return dataLoader;
}

CameraModelType CameraDataLoader::GetCameraModel() const { return _model; }

void CameraDataLoader::RefineImgMsgWrongEncoding(const sensor_msgs::Image::Ptr &msg) {
    if (msg->encoding == sensor_msgs::image_encodings::TYPE_8UC1) {
        msg->encoding = sensor_msgs::image_encodings::MONO8;
        static bool warn = false;
        if (!warn) {
            spdlog::warn("encoding type of images is wrong: '{}', change encoding type to '{}'",
                         sensor_msgs::image_encodings::TYPE_8UC1, msg->encoding);
            warn = true;
        }
    } else if (msg->encoding == sensor_msgs::image_encodings::TYPE_16UC1) {
        msg->encoding = sensor_msgs::image_encodings::MONO16;
        static bool warn = false;
        if (!warn) {
            spdlog::warn("encoding type of images is wrong: '{}', change encoding type to '{}'",
                         sensor_msgs::image_encodings::TYPE_16UC1, msg->encoding);
            warn = true;
        }
    } else if (msg->encoding == sensor_msgs::image_encodings::TYPE_8UC3) {
        msg->encoding = sensor_msgs::image_encodings::BGR8;
        static bool warn = false;
        if (!warn) {
            spdlog::warn("encoding type of images is wrong: '{}', change encoding type to '{}'",
                         sensor_msgs::image_encodings::TYPE_8UC3, msg->encoding);
            warn = true;
        }
    } else if (msg->encoding == sensor_msgs::image_encodings::TYPE_8UC4) {
        msg->encoding = sensor_msgs::image_encodings::BGRA8;
        static bool warn = false;
        if (!warn) {
            spdlog::warn("encoding type of images is wrong: '{}', change encoding type to '{}'",
                         sensor_msgs::image_encodings::TYPE_8UC4, msg->encoding);
            warn = true;
        }
    } else if (msg->encoding == sensor_msgs::image_encodings::TYPE_16UC3) {
        msg->encoding = sensor_msgs::image_encodings::BGR16;
        static bool warn = false;
        if (!warn) {
            spdlog::warn("encoding type of images is wrong: '{}', change encoding type to '{}'",
                         sensor_msgs::image_encodings::TYPE_16UC3, msg->encoding);
            warn = true;
        }
    } else if (msg->encoding == sensor_msgs::image_encodings::TYPE_16UC4) {
        msg->encoding = sensor_msgs::image_encodings::BGRA16;
        static bool warn = false;
        if (!warn) {
            spdlog::warn("encoding type of images is wrong: '{}', change encoding type to '{}'",
                         sensor_msgs::image_encodings::TYPE_16UC4, msg->encoding);
            warn = true;
        }
    }
}

// -----------------
// SensorImageLoader
// -----------------
SensorImageLoader::SensorImageLoader(CameraModelType model)
    : CameraDataLoader(model) {}

SensorImageLoader::Ptr SensorImageLoader::Create(CameraModelType model) {
    return std::make_shared<SensorImageLoader>(model);
}

CameraFrame::Ptr SensorImageLoader::UnpackFrame(const rosbag::MessageInstance &msgInstance) {
    sensor_msgs::ImagePtr msg = msgInstance.instantiate<sensor_msgs::Image>();

    CheckMessage<sensor_msgs::Image>(msg);
    RefineImgMsgWrongEncoding(msg);

    cv::Mat cImg, gImg;
    cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::BGR8)->image.copyTo(cImg);
    cv::cvtColor(cImg, gImg, cv::COLOR_BGR2GRAY);

    return CameraFrame::Create(msg->header.stamp.toSec(), gImg, cImg);
}

// ---------------------
// SensorImageCompLoader
// ---------------------
SensorImageCompLoader::SensorImageCompLoader(CameraModelType model)
    : CameraDataLoader(model) {}

SensorImageCompLoader::Ptr SensorImageCompLoader::Create(CameraModelType model) {
    return std::make_shared<SensorImageCompLoader>(model);
}

CameraFrame::Ptr SensorImageCompLoader::UnpackFrame(const rosbag::MessageInstance &msgInstance) {
    sensor_msgs::CompressedImageConstPtr msg =
        msgInstance.instantiate<sensor_msgs::CompressedImage>();

    CheckMessage<sensor_msgs::CompressedImage>(msg);

    cv::Mat cImg, gImg;
    cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::BGR8)->image.copyTo(cImg);
    cv::cvtColor(cImg, gImg, cv::COLOR_BGR2GRAY);

    return CameraFrame::Create(msg->header.stamp.toSec(), gImg, cImg);
}
}  // namespace ns_ikalibr