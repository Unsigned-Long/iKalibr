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

#include "sensor/depth_data_loader.h"
#include "sensor_msgs/Image.h"
#include "sensor_msgs/CompressedImage.h"
#include "cv_bridge/cv_bridge.h"
#include "util/status.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

DepthDataLoader::DepthDataLoader(CameraModelType model, bool isInverse)
    : _model(model),
      _isInverse(isInverse) {}

DepthDataLoader::Ptr DepthDataLoader::GetLoader(const std::string &modelStr, bool isInverse) {
    // try extract radar model
    CameraModelType model;
    try {
        model = EnumCast::stringToEnum<CameraModelType>(modelStr);
    } catch (...) {
        throw Status(Status::ERROR, CameraModel::UnsupportedCameraModelMsg(modelStr));
    }
    DepthDataLoader::Ptr dataLoader;
    switch (model) {
        case CameraModelType::SENSOR_IMAGE_GS:
        case CameraModelType::SENSOR_IMAGE_RS_FIRST:
        case CameraModelType::SENSOR_IMAGE_RS_MID:
        case CameraModelType::SENSOR_IMAGE_RS_LAST:
            dataLoader = DepthSensorImageLoader::Create(model, isInverse);
            break;
        case CameraModelType::SENSOR_IMAGE_COMP_GS:
        case CameraModelType::SENSOR_IMAGE_COMP_RS_FIRST:
        case CameraModelType::SENSOR_IMAGE_COMP_RS_MID:
        case CameraModelType::SENSOR_IMAGE_COMP_RS_LAST:
            dataLoader = DepthSensorImageCompLoader::Create(model, isInverse);
            break;
        default:
            throw Status(Status::ERROR, CameraModel::UnsupportedCameraModelMsg(modelStr));
    }
    return dataLoader;
}

CameraModelType DepthDataLoader::GetCameraModel() const { return _model; }

void DepthDataLoader::InverseMat(cv::Mat &floatMat) {
    int rowCnt = floatMat.rows;
    int colCnt = floatMat.cols;

    for (int row = 0; row < rowCnt; ++row) {
        auto dData = floatMat.ptr<float>(row);
        for (int col = 0; col < colCnt; ++col) {
            auto depth = dData[0];
            if (depth > 1E-6) {
                // non-zero values
                dData[0] = 1.0f / depth;
            } else {
                // zero values, keep zeros
                dData[0] = 0.0f;
            }
            dData += 1;
        }
    }
}

// ----------------------
// DepthSensorImageLoader
// ----------------------
DepthSensorImageLoader::DepthSensorImageLoader(CameraModelType model, bool isInverse)
    : DepthDataLoader(model, isInverse) {}

DepthSensorImageLoader::Ptr DepthSensorImageLoader::Create(CameraModelType model, bool isInverse) {
    return std::make_shared<DepthSensorImageLoader>(model, isInverse);
}

DepthFrame::Ptr DepthSensorImageLoader::UnpackFrame(const rosbag::MessageInstance &msgInstance) {
    sensor_msgs::ImagePtr msg = msgInstance.instantiate<sensor_msgs::Image>();

    CheckMessage<sensor_msgs::Image>(msg);

    cv::Mat dImg;
    cv_bridge::toCvCopy(msg)->image.copyTo(dImg);
    if (dImg.channels() != 1) {
        throw Status(Status::ERROR, "the channel of depth image dose not equal to 1!!!");
    }
    if (dImg.type() != CV_32FC1) {
        // todo: the alpha=1.0, beta=0.0, assume that the depth factor is not provided
        dImg.convertTo(dImg, CV_32F, 1.0, 0.0);
    }
    if (_isInverse) {
        InverseMat(dImg);
    }

    return DepthFrame::Create(msg->header.stamp.toSec(), dImg);
}

// --------------------------
// DepthSensorImageCompLoader
// --------------------------
DepthSensorImageCompLoader::DepthSensorImageCompLoader(CameraModelType model, bool isInverse)
    : DepthDataLoader(model, isInverse) {}

DepthSensorImageCompLoader::Ptr DepthSensorImageCompLoader::Create(CameraModelType model,
                                                                   bool isInverse) {
    return std::make_shared<DepthSensorImageCompLoader>(model, isInverse);
}

DepthFrame::Ptr DepthSensorImageCompLoader::UnpackFrame(
    const rosbag::MessageInstance &msgInstance) {
    sensor_msgs::CompressedImageConstPtr msg =
        msgInstance.instantiate<sensor_msgs::CompressedImage>();

    CheckMessage<sensor_msgs::CompressedImage>(msg);

    cv::Mat dImg;
    cv_bridge::toCvCopy(msg)->image.copyTo(dImg);
    if (dImg.channels() != 1) {
        throw Status(Status::ERROR, "the channel of depth image dose not equal to 1!!!");
    }
    if (dImg.type() != CV_32FC1) {
        // todo: the alpha=1.0, beta=0.0, assume that the depth factor is not provided
        dImg.convertTo(dImg, CV_32F, 1.0, 0.0);
    }
    if (_isInverse) {
        InverseMat(dImg);
    }
    return DepthFrame::Create(msg->header.stamp.toSec(), dImg);
}

}  // namespace ns_ikalibr
