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

#include "sensor/rgbd.h"
#include "spdlog/spdlog.h"
#include "opencv2/imgproc.hpp"
#include "sensor/rgbd_intrinsic.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

// ---------
// RGBDFrame
// ---------

RGBDFrame::RGBDFrame(double timestamp,
                     const cv::Mat& greyImg,
                     const cv::Mat& colorImg,
                     cv::Mat depthImg,
                     ns_veta::IndexT id)
    : CameraFrame(timestamp, greyImg, colorImg, id),
      _depthImg(std::move(depthImg)) {
    if (!greyImg.empty() && !depthImg.empty() && greyImg.size() != depthImg.size()) {
        spdlog::warn(
            "the size of grey image ({}x{}) is not the same as the one of depth image ({}x{})!",
            greyImg.size().width, greyImg.size().height, depthImg.size().width,
            depthImg.size().height);
    } else if (!colorImg.empty() && !depthImg.empty() && colorImg.size() != depthImg.size()) {
        spdlog::warn(
            "the size of color image ({}x{}) is not the same as the one of depth image ({}x{})!",
            colorImg.size().width, colorImg.size().height, depthImg.size().width,
            depthImg.size().height);
    }
}

RGBDFrame::Ptr RGBDFrame::Create(double timestamp,
                                 const cv::Mat& greyImg,
                                 const cv::Mat& colorImg,
                                 const cv::Mat& depthImg,
                                 ns_veta::IndexT id) {
    return std::make_shared<RGBDFrame>(timestamp, greyImg, colorImg, depthImg, id);
}

cv::Mat& RGBDFrame::GetDepthImage() { return _depthImg; }

void RGBDFrame::ReleaseMat() {
    CameraFrame::ReleaseMat();
    _depthImg.release();
}

cv::Mat RGBDFrame::CreateColorDepthMap(const RGBDIntrinsics::Ptr& intri,
                                       bool withColorMat,
                                       float zMin,
                                       float zMax) const {
    int rowCnt = _depthImg.rows;
    int colCnt = _depthImg.cols;

    cv::Mat invDepthImg(rowCnt, colCnt, CV_32FC1, cv::Scalar(0.0f));
    float min = std::numeric_limits<float>::max(), max = std::numeric_limits<float>::min();

    for (int row = 0; row < rowCnt; ++row) {
        auto iData = invDepthImg.ptr<float>(row);
        auto dData = _depthImg.ptr<float>(row);
        for (int col = 0; col < colCnt; ++col) {
            auto depth = (float)intri->ActualDepth(dData[0]);
            if (depth > zMin && depth < zMax) {
                iData[0] = 1.0f / depth;
            } else {
                iData[0] = 0.0f;
            }
            if (iData[0] < min) {
                min = iData[0];
            }
            if (iData[0] > max) {
                max = iData[0];
            }
            // color mat ptr
            iData += 1;
            // depth mat ptr
            dData += 1;
        }
    }
    float alpha = 255.0f / (max - min), beta = -min * alpha;
    cv::Mat uCharImg, colorImg;

    cv::convertScaleAbs(invDepthImg, uCharImg, alpha, beta);
    cv::applyColorMap(uCharImg, colorImg, cv::COLORMAP_PLASMA);

    // double min = 0.0, max = 0.0;
    // cv::minMaxIdx(_depthImg, &min, &max);
    // new value = (value - min) * (255.0f / (max - min));
    // double alpha = 255.0 / (max - min), beta = -min * alpha;
    // cv::Mat uCharImg, colorImg;
    // cv::convertScaleAbs(_depthImg, uCharImg, alpha, beta);
    // cv::applyColorMap(uCharImg, colorImg, cv::COLORMAP_PLASMA);

    if (withColorMat) {
        cv::Mat rgbdMap;
        cv::hconcat(_colorImg, colorImg, rgbdMap);
        return rgbdMap;
    } else {
        return colorImg;
    }
}

ColorPointCloud::Ptr RGBDFrame::CreatePointCloud(const RGBDIntrinsicsPtr& intri,
                                                 float zMin,
                                                 float zMax) {
    auto cMat = _colorImg;
    auto dMat = _depthImg;
    int rowCnt = cMat.rows;
    int colCnt = cMat.cols;

    if (cMat.empty() || dMat.empty() || cMat.size != dMat.size) {
        return nullptr;
    }

    ColorPointCloud::Ptr cloud(new ColorPointCloud);
    cloud->reserve(rowCnt * colCnt);
    for (int row = 0; row < rowCnt; ++row) {
        auto cData = cMat.ptr<uchar>(row);
        auto dData = dMat.ptr<float>(row);
        for (int col = 0; col < colCnt; ++col) {
            auto depth = (float)intri->ActualDepth(dData[0]);
            if (depth > zMin && depth < zMax) {
                Eigen::Vector2d lmInDnPlane = intri->intri->ImgToCam({col, row});
                Eigen::Vector3d lmInDn(lmInDnPlane(0) * depth, lmInDnPlane(1) * depth, depth);

                ColorPoint p;
                p.x = (float)lmInDn(0);
                p.y = (float)lmInDn(1);
                p.z = (float)lmInDn(2);
                p.b = cData[0];
                p.g = cData[1];
                p.r = cData[2];
                p.a = 255;
                cloud->push_back(p);
            }
            // color mat ptr
            cData += 3;
            // depth mat ptr
            dData += 1;
        }
    }
    return cloud;
}

IKalibrPointCloud::Ptr RGBDFrame::CreatePointCloud(
    double rsExpFactor, double readout, const RGBDIntrinsicsPtr& intri, float zMin, float zMax) {
    auto cMat = _colorImg;
    auto dMat = _depthImg;
    int rowCnt = cMat.rows;
    int colCnt = cMat.cols;

    if (cMat.empty() || dMat.empty() || cMat.size != dMat.size) {
        return nullptr;
    }

    IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud);
    cloud->reserve(rowCnt * colCnt);
    const int imgHeight = _greyImg.rows;
    for (int row = 0; row < rowCnt; ++row) {
        auto dData = dMat.ptr<float>(row);
        const double rdFactorAry = row / (double)imgHeight - rsExpFactor;
        for (int col = 0; col < colCnt; ++col) {
            auto depth = (float)intri->ActualDepth(dData[0]);
            if (depth > zMin && depth < zMax) {
                Eigen::Vector2d lmInDnPlane = intri->intri->ImgToCam({col, row});
                Eigen::Vector3d lmInDn(lmInDnPlane(0) * depth, lmInDnPlane(1) * depth, depth);

                IKalibrPoint p;
                p.timestamp = _timestamp + rdFactorAry * readout;
                p.x = (float)lmInDn(0);
                p.y = (float)lmInDn(1);
                p.z = (float)lmInDn(2);
                cloud->push_back(p);
            }
            // depth mat ptr
            dData += 1;
        }
    }
    return cloud;
}

// ----------
// DepthFrame
// ----------

DepthFrame::DepthFrame(double timestamp, cv::Mat depthImg, ns_veta::IndexT id)
    : _timestamp(timestamp),
      _depthImg(std::move(depthImg)),
      _id(id) {}

DepthFrame::Ptr DepthFrame::Create(double timestamp, const cv::Mat& depthImg, ns_veta::IndexT id) {
    return std::make_shared<DepthFrame>(timestamp, depthImg, id);
}

double DepthFrame::GetTimestamp() const { return _timestamp; }

void DepthFrame::SetTimestamp(double timestamp) { _timestamp = timestamp; }

std::ostream& operator<<(std::ostream& os, const DepthFrame& frame) {
    os << "image: " << frame._depthImg.size << ", timestamp: " << frame._timestamp;
    return os;
}

void DepthFrame::ReleaseMat() { _depthImg.release(); }

ns_veta::IndexT DepthFrame::GetId() const { return _id; }

void DepthFrame::SetId(ns_veta::IndexT id) { _id = id; }

cv::Mat& DepthFrame::GetDepthImage() { return _depthImg; }

}  // namespace ns_ikalibr