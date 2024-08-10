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