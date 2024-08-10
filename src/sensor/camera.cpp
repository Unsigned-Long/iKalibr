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

#include "sensor/camera.h"
#include "spdlog/spdlog.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

CameraFrame::CameraFrame(double timestamp, cv::Mat greyImg, cv::Mat colorImg, ns_veta::IndexT id)
    : _timestamp(timestamp),
      _greyImg(std::move(greyImg)),
      _colorImg(std::move(colorImg)),
      _id(id) {
    if (!greyImg.empty() && !colorImg.empty() && greyImg.size() != colorImg.size()) {
        spdlog::warn(
            "the size of grey image ({}x{}) is not the same as the one of color image ({}x{})!",
            greyImg.size().width, greyImg.size().height, colorImg.size().width,
            colorImg.size().height);
    }
}

CameraFrame::Ptr CameraFrame::Create(double timestamp,
                                     const cv::Mat &greyImg,
                                     const cv::Mat &colorImg,
                                     ns_veta::IndexT id) {
    return std::make_shared<CameraFrame>(timestamp, greyImg, colorImg, id);
}

cv::Mat &CameraFrame::GetImage() { return _greyImg; }

double CameraFrame::GetTimestamp() const { return _timestamp; }

void CameraFrame::SetTimestamp(double timestamp) { _timestamp = timestamp; }

std::ostream &operator<<(std::ostream &os, const CameraFrame &frame) {
    os << "image: " << frame._greyImg.size << ", timestamp: " << frame._timestamp;
    return os;
}

void CameraFrame::ReleaseMat() {
    _greyImg.release();
    _colorImg.release();
}

ns_veta::IndexT CameraFrame::GetId() const { return _id; }

void CameraFrame::SetId(ns_veta::IndexT id) { _id = id; }

cv::Mat &CameraFrame::GetColorImage() { return _colorImg; }
}  // namespace ns_ikalibr