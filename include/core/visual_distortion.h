// iKalibr: Unified Targetless Spatiotemporal Calibration Framework
// Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China
// https://github.com/Unsigned-Long/iKalibr.git
// Author: Shuolong Chen (shlchen@whu.edu.cn)
// GitHub: https://github.com/Unsigned-Long
//  ORCID: 0000-0002-5283-9057
// Purpose: See .h/.hpp file.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * The names of its contributors can not be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
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

#ifndef VISUAL_DISTORTION_H
#define VISUAL_DISTORTION_H

#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_veta {
struct PinholeIntrinsic;
using PinholeIntrinsicPtr = std::shared_ptr<PinholeIntrinsic>;
struct PinholeIntrinsicBrownT2;
struct PinholeIntrinsicFisheye;
}  // namespace ns_veta

namespace ns_ikalibr {

struct VisualUndistortionMap {
public:
    using Ptr = std::shared_ptr<VisualUndistortionMap>;

private:
    // operator on single pixel
    Eigen::Array<Eigen::Vector2d, Eigen::Dynamic, Eigen::Dynamic> _dataRemoveDisto, _dataAddDisto;
    // operate on entire image (remove disto)
    cv::Mat _map1, _map2;

public:
    VisualUndistortionMap(const ns_veta::PinholeIntrinsicPtr& intri);

    static Ptr Create(const ns_veta::PinholeIntrinsicPtr& intri);

    // for events
    std::pair<double, double> RemoveDistortion(const int& x, const int& y) const;

    // for events
    std::pair<double, double> AddDistortion(const int& x, const int& y) const;

    cv::Mat RemoveDistortion(const cv::Mat& distoImg, int interpolation = cv::INTER_LINEAR) const;

protected:
    static std::pair<cv::Mat, cv::Mat> ObtainKDMatForUndisto(
        const ns_veta::PinholeIntrinsicPtr& intri);

    static std::pair<cv::Mat, cv::Mat> InitUndistortRectifyMap(
        const ns_veta::PinholeIntrinsicPtr& intri);
};
}  // namespace ns_ikalibr

#endif  // VISUAL_DISTORTION_H
