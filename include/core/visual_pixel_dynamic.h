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

#ifndef IKALIBR_VISUAL_PIXEL_DYNAMIC_H
#define IKALIBR_VISUAL_PIXEL_DYNAMIC_H

#include "util/utils.h"
#include "opencv4/opencv2/core.hpp"
#include "veta/camera/pinhole.h"
#include "sensor/camera_data_loader.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct RGBDVelocityCorr;
using RGBDVelocityCorrPtr = std::shared_ptr<RGBDVelocityCorr>;
struct RGBDIntrinsics;
using RGBDIntrinsicsPtr = std::shared_ptr<RGBDIntrinsics>;

class VisualPixelDynamic {
public:
    using Ptr = std::shared_ptr<VisualPixelDynamic>;
    static constexpr int MID = 1;

protected:
    // camera frame, pixel (undistorted)
    std::array<std::pair<CameraFrame::Ptr, Eigen::Vector2d>, 3> _movement;

public:
    explicit VisualPixelDynamic(
        const std::array<std::pair<CameraFrame::Ptr, Eigen::Vector2d>, 3>& movement);

    static Ptr Create(const std::array<std::pair<CameraFrame::Ptr, Eigen::Vector2d>, 3>& movement);

    [[nodiscard]] const CameraFrame::Ptr& GetMidCameraFrame() const;

    // for rgbd cameras whose have depth images
    [[nodiscard]] RGBDVelocityCorrPtr CreateRGBDVelocityCorr(const RGBDIntrinsicsPtr& intri,
                                                             const CameraModelType& type) const;

    // visualization
    [[nodiscard]] cv::Mat CreatePixelDynamicMat(const ns_veta::PinholeIntrinsic::Ptr& intri,
                                                const Eigen::Vector2d& midVel) const;

protected:
    static cv::Mat GetInRangeSubMat(const cv::Mat& img, const Eigen::Vector2d& p, int padding);

    void DrawTrace(cv::Mat& img, const Eigen::Vector2d& midVel, int pixelDist = 5) const;
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_VISUAL_PIXEL_DYNAMIC_H
