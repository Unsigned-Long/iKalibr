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

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;

class VisualPixelDynamic {
public:
    using Ptr = std::shared_ptr<VisualPixelDynamic>;
    static constexpr int MID = 1;

protected:
    // camera frame, pixel (undistorted)
    std::array<std::pair<CameraFramePtr, Eigen::Vector2d>, 3> _movement;
    Eigen::Vector2d _midPointVel;

public:
    explicit VisualPixelDynamic(
        const std::array<std::pair<CameraFramePtr, Eigen::Vector2d>, 3>& movement);

    static Ptr Create(const std::array<std::pair<CameraFramePtr, Eigen::Vector2d>, 3>& movement);

    [[nodiscard]] cv::Mat CreatePixelDynamicMat(const ns_veta::PinholeIntrinsic::Ptr& intri) const;

    [[nodiscard]] const Eigen::Vector2d& GetMidPointVel() const;

    [[nodiscard]] const Eigen::Vector2d& GetMidPoint() const;

    [[nodiscard]] const CameraFramePtr& GetMidCameraFrame() const;

protected:
    // given three points, compute the first order of the middle point using lagrange polynomial
    static double MidLagrangePolynomialFOD(const std::array<std::pair<double, double>, 3>& data);

    [[nodiscard]] Eigen::Vector2d MidLagrangePolynomialFOD() const;

    static cv::Mat GetInRangeSubMat(const cv::Mat& img, const Eigen::Vector2d& p, int padding);

    static cv::Mat DrawKeypoint(cv::Mat img, const Eigen::Vector2d& p);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_VISUAL_PIXEL_DYNAMIC_H
