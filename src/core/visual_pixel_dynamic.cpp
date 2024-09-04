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

#include "core/visual_pixel_dynamic.h"
#include "sensor/camera.h"
#include "sensor/rgbd.h"
#include "opencv2/imgproc.hpp"
#include "calib/calib_param_manager.h"
#include "factor/rgbd_velocity_factor.hpp"
#include "sensor/rgbd_intrinsic.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
VisualPixelDynamic::VisualPixelDynamic(
    const std::array<std::pair<CameraFrame::Ptr, Eigen::Vector2d>, 3>& movement)
    : _movement(movement) {}

VisualPixelDynamic::Ptr VisualPixelDynamic::Create(
    const std::array<std::pair<CameraFrame::Ptr, Eigen::Vector2d>, 3>& movement) {
    return std::make_shared<VisualPixelDynamic>(movement);
}

const CameraFrame::Ptr& VisualPixelDynamic::GetMidCameraFrame() const {
    return _movement.at(MID).first;
}

RGBDVelocityCorr::Ptr VisualPixelDynamic::CreateRGBDVelocityCorr(const RGBDIntrinsics::Ptr& intri,
                                                                 double rsExposureFactor,
                                                                 bool rawDepth) const {
    std::array<double, 3> timeAry{}, xAry{}, yAry{};
    for (int i = 0; i < 3; ++i) {
        timeAry[i] = _movement[i].first->GetTimestamp();
        xAry[i] = _movement[i].second(0);
        yAry[i] = _movement[i].second(1);
    }
    const auto& midFrame = _movement.at(MID).first;
    const Eigen::Vector2d& midPoint = _movement.at(MID).second;
    double depth = -1.0;
    if (auto rgbdFrame = std::dynamic_pointer_cast<RGBDFrame>(midFrame); rgbdFrame) {
        if (auto depthMat = rgbdFrame->GetDepthImage(); !depthMat.empty()) {
            if (rawDepth) {
                depth = depthMat.at<float>((int)midPoint(1), (int)midPoint(0));
            } else {
                depth = intri->ActualDepth(depthMat.at<float>((int)midPoint(1), (int)midPoint(0)));
            }
        }
    }
    auto depthMat = std::dynamic_pointer_cast<RGBDFrame>(midFrame)->GetDepthImage();
    return RGBDVelocityCorr::Create(timeAry, xAry, yAry, depth, midFrame, rsExposureFactor);
}

// -------------
// visualization
// -------------
cv::Mat VisualPixelDynamic::CreatePixelDynamicMat(const ns_veta::PinholeIntrinsic::Ptr& intri,
                                                  const Eigen::Vector2d& midVel) const {
    std::array<cv::Mat, 3> imgs;
    // obtain images
    for (int i = 0; i < static_cast<int>(_movement.size()); ++i) {
        imgs[i] = CalibParamManager::ParIntri::UndistortImage(
            intri, _movement.at(i).first->GetColorImage());
    }
    // trace of point
    DrawTrace(imgs[MID], midVel, 2);

    // draw point in each image
    for (int i = 0; i < static_cast<int>(_movement.size()); ++i) {
        const auto& [frame, feat] = _movement.at(i);
        // draw point
        DrawKeypointOnCVMat(imgs[i], feat);
        // put index text
        PutTextOnCVMat(imgs[i], fmt::format("{}", i), feat);
        // draw point in the middle image
        DrawKeypointOnCVMat(imgs[MID], feat);
        // put index text in the middle image
        PutTextOnCVMat(imgs[MID], fmt::format("{}", i), feat);
    }
    // draw mid-point dynamics (points in mid image)
    const auto& [midFrame, midFeat] = _movement.at(MID);
    Eigen::Vector2d endPixel = midFeat + 0.1 * midVel;
    DrawLineOnCVMat(imgs[MID], midFeat, endPixel);

    // concat images
    for (int i = 0; i < static_cast<int>(_movement.size()); ++i) {
        imgs[i] = GetInRangeSubMat(imgs[i], _movement.at(MID).second, 120);
    }
    cv::Mat img;
    cv::hconcat(imgs[0], imgs[1], img);
    cv::hconcat(img, imgs[2], img);
    return img;
}

cv::Mat VisualPixelDynamic::GetInRangeSubMat(const cv::Mat& image,
                                             const Eigen::Vector2d& p,
                                             int padding) {
    int x = (int)p(0);
    int y = (int)p(1);

    int x1 = std::max(0, x - padding);
    int y1 = std::max(0, y - padding);
    int x2 = std::min(image.cols, x + padding);
    int y2 = std::min(image.rows, y + padding);

    cv::Mat sub_image = image(cv::Rect(x1, y1, x2 - x1, y2 - y1));

    int dx = (2 * padding - (x2 - x1)) / 2;
    int dy = (2 * padding - (y2 - y1)) / 2;

    cv::Mat result(2 * padding, 2 * padding, image.type(), cv::Scalar(0, 0, 0));
    sub_image.copyTo(result(cv::Rect(dx, dy, x2 - x1, y2 - y1)));

    return result;
}

void VisualPixelDynamic::DrawTrace(cv::Mat& img,
                                   const Eigen::Vector2d& midVel,
                                   const int pixelDist) const {
    const double duration =
        _movement.back().first->GetTimestamp() - _movement.front().first->GetTimestamp();
    const double sTime = _movement.front().first->GetTimestamp() - duration * 0.25;
    const double eTime = _movement.back().first->GetTimestamp() + duration * 0.25;
    const double deltaTime = pixelDist / midVel.norm();
    std::array<double, 3> tData{};
    std::array<double, 3> xData{};
    std::array<double, 3> yData{};
    for (int i = 0; i < 3; ++i) {
        tData[i] = _movement[i].first->GetTimestamp();
        xData[i] = _movement[i].second(0);
        yData[i] = _movement[i].second(1);
    }
    double xLast = LagrangePolynomial<double, 3>(sTime, tData, xData);
    double yLast = LagrangePolynomial<double, 3>(sTime, tData, yData);
    for (double t = sTime + deltaTime; t < eTime;) {
        double x = LagrangePolynomial<double, 3>(t, tData, xData);
        double y = LagrangePolynomial<double, 3>(t, tData, yData);
        DrawLineOnCVMat(img, cv::Point2d(xLast, yLast), cv::Point2d(x, y), cv::Scalar(0, 0, 255));
        t += deltaTime;
        xLast = x;
        yLast = y;
    }
}
}  // namespace ns_ikalibr