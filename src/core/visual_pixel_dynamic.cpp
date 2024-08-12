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

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
VisualPixelDynamic::VisualPixelDynamic(
    const std::array<std::pair<CameraFrame::Ptr, Eigen::Vector2d>, 3>& movement)
    : _movement(movement),
      _midPointVel(MidLagrangePolynomialFOD()) {}

VisualPixelDynamic::Ptr VisualPixelDynamic::Create(
    const std::array<std::pair<CameraFramePtr, Eigen::Vector2d>, 3>& movement) {
    return std::make_shared<VisualPixelDynamic>(movement);
}

double VisualPixelDynamic::MidLagrangePolynomialFOD(
    const std::array<std::pair<double, double>, 3>& data) {
    double x0 = data[0].first;
    double x1 = data[1].first;
    double x2 = data[2].first;
    double y0 = data[0].second;
    double y1 = data[1].second;
    double y2 = data[2].second;
    double v1 = y0 * (x1 - x2) / (x0 - x1) / (x0 - x2);
    double v2 = y1 * (1 / (x1 - x0) + 1 / (x1 - x2));
    double v3 = y2 * (x1 - x0) / (x2 - x0) / (x2 - x1);
    return v1 + v2 + v3;
}

Eigen::Vector2d VisualPixelDynamic::MidLagrangePolynomialFOD() const {
    std::array<std::pair<double, double>, 3> xAry;
    std::array<std::pair<double, double>, 3> yAry;
    for (int i = 0; i < 3; ++i) {
        xAry[i].first = _movement[i].first->GetTimestamp();
        xAry[i].second = _movement[i].second(0);

        yAry[i].first = _movement[i].first->GetTimestamp();
        yAry[i].second = _movement[i].second(1);
    }
    return {MidLagrangePolynomialFOD(xAry), MidLagrangePolynomialFOD(yAry)};
}

cv::Mat VisualPixelDynamic::CreatePixelDynamicMat(
    const ns_veta::PinholeIntrinsic::Ptr& intri) const {
    std::array<cv::Mat, 3> imgs;
    for (int i = 0; i < static_cast<int>(_movement.size()); ++i) {
        imgs[i] = CalibParamManager::ParIntri::UndistortImage(
            intri, _movement.at(i).first->GetColorImage());
    }
    for (int i = 0; i < static_cast<int>(_movement.size()); ++i) {
        const auto& [frame, feat] = _movement.at(i);
        // pixel
        DrawKeypoint(imgs[i], feat);
        // index
        cv::putText(imgs[i], fmt::format("{}", i), cv::Point2d(feat(0) + 10, feat(1)),
                    cv::HersheyFonts::FONT_HERSHEY_PLAIN, 1.0, cv::Scalar(255, 0, 0), 2);
        // depth
        if (auto depthFrame = std::dynamic_pointer_cast<RGBDFrame>(frame); depthFrame) {
            // img.at<type>(row, col)
            float depth = depthFrame->GetDepthImage().at<float>((int)feat(1), (int)feat(0));
            cv::putText(imgs[i], fmt::format("d: {:.1f}", depth),
                        cv::Point2d(feat(0) + 10, feat(1) - 20),
                        cv::HersheyFonts::FONT_HERSHEY_PLAIN, 1.0, cv::Scalar(255, 0, 0), 2);
        }
    }
    for (int i = 0; i < static_cast<int>(_movement.size()); ++i) {
        const auto& [frame, feat] = _movement.at(i);
        // pixel
        DrawKeypoint(imgs[MID], feat);
        // index
        cv::putText(imgs[MID], std::to_string(i), cv::Point2d(feat(0) + 10, feat(1)),
                    cv::HersheyFonts::FONT_HERSHEY_PLAIN, 1.0, cv::Scalar(255, 0, 0), 2);
        if (i == MID) {
            // mid-point velocity
            Eigen::Vector2d endPixel = feat - 0.1 * _midPointVel;
            cv::line(imgs[MID], cv::Point2d(feat(0), feat(1)),
                     cv::Point2d(endPixel(0), endPixel(1)), cv::Scalar(0, 255, 0), 1, cv::LINE_AA);
        }
    }
    for (int i = 0; i < static_cast<int>(_movement.size()); ++i) {
        imgs[i] = GetInRangeSubMat(imgs[i], _movement.at(i).second, 200);
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

cv::Mat VisualPixelDynamic::DrawKeypoint(cv::Mat img, const Eigen::Vector2d& feat) {
    // square
    cv::drawMarker(img, cv::Point2d(feat(0), feat(1)), cv::Scalar(0, 255, 0),
                   cv::MarkerTypes::MARKER_SQUARE, 10, 1);
    // key point
    cv::drawMarker(img, cv::Point2d(feat(0), feat(1)), cv::Scalar(0, 255, 0),
                   cv::MarkerTypes::MARKER_SQUARE, 2, 2);
    return img;
}

const Eigen::Vector2d& VisualPixelDynamic::GetMidPointVel() const { return _midPointVel; }

const CameraFramePtr& VisualPixelDynamic::GetMidCameraFrame() const {
    return _movement.at(MID).first;
}
}  // namespace ns_ikalibr