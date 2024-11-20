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

#include "core/visual_distortion.h"
#include "veta/camera/pinhole.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

/**
 * EventUndistortionMap
 */
VisualUndistortionMap::VisualUndistortionMap(const ns_veta::PinholeIntrinsic::Ptr &intri) {
    _eventRemoveDisto.resize(intri->imgWidth, intri->imgHeight);
    for (size_t x = 0; x < intri->imgWidth; x++) {
        for (size_t y = 0; y < intri->imgHeight; y++) {
            _eventRemoveDisto(x, y) = intri->GetUndistoPixel({x, y});
        }
    }
    std::tie(_map1, _map2) = InitUndistortRectifyMap(intri);
}

VisualUndistortionMap::Ptr VisualUndistortionMap::Create(
    const ns_veta::PinholeIntrinsicPtr &intri) {
    return std::make_shared<VisualUndistortionMap>(intri);
}

cv::Mat VisualUndistortionMap::RemoveDistortion(const cv::Mat &distoImg, int interpolation) const {
    cv::Mat undistImg;
    cv::remap(distoImg, undistImg, _map1, _map2, interpolation, CV_HAL_BORDER_CONSTANT);
    return undistImg;
}

cv::Mat VisualUndistortionMap::ObtainKMat(const ns_veta::PinholeIntrinsicPtr &intri) {
    const auto &par = intri->GetParams();
    cv::Mat K = (cv::Mat_<double>(3, 3) << par[0], 0.0, par[2], 0.0, par[1], par[3], 0.0, 0.0, 1.0);
    return K;
}

Eigen::Matrix3d VisualUndistortionMap::ObtainKMatEigen(const ns_veta::PinholeIntrinsicPtr &intri) {
    const auto &par = intri->GetParams();
    Eigen::Matrix3d K;
    K << par[0], 0.0, par[2], 0.0, par[1], par[3], 0.0, 0.0, 1.0;
    return K;
}

cv::Mat VisualUndistortionMap::ObtainDMat(const ns_veta::PinholeIntrinsicPtr &intri) {
    const auto &par = intri->GetParams();
    cv::Mat D;

    if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri)) {
        // k_1, k_2, p_1, p_2[, k_3[, k_4, k_5, k_6[, s_1, s_2, s_3, s_4[, t_x, t_y]]]]
        D = (cv::Mat_<double>(5, 1) << par[4], par[5], par[7], par[8], par[6]);
    } else if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicFisheye>(intri)) {
        // k_1, k_2, k_3, k_4
        D = (cv::Mat_<double>(4, 1) << par[4], par[5], par[6], par[7]);
    } else {
        throw Status(Status::CRITICAL,
                     "unknown camera intrinsic model! supported models:\n"
                     "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                     "(b)  pinhole_fisheye (k1, k2, k3, k4)");
    }
    return D;
}

std::pair<cv::Mat, cv::Mat> VisualUndistortionMap::InitUndistortRectifyMap(
    const ns_veta::PinholeIntrinsic::Ptr &intri) {
    cv::Mat map1, map2;

    auto [K, D] = ObtainKDMatForUndisto(intri);
    cv::Mat E = cv::Mat::eye(3, 3, cv::DataType<double>::type);
    cv::Size size(intri->imgWidth, intri->imgHeight);

    if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri)) {
        cv::initUndistortRectifyMap(K, D, E, K, size, CV_16SC2, map1, map2);
    } else if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicFisheye>(intri)) {
        cv::fisheye::initUndistortRectifyMap(K, D, E, K, size, CV_16SC2, map1, map2);
        // cv::remap(src, undistImg, map1, map2, cv::INTER_LINEAR, CV_HAL_BORDER_CONSTANT);
    } else {
        throw Status(Status::CRITICAL,
                     "unknown camera intrinsic model! supported models:\n"
                     "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                     "(b)  pinhole_fisheye (k1, k2, k3, k4)");
    }
    return {map1, map2};
}

std::pair<cv::Mat, cv::Mat> VisualUndistortionMap::ObtainKDMatForUndisto(
    const ns_veta::PinholeIntrinsic::Ptr &intri) {
    return {ObtainKMat(intri), ObtainDMat(intri)};
}

}  // namespace ns_ikalibr
