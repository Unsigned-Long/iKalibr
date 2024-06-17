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

#include "viewer/visual_lidar_covisibility.h"
#include "pcl/common/transforms.h"
#include "opencv2/imgproc.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

VisualLiDARCovisibility::VisualLiDARCovisibility(IKalibrPointCloud::Ptr cloudMap)
    : _cloudMap(std::move(cloudMap)) {}

VisualLiDARCovisibility::Ptr VisualLiDARCovisibility::Create(
    const IKalibrPointCloud::Ptr &cloudMap) {
    return std::make_shared<VisualLiDARCovisibility>(cloudMap);
}

std::pair<cv::Mat, cv::Mat> VisualLiDARCovisibility::CreateCovisibility(
    const Sophus::SE3d &SE3_CurCmToW,
    const ns_veta::PinholeIntrinsic::Ptr &intri,
    float zMin,
    float zMax) {
    IKalibrPointCloud::Ptr mapInCam(new IKalibrPointCloud);
    pcl::transformPointCloud(*_cloudMap, *mapInCam, SE3_CurCmToW.inverse().matrix().cast<float>());

    const int width = (int)intri->imgWidth, height = (int)intri->imgHeight, padding = 1;
    Eigen::Vector2d leftTop = intri->ImgToCam(Eigen::Vector2d(0.0 + padding, 0.0 + padding));
    Eigen::Vector2d rightBottom =
        intri->ImgToCam(Eigen::Vector2d(width - padding, height - padding));

    float min = std::numeric_limits<float>::max(), max = std::numeric_limits<float>::min();
    cv::Mat invDepthImg(height, width, CV_32FC1, cv::Scalar(0.0f));
    const float zMaxInv = 1.0f / zMax;

    for (const auto &p : mapInCam->points) {
        if (IS_POS_NAN(p) || p.z < zMin || p.z > zMax) {
            continue;
        }

        const float zInv = 1.0f / p.z;
        Eigen::Vector2d pInCamPlane(p.x * zInv, p.y * zInv);

        // invalid
        if (pInCamPlane(0) < leftTop(0) || pInCamPlane(0) > rightBottom(0) ||
            pInCamPlane(1) < leftTop(1) || pInCamPlane(1) > rightBottom(1)) {
            continue;
        }

        Eigen::Vector2i pixel = intri->CamToImg(pInCamPlane).cast<int>();

        // invalid
        if (pixel(0) < 0 || pixel(1) < 0 || pixel(0) > width - 1 || pixel(1) > height - 1) {
            continue;
        }

        // row: pixel(1), col: pixel(0)
        auto &val = invDepthImg.at<float>(pixel(1), pixel(0));

        if (val < zMaxInv) {
            // this pos has no value
            val = zInv;
            if (val < min) {
                min = val;
            }
            if (val > max) {
                max = val;
            }
        } else {
            // this pos has value
            if (val < zInv) {
                // the new pose is closer
                val = zInv;
                if (val < min) {
                    min = val;
                }
                if (val > max) {
                    max = val;
                }
            }
        }
    }

    // invDepthImg = (invDepthImg - min) * (255.0f / (max - min));
    float alpha = 255.0f / (max - min), beta = -min * alpha;
    cv::Mat uCharImg, colorImg;
    cv::convertScaleAbs(invDepthImg, uCharImg, alpha, beta);
    cv::applyColorMap(uCharImg, colorImg, cv::COLORMAP_PLASMA);

    return {invDepthImg, colorImg};
}
}  // namespace ns_ikalibr