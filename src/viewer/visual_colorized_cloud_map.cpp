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

#include "viewer/visual_colorized_cloud_map.h"
#include "calib/calib_param_manager.h"
#include "pcl/kdtree/kdtree_flann.h"
#include "sensor/camera.h"
#include "sensor/rgbd.h"
#include "opencv2/imgproc.hpp"
#include "spdlog/spdlog.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

// -----------------
// ColorizedCloudMap
// -----------------
ColorizedCloudMap::ColorizedCloudMap(const std::vector<CameraFrame::Ptr> &frames,
                                     ns_veta::Veta::Ptr veta,
                                     SplineBundleType::Ptr splines,
                                     ns_veta::PinholeIntrinsic::Ptr intri,
                                     const Sophus::SE3d &SE3_SenToBr,
                                     const double &TO_SenToBr)
    : _frames(frames),
      _veta(std::move(veta)),
      _splines(std::move(splines)),
      _intri(std::move(intri)),
      SE3_CmToBr(SE3_SenToBr),
      TO_CmToBr(TO_SenToBr) {
    for (const auto &frame : _frames) {
        // current frame is not rebuild by SfM
        if (_veta->views.find(frame->GetId()) == _veta->views.cend()) {
            continue;
        }
        // undistorted gray image
        cv::Mat undistImg =
            CalibParamManager::ParIntri::UndistortImage(_intri, frame->GetColorImage());
        _viewIdToFrame.insert({frame->GetId(), {CurCmToW(frame->GetTimestamp()), undistImg}});
    }
}

ColorizedCloudMap::Ptr ColorizedCloudMap::CreateForCameras(
    const std::string &topic,
    const std::vector<CameraFrame::Ptr> &frames,
    const ns_veta::Veta::Ptr &veta,
    const SplineBundleType::Ptr &splines,
    const CalibParamManager::Ptr &parMagr) {
    return std::make_shared<ColorizedCloudMap>(
        frames, veta, splines, parMagr->INTRI.Camera.at(topic), parMagr->EXTRI.SE3_CmToBr(topic),
        parMagr->TEMPORAL.TO_CmToBr.at(topic));
}

ColorizedCloudMap::Ptr ColorizedCloudMap::CreateForRGBDs(const std::string &topic,
                                                         const std::vector<CameraFramePtr> &frames,
                                                         const ns_veta::Veta::Ptr &veta,
                                                         const SplineBundleType::Ptr &splines,
                                                         const CalibParamManagerPtr &parMagr) {
    return std::make_shared<ColorizedCloudMap>(
        frames, veta, splines, parMagr->INTRI.RGBD.at(topic)->intri,
        parMagr->EXTRI.SE3_DnToBr(topic), parMagr->TEMPORAL.TO_DnToBr.at(topic));
}

ColorPointCloud::Ptr ColorizedCloudMap::Colorize(const IKalibrPointCloud::Ptr &cloudMap, int K) {
    PosPointCloud::Ptr lmCloud(new PosPointCloud);
    lmCloud->resize(_veta->structure.size());
    std::map<std::size_t, ns_veta::IndexT> cloudIdxToLMIdx;
    std::size_t index = 0;
    for (const auto &[lmId, lm] : _veta->structure) {
        // store landmark
        PosPoint &p = lmCloud->at(index);
        p.x = static_cast<float>(lm.X(0));
        p.y = static_cast<float>(lm.X(1));
        p.z = static_cast<float>(lm.X(2));
        // record index
        cloudIdxToLMIdx.insert({index, lmId});

        ++index;
    }

    const int width = (int)_intri->imgWidth, height = (int)_intri->imgHeight, padding = 1;
    Eigen::Vector2d leftTop = _intri->ImgToCam(Eigen::Vector2d(0.0 + padding, 0.0 + padding));
    Eigen::Vector2d rightBottom =
        _intri->ImgToCam(Eigen::Vector2d(width - padding, height - padding));

    pcl::KdTreeFLANN<PosPoint> kdtree;
    kdtree.setInputCloud(lmCloud);
    std::vector<int> pointIdxKNNSearch(K);
    std::vector<float> pointKNNSquaredDistance(K);

    ColorPointCloud::Ptr colorMap(new ColorPointCloud);
    colorMap->reserve(cloudMap->size());
    spdlog::info("performing colorizing, this would cost some time...");
#pragma omp parallel for num_threads(omp_get_max_threads()) default(none)                    \
    shared(K, pointIdxKNNSearch, pointKNNSquaredDistance, cloudMap, cloudIdxToLMIdx, kdtree, \
               colorMap, leftTop, rightBottom, width, height)
    for (int i = 0; i < static_cast<int>(cloudMap->points.size()); ++i) {
        // if (i % 100000 == 0) { std::cout << i << '/' << cloudMap->points.size() << std::endl; }
        const auto &ip = cloudMap->points.at(i);

        if (IS_POS_NAN(ip)) {
            continue;
        }

        if (kdtree.nearestKSearchT(ip, K, pointIdxKNNSearch, pointKNNSquaredDistance) == 0) {
            continue;
        }

        Eigen::Vector3d totalBGR = Eigen::Vector3d::Zero();
        int count = 0;

        for (const auto &pId : pointIdxKNNSearch) {
            auto lmId = cloudIdxToLMIdx.at(pId);
            for (const auto &[viewId, feat] : _veta->structure.at(lmId).obs) {
                const auto &[SE3_CurCmToW, undistImg] = _viewIdToFrame.at(viewId);
                if (!SE3_CurCmToW) {
                    continue;
                }

                // transform point to camera frame
                Eigen::Vector3d pInCm =
                    SE3_CurCmToW->inverse() * Eigen::Vector3d(ExpandPCLPointXYZ(ip));
                if (pInCm(2) < 0.1) {
                    continue;
                }

                // project to camera plane
                const double zInv = 1.0 / pInCm(2);
                Eigen::Vector2d pInCamPlane(pInCm(0) * zInv, pInCm(1) * zInv);

                // invalid
                if (pInCamPlane(0) < leftTop(0) || pInCamPlane(0) > rightBottom(0) ||
                    pInCamPlane(1) < leftTop(1) || pInCamPlane(1) > rightBottom(1)) {
                    continue;
                }

                Eigen::Vector2i pixel = _intri->CamToImg(pInCamPlane).cast<int>();

                // invalid
                if (pixel(0) < 0 || pixel(1) < 0 || pixel(0) > width - 1 || pixel(1) > height - 1) {
                    continue;
                }

                // row: pixel(1), col: pixel(0)
                auto bgr = undistImg.at<cv::Vec3b>(pixel(1), pixel(0));
                totalBGR += Eigen::Vector3d(bgr(0), bgr(1), bgr(2));
                ++count;

                // cv::imshow(std::to_string(viewId), DrawPoint(undistImg, pixel));
            }
        }

        if (count == 0) {
            continue;
        }
        ColorPoint op;
        op.x = ip.x, op.y = ip.y, op.z = ip.z;
        // although we should perform color average in the hsv space, to reduce the computation, we
        // directly perform it in BGR space
        Eigen::Vector3d avgBGR = totalBGR / count;
        op.r = static_cast<uchar>(avgBGR(2));
        op.g = static_cast<uchar>(avgBGR(1));
        op.b = static_cast<uchar>(avgBGR(0));
        op.a = 255;
#pragma omp critical
        { colorMap->push_back(op); }
        // cv::waitKey(0);
    }
    return colorMap;
}

std::optional<Sophus::SE3d> ColorizedCloudMap::CurCmToW(double timeByCm) {
    double timeByBr = timeByCm + TO_CmToBr;
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &posSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

    if (!so3Spline.TimeStampInRange(timeByBr) || !posSpline.TimeStampInRange(timeByBr)) {
        return {};
    } else {
        Sophus::SE3d curBrToW(so3Spline.Evaluate(timeByBr), posSpline.Evaluate(timeByBr));
        return curBrToW * SE3_CmToBr;
    }
}

cv::Mat ColorizedCloudMap::DrawPoint(const cv::Mat &img,
                                     const Eigen::Vector2i &feat,
                                     const cv::Scalar &color) {
    cv::Mat res;
    cv::cvtColor(img, res, cv::COLOR_GRAY2BGR);
    // square
    cv::drawMarker(res, cv::Point2d(feat(0), feat(1)), color, cv::MarkerTypes::MARKER_SQUARE, 10,
                   1);
    // key point
    cv::drawMarker(res, cv::Point2d(feat(0), feat(1)), color, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
    return res;
}
}  // namespace ns_ikalibr