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

#include "viewer/visual_gravity.h"
#include "calib/calib_param_manager.h"
#include "opencv2/imgproc.hpp"
#include "sensor/camera.h"
#include "factor/data_correspondence.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

// -------------------
// VisualGravityDrawer
// -------------------

VisualGravityDrawer::VisualGravityDrawer(std::string topic,
                                         ns_veta::Veta::Ptr veta,
                                         SplineBundleType::Ptr splines,
                                         const CalibParamManager::Ptr &parMagr)
    : _topic(std::move(topic)),
      _veta(std::move(veta)),
      _splines(std::move(splines)) {
    _intri = parMagr->INTRI.Camera.at(_topic);
    SE3_CmToBr = parMagr->EXTRI.SE3_CmToBr(_topic);
    TO_CmToBr = parMagr->TEMPORAL.TO_CmToBr.at(_topic);
    GRAVITY = parMagr->GRAVITY;
}

VisualGravityDrawer::Ptr VisualGravityDrawer::Create(const std::string &topic,
                                                     const ns_veta::Veta::Ptr &veta,
                                                     const SplineBundleType::Ptr &splines,
                                                     const CalibParamManager::Ptr &parMagr) {
    return std::make_shared<VisualGravityDrawer>(topic, veta, splines, parMagr);
}

cv::Mat VisualGravityDrawer::CreateGravityImg(const CameraFrame::Ptr &frame, float scale) {
    // undistorted gray image
    cv::Mat undistImgColor, res;
    undistImgColor = CalibParamManager::ParIntri::UndistortImage(_intri, frame->GetColorImage());

    // compute timestamp by reference IMU, we do not consider the readout time for RS cameras here
    double timeByBr = frame->GetTimestamp() + TO_CmToBr;
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &posSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    if (!so3Spline.TimeStampInRange(timeByBr) || !posSpline.TimeStampInRange(timeByBr)) {
        return undistImgColor;
    }

    auto SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
    Eigen::Vector3d POS_BrInBr0 = posSpline.Evaluate(timeByBr);
    Sophus::SE3d SE3_BrToBr0(SO3_BrToBr0, POS_BrInBr0);
    auto SE3_Br0ToCm = (SE3_BrToBr0 * SE3_CmToBr).inverse();

    for (const auto &[lmId, lm] : _veta->structure) {
        auto iter = lm.obs.find(frame->GetId());
        if (iter == lm.obs.cend()) {
            continue;
        }

        Eigen::Vector3d lmInCm = SE3_Br0ToCm * lm.X;
        Eigen::Vector3d gravityInCm = SE3_Br0ToCm.so3() * GRAVITY;

        Eigen::Vector3d end = lmInCm + scale * gravityInCm;

        if (end(2) < 1E-3 || lmInCm(2) < 1E-3) {
            continue;
        }

        Eigen::Vector2d endPixel = _intri->CamToImg({end(0) / end(2), end(1) / end(2)});

        // we do not use extracted raw feature here to keep better consistency
        Eigen::Vector2d feat = _intri->CamToImg({lmInCm(0) / lmInCm(2), lmInCm(1) / lmInCm(2)});

        DrawKeypointOnCVMat(undistImgColor, feat);
        DrawLineOnCVMat(undistImgColor, feat, endPixel);
    }
    return undistImgColor;
}

// -----------------------
// RGBDVisualGravityDrawer
// -----------------------

RGBDVisualGravityDrawer::RGBDVisualGravityDrawer(const std::string &topic,
                                                 const std::vector<OpticalFlowCorr::Ptr> &corrs,
                                                 SplineBundleType::Ptr splines,
                                                 const CalibParamManager::Ptr &parMagr)
    : _splines(std::move(splines)) {
    _intri = parMagr->INTRI.RGBD.at(topic);
    SE3_DnToBr = parMagr->EXTRI.SE3_DnToBr(topic);
    TO_DnToBr = parMagr->TEMPORAL.TO_DnToBr.at(topic);
    GRAVITY = parMagr->GRAVITY;

    for (const auto &corr : corrs) {
        if (_intri->ActualDepth(corr->depth) > 1E-3 /* 1 mm */) {
            // a valid depth
            this->_velCorrs[corr->frame->GetId()].emplace_back(corr);
        }
    }
}

RGBDVisualGravityDrawer::Ptr RGBDVisualGravityDrawer::Create(
    const std::string &topic,
    const std::vector<OpticalFlowCorr::Ptr> &corrs,
    const SplineBundleType::Ptr &splines,
    const CalibParamManager::Ptr &parMagr) {
    return std::make_shared<RGBDVisualGravityDrawer>(topic, corrs, splines, parMagr);
}

cv::Mat RGBDVisualGravityDrawer::CreateGravityImg(const CameraFrame::Ptr &frame, float scale) {
    // undistorted gray image
    cv::Mat undistImgColor, res;
    undistImgColor =
        CalibParamManager::ParIntri::UndistortImage(_intri->intri, frame->GetColorImage());

    // compute timestamp by reference IMU, we do not consider the readout time for RS cameras here
    double timeByBr = frame->GetTimestamp() + TO_DnToBr;
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    if (!so3Spline.TimeStampInRange(timeByBr)) {
        return undistImgColor;
    }

    auto SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
    auto SO3_Br0ToDn = (SO3_BrToBr0 * SE3_DnToBr.so3()).inverse();

    for (const auto &velCorr : _velCorrs[frame->GetId()]) {
        // map depth using alpha and beta
        const double depth = _intri->ActualDepth(velCorr->depth);

        // obtain the landmark
        Eigen::Vector2d lmInDnPlane = _intri->intri->ImgToCam(velCorr->MidPoint());
        Eigen::Vector3d lmInDn(lmInDnPlane(0) * depth, lmInDnPlane(1) * depth, depth);
        Eigen::Vector3d gravityInDn = SO3_Br0ToDn * GRAVITY;

        Eigen::Vector3d end = lmInDn + scale * gravityInDn;

        if (end(2) < 1E-3 || lmInDn(2) < 1E-3) {
            continue;
        }

        Eigen::Vector2d endPixel = _intri->intri->CamToImg({end(0) / end(2), end(1) / end(2)});

        Eigen::Vector2d feat = velCorr->MidPoint();

        DrawKeypointOnCVMat(undistImgColor, feat);
        DrawLineOnCVMat(undistImgColor, feat, endPixel);
    }
    return undistImgColor;
}
}  // namespace ns_ikalibr