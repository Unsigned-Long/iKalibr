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

#include "viewer/visual_lin_vel_drawer.h"
#include "calib/calib_param_manager.h"
#include "opencv2/imgproc.hpp"
#include "sensor/rgbd.h"
#include "factor/data_correspondence.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
// ------------------
// VisualLinVelDrawer
// ------------------
VisualLinVelDrawer::VisualLinVelDrawer(const std::string &topic,
                                       ns_veta::Veta::Ptr veta,
                                       SplineBundleType::Ptr splines,
                                       const CalibParamManager::Ptr &parMagr)
    : _veta(std::move(veta)),
      _splines(std::move(splines)) {
    _intri = parMagr->INTRI.Camera.at(topic);

    SE3_CmToBr = parMagr->EXTRI.SE3_CmToBr(topic);
    TO_CmToBr = parMagr->TEMPORAL.TO_CmToBr.at(topic);
}

VisualLinVelDrawer::Ptr VisualLinVelDrawer::Create(const std::string &topic,
                                                   const ns_veta::Veta::Ptr &veta,
                                                   const SplineBundleType::Ptr &splines,
                                                   const CalibParamManager::Ptr &parMagr) {
    return std::make_shared<VisualLinVelDrawer>(topic, veta, splines, parMagr);
}

cv::Mat VisualLinVelDrawer::CreateLinVelImg(const CameraFrame::Ptr &frame, float scale) {
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
    auto SE3_CmToBr0 = SE3_BrToBr0 * SE3_CmToBr;

    Eigen::Vector3d ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * so3Spline.VelocityBody(timeByBr);
    Eigen::Vector3d LIN_VEL_BrToBr0InBr0 = posSpline.Evaluate<1>(timeByBr);
    Eigen::Vector3d LIN_VEL_CmToBr0InBr0 =
        LIN_VEL_BrToBr0InBr0 -
        Sophus::SO3d::hat(SO3_BrToBr0 * SE3_CmToBr.translation()) * ANG_VEL_BrToBr0InBr0;

    for (const auto &[lmId, lm] : _veta->structure) {
        auto iter = lm.obs.find(frame->GetId());
        if (iter == lm.obs.cend()) {
            continue;
        }

        Eigen::Vector3d lmInCm = SE3_CmToBr0.inverse() * lm.X;
        Eigen::Vector3d val1 =
            Sophus::SO3d::hat(SE3_CmToBr0.so3() * lmInCm) * ANG_VEL_BrToBr0InBr0 -
            LIN_VEL_CmToBr0InBr0;
        Eigen::Vector3d val2 = SE3_CmToBr0.so3().inverse() * val1;

        Eigen::Vector3d end = lmInCm + scale * val2;

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

// ----------------------
// RGBDVisualLinVelDrawer
// ----------------------

VisualOpticalFlowLinVelDrawer::VisualOpticalFlowLinVelDrawer(
    const std::vector<OpticalFlowCorr::Ptr> &corrs,
    SplineBundleType::Ptr splines,
    ns_veta::PinholeIntrinsic::Ptr intri,
    const Sophus::SE3d &SE3_SenToBr,
    const double &TO_SenToBr,
    const double &RS_READOUT)
    : _splines(std::move(splines)),
      _intri(std::move(intri)),
      SE3_SenToBr(SE3_SenToBr),
      TO_SenToBr(TO_SenToBr),
      RS_READOUT(RS_READOUT) {
    for (const auto &corr : corrs) {
        if (corr->depth > 1E-3 /* 1 mm */) {
            // a valid depth
            this->_velCorrs[corr->frame->GetId()].emplace_back(corr);
        }
    }
}

VisualOpticalFlowLinVelDrawer::Ptr VisualOpticalFlowLinVelDrawer::CreateDrawerForRGBDs(
    const std::string &topic,
    const std::vector<OpticalFlowCorr::Ptr> &corrs,
    const SplineBundleType::Ptr &splines,
    const CalibParamManager::Ptr &parMagr) {
    return std::make_shared<VisualOpticalFlowLinVelDrawer>(
        corrs, splines, parMagr->INTRI.RGBD.at(topic)->intri, parMagr->EXTRI.SE3_DnToBr(topic),
        parMagr->TEMPORAL.TO_DnToBr.at(topic), parMagr->TEMPORAL.RS_READOUT.at(topic));
}

VisualOpticalFlowLinVelDrawer::Ptr VisualOpticalFlowLinVelDrawer::CreateDrawerForVelCameras(
    const std::string &topic,
    const std::vector<OpticalFlowCorrPtr> &corrs,
    const SplineBundleType::Ptr &splines,
    const CalibParamManagerPtr &parMagr) {
    return std::make_shared<VisualOpticalFlowLinVelDrawer>(
        corrs, splines, parMagr->INTRI.Camera.at(topic), parMagr->EXTRI.SE3_CmToBr(topic),
        parMagr->TEMPORAL.TO_CmToBr.at(topic), parMagr->TEMPORAL.RS_READOUT.at(topic));
}

cv::Mat VisualOpticalFlowLinVelDrawer::CreateLinVelImg(
    const CameraFrame::Ptr &frame, const TimeDeriv::ScaleSplineType &scaleSplineType, float scale) {
    // undistorted gray image
    cv::Mat undistImgColor, res;
    undistImgColor = CalibParamManager::ParIntri::UndistortImage(_intri, frame->GetColorImage());

    // compute timestamp by reference IMU, we do not consider the readout time for RS cameras here
    double timeByBr = frame->GetTimestamp() + TO_SenToBr;
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    if (!so3Spline.TimeStampInRange(timeByBr) || !scaleSpline.TimeStampInRange(timeByBr)) {
        return undistImgColor;
    }

    auto SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
    auto SO3_SenToBr0 = SO3_BrToBr0 * SE3_SenToBr.so3();
    Sophus::SO3d SO3_BrToSen = SE3_SenToBr.so3().inverse();

    Eigen::Vector3d ANG_VEL_BrToBr0InBr = so3Spline.VelocityBody(timeByBr);
    Eigen::Vector3d ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * ANG_VEL_BrToBr0InBr;
    Eigen::Vector3d ANG_VEL_SenToBr0InSen = SO3_BrToSen * ANG_VEL_BrToBr0InBr;

    Eigen::Vector3d LIN_VEL_BrToBr0InBr0;
    switch (scaleSplineType) {
        case TimeDeriv::LIN_ACCE_SPLINE:
            // this would not happen
            return undistImgColor;
        case TimeDeriv::LIN_VEL_SPLINE:
            LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<0>(timeByBr);
            break;
        case TimeDeriv::LIN_POS_SPLINE:
            LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<1>(timeByBr);
            break;
    }

    Eigen::Vector3d LIN_VEL_SenToBr0InBr0 =
        LIN_VEL_BrToBr0InBr0 -
        Sophus::SO3d::hat(SO3_BrToBr0 * SE3_SenToBr.translation()) * ANG_VEL_BrToBr0InBr0;
    Eigen::Vector3d LIN_VEL_SenToBr0InSen =
        SO3_BrToSen * SO3_BrToBr0.inverse() * LIN_VEL_SenToBr0InBr0;

    const double FX = _intri->FocalX(), FY = _intri->FocalY();
    const double CX = _intri->PrincipalPoint()(0), CY = _intri->PrincipalPoint()(1);

    for (const auto &velCorr : _velCorrs[frame->GetId()]) {
        // map depth using alpha and beta
        const double depth = velCorr->depth;

        // draw linear velocity of the landmark with respect to the rgbd camera
        {  // obtain the landmark
            Eigen::Vector2d lmInCamPlane = _intri->ImgToCam(velCorr->MidPoint());
            Eigen::Vector3d lmInCam(lmInCamPlane(0) * depth, lmInCamPlane(1) * depth, depth);

            Eigen::Vector3d val1 =
                Sophus::SO3d::hat(SO3_SenToBr0 * lmInCam) * ANG_VEL_BrToBr0InBr0 -
                LIN_VEL_SenToBr0InBr0;
            Eigen::Vector3d val2 = SO3_SenToBr0.inverse() * val1;

            Eigen::Vector3d end = lmInCam + scale * val2;

            if (end(2) < 1E-3 || lmInCam(2) < 1E-3) {
                continue;
            }

            Eigen::Vector2d endPixel = _intri->CamToImg({end(0) / end(2), end(1) / end(2)});

            Eigen::Vector2d feat = velCorr->MidPoint();

            DrawKeypointOnCVMat(undistImgColor, feat);
            DrawLineOnCVMat(undistImgColor, feat, endPixel);
        }

        // draw computed pixel velocity
        {
            Eigen::Vector2d endPixel =
                velCorr->MidPoint() - scale * 0.5 * velCorr->MidPointVel(RS_READOUT);
            DrawLineOnCVMat(undistImgColor, velCorr->MidPoint(), endPixel, cv::Scalar(255, 0, 0));
        }

        // draw estimated pixel velocity
        {
            Eigen::Matrix<double, 2, 3> subAMat, subBMat;
            OpticalFlowCorr::SubMats<double>(&FX, &FY, &CX, &CY, velCorr->MidPoint(), &subAMat,
                                             &subBMat);
            Eigen::Vector2d pred =
                1.0 / depth * subAMat * LIN_VEL_SenToBr0InSen + subBMat * ANG_VEL_SenToBr0InSen;

            Eigen::Vector2d endPixel = velCorr->MidPoint() - scale * 0.5 * pred;
            DrawLineOnCVMat(undistImgColor, velCorr->MidPoint(), endPixel, cv::Scalar(0, 0, 255));
        }
    }
    return undistImgColor;
}

}  // namespace ns_ikalibr