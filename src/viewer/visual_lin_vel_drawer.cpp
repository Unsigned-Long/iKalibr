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
#include "sensor/camera.h"

_3_

namespace ns_ikalibr {

    VisualLinVelDrawer::VisualLinVelDrawer(std::string topic, ns_veta::Veta::Ptr veta,
                                           SplineBundleType::Ptr splines, CalibParamManager::Ptr parMagr)
            : _topic(std::move(topic)), _veta(std::move(veta)), _splines(std::move(splines)),
              _parMagr(std::move(parMagr)) {

        _intri = _parMagr->INTRI.Camera.at(_topic);

        SE3_CmToBr = _parMagr->EXTRI.SE3_CmToBr(_topic);
        TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(_topic);
    }

    VisualLinVelDrawer::Ptr
    VisualLinVelDrawer::Create(const std::string &topic, const ns_veta::Veta::Ptr &veta,
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
        if (!so3Spline.TimeStampInRange(timeByBr) || !posSpline.TimeStampInRange(timeByBr)) { return undistImgColor; }

        auto SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
        Eigen::Vector3d POS_BrInBr0 = posSpline.Evaluate(timeByBr);
        Sophus::SE3d SE3_BrToBr0(SO3_BrToBr0, POS_BrInBr0);
        auto SE3_CmToBr0 = SE3_BrToBr0 * SE3_CmToBr;

        Eigen::Vector3d ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * so3Spline.VelocityBody(timeByBr);
        Eigen::Vector3d LIN_VEL_BrToBr0InBr0 = posSpline.Evaluate<1>(timeByBr);
        Eigen::Vector3d LIN_VEL_CmToBr0InBr0 =
                LIN_VEL_BrToBr0InBr0 - Sophus::SO3d::hat(SO3_BrToBr0 * SE3_CmToBr.translation()) * ANG_VEL_BrToBr0InBr0;

        for (const auto &[lmId, lm]: _veta->structure) {
            auto iter = lm.obs.find(frame->GetId());
            if (iter == lm.obs.cend()) { continue; }

            Eigen::Vector3d lmInCm = SE3_CmToBr0.inverse() * lm.X;
            Eigen::Vector3d val1 = Sophus::SO3d::hat(SE3_CmToBr0.so3() * lmInCm) * ANG_VEL_BrToBr0InBr0
                                   - LIN_VEL_CmToBr0InBr0;
            Eigen::Vector3d val2 = SE3_CmToBr0.so3().inverse() * val1;

            Eigen::Vector3d end = lmInCm + scale * val2;
            Eigen::Vector2d endPixel = _intri->CamToImg({end(0) / end(2), end(1) / end(2)});

            // we do not use extracted raw feature here to keep better consistency
            Eigen::Vector2d feat = _intri->CamToImg({lmInCm(0) / lmInCm(2), lmInCm(1) / lmInCm(2)});

            // square
            cv::drawMarker(
                    undistImgColor, cv::Point2d(feat(0), feat(1)), cv::Scalar(0, 255, 0),
                    cv::MarkerTypes::MARKER_SQUARE, 10, 1
            );

            // key point
            cv::drawMarker(
                    undistImgColor, cv::Point2d(feat(0), feat(1)), cv::Scalar(0, 255, 0),
                    cv::MarkerTypes::MARKER_SQUARE, 2, 2);

            // tail
            cv::line(
                    undistImgColor, cv::Point2d(feat(0), feat(1)),
                    cv::Point2d(endPixel(0), endPixel(1)), cv::Scalar(0, 255, 0), 1, cv::LINE_AA
            );
        }
        return undistImgColor;
    }
}