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

#include "core/visual_velocity_estimator.h"
#include "opencv2/imgproc.hpp"
#include "sensor/camera.h"
#include "calib/calib_param_manager.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

VisualVelocityEstimator::VisualVelocityEstimator(
    const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>>& dynamics,
    ns_veta::PinholeIntrinsic::Ptr intri)
    : _dynamics(dynamics),
      _intri(std::move(intri)) {}

VisualVelocityEstimator::Ptr VisualVelocityEstimator::Create(
    const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, double>>& dynamics,
    const ns_veta::PinholeIntrinsic::Ptr& intri) {
    return std::make_shared<VisualVelocityEstimator>(dynamics, intri);
}

std::optional<Eigen::Vector3d> VisualVelocityEstimator::Estimate(
    double timeByBr,
    const VisualVelocityEstimator::So3SplineType& spline,
    const Sophus::SO3d& SO3_DnToBr) const {
    if (_dynamics.size() < 2) {
        return {};
    }
    const double fx = _intri->FocalX(), fy = _intri->FocalY();
    const double cx = _intri->PrincipalPoint()(0), cy = _intri->PrincipalPoint()(1);

    int size = static_cast<int>(_dynamics.size());
    Eigen::MatrixXd AMat = Eigen::MatrixXd ::Zero(2 * size, 3);
    Eigen::MatrixXd BMat = Eigen::MatrixXd ::Zero(2 * size, 1);
    Eigen::MatrixXd VMat = Eigen::MatrixXd ::Zero(2 * size, 1);
    for (int i = 0; i < size; ++i) {
        const auto& [pixel, vel, depth] = _dynamics.at(i);
        const double up = pixel(0) - cx, vp = pixel(1) - cy;

        Eigen::Matrix<double, 2, 3> subAMat = Eigen::Matrix<double, 2, 3>::Zero();
        subAMat(0, 0) = -fx;
        subAMat(1, 1) = -fy;
        subAMat(0, 2) = up;
        subAMat(1, 2) = vp;
        AMat.block<2, 3>(i * 2, 0) = 1 / depth * subAMat;

        Eigen::Matrix<double, 2, 3> subBMat = Eigen::Matrix<double, 2, 3>::Zero();
        subBMat(0, 0) = up * vp / fy;
        subBMat(0, 1) = -fx - up * up / fx;
        subBMat(0, 2) = fx * vp / fy;
        subBMat(1, 0) = fy + vp * vp / fy;
        subBMat(1, 1) = -up * vp / fx;
        subBMat(1, 2) = -fy * up / fx;

        Eigen::Vector3d ANG_VEL_BrToWInBr = spline.VelocityBody(timeByBr);
        Eigen::Vector3d ANG_VEL_DnToWInDn = SO3_DnToBr.inverse() * ANG_VEL_BrToWInBr;

        BMat.block<2, 1>(i * 2, 0) = subBMat * ANG_VEL_DnToWInDn;

        VMat.block<2, 1>(i * 2, 0) = vel;
    }
    Eigen::MatrixXd lVec = VMat - BMat;
    Eigen::Matrix3d HMat = (AMat.transpose() * AMat).inverse();

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(HMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d cov = svd.singularValues();
    if (cov(2) < 0.01) {
        return {};
    }
    Eigen::Vector3d LIN_VEL_DnToWInDn = HMat * AMat.transpose() * lVec;

    return {LIN_VEL_DnToWInDn};
}

cv::Mat VisualVelocityEstimator::DrawVisualVelocityMat(const Eigen::Vector3d& LIN_VEL_DnToWInDn,
                                                       const CameraFrame::Ptr& frame,
                                                       double factor) {
    cv::Mat img = CalibParamManager::ParIntri::UndistortImage(_intri, frame->GetColorImage());
    for (const auto& [feat, vel, depth] : _dynamics) {
        Eigen::Vector2d lmInPlane = _intri->ImgToCam(feat);
        Eigen::Vector3d lmInCm(lmInPlane(0) * depth, lmInPlane(1) * depth, depth);

        // draw negative velocity
        Eigen::Vector3d end = lmInCm - factor * LIN_VEL_DnToWInDn;

        if (end(2) < 1E-3 || lmInCm(2) < 1E-3) {
            continue;
        }

        Eigen::Vector2d endPixel = _intri->CamToImg({end(0) / end(2), end(1) / end(2)});

        // square
        cv::drawMarker(img, cv::Point2d(feat(0), feat(1)), cv::Scalar(0, 255, 0),
                       cv::MarkerTypes::MARKER_SQUARE, 10, 1);

        // key point
        cv::drawMarker(img, cv::Point2d(feat(0), feat(1)), cv::Scalar(0, 255, 0),
                       cv::MarkerTypes::MARKER_SQUARE, 2, 2);

        // tail
        cv::line(img, cv::Point2d(feat(0), feat(1)), cv::Point2d(endPixel(0), endPixel(1)),
                 cv::Scalar(0, 255, 0), 1, cv::LINE_AA);
    }
    return img;
}

}  // namespace ns_ikalibr