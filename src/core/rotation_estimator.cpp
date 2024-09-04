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

#include "core/rotation_estimator.h"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

RotationEstimator::RotationEstimator()
    : _solveFlag(false),
      _sensorToSpline() {}

RotationEstimator::Ptr RotationEstimator::Create() { return std::make_shared<RotationEstimator>(); }

void RotationEstimator::Estimate(const So3SplineType &spline, const RotationSequence &rotSeq) {
    _solveFlag = false;

    std::vector<Eigen::Matrix4d> AMatSeq = OrganizeCoeffMatSeq(spline, rotSeq);

    if (AMatSeq.size() < 15) {
        return;
    }

    Eigen::MatrixXd AMat(AMatSeq.size() * 4, 4);
    for (int i = 0; i < static_cast<int>(AMatSeq.size()); ++i) {
        AMat.block<4, 4>(i * 4, 0) = AMatSeq.at(i);
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(AMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector4d cov = svd.singularValues();

    if (cov(2) > 0.25) {
        // get result
        Eigen::Matrix<double, 4, 1> x = svd.matrixV().col(3);
        Eigen::Quaterniond quat(x);
        Sophus::SO3d splineToSensor(quat);

        _solveFlag = true;
        _sensorToSpline = splineToSensor.inverse();
    }
}

void RotationEstimator::Estimate(const RotationEstimator::So3SplineType &spline,
                                 const std::vector<ns_ctraj::Posed> &poseSeq) {
    RotationSequence rotSeq(poseSeq.size());
    for (int i = 0; i < static_cast<int>(poseSeq.size()); ++i) {
        rotSeq.at(i).first = poseSeq.at(i).timeStamp;
        rotSeq.at(i).second = poseSeq.at(i).so3;
    }
    return Estimate(spline, rotSeq);
}

bool RotationEstimator::SolveStatus() const { return _solveFlag; }

const Sophus::SO3d &RotationEstimator::GetSO3SensorToSpline() const { return _sensorToSpline; }

std::vector<Eigen::Matrix4d> RotationEstimator::OrganizeCoeffMatSeq(
    const So3SplineType &spline, const RotationSequence &rotSeq) {
    std::vector<Eigen::Matrix4d> AMatSeq;

    for (int i = 1; i < static_cast<int>(rotSeq.size()); ++i) {
        int curIdx = i, lastIdx = i - 1;
        auto curTime = rotSeq[curIdx].first, lastTime = rotSeq[lastIdx].first;

        // check time stamp
        if (!spline.TimeStampInRange(curTime) || !spline.TimeStampInRange(lastTime)) {
            continue;
        }

        Eigen::Quaterniond curToLast, trajCurToLast;
        double factor;

        // sensor
        const auto &curRot = rotSeq[curIdx].second, lastRot = rotSeq[lastIdx].second;
        curToLast = (lastRot.inverse() * curRot).unit_quaternion();

        // spline
        auto curToRef = spline.Evaluate(curTime), lastToRef = spline.Evaluate(lastTime);
        trajCurToLast = (lastToRef.inverse() * curToRef).unit_quaternion();

        {
            Eigen::AngleAxisd sensorAngleAxis(curToLast.toRotationMatrix());
            Eigen::AngleAxisd trajAngleAxis(trajCurToLast.toRotationMatrix());

            // compute weight factor
            constexpr static double RAD_TO_DEG = 180 / M_PI;
            double deltaAngle =
                RAD_TO_DEG * std::fabs(sensorAngleAxis.angle() - trajAngleAxis.angle());
            factor = deltaAngle > 1.0 ? 1.0 / deltaAngle : 1.0;
        }

        Eigen::Matrix4d lqMat = LeftQuatMatrix(curToLast);
        Eigen::Matrix4d rqMat = RightQuatMatrix(trajCurToLast);
        AMatSeq.emplace_back(factor * (lqMat - rqMat));
    }
    return AMatSeq;
}
}  // namespace ns_ikalibr