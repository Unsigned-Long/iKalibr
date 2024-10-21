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

#include "factor/data_correspondence.h"
#include "sensor/camera.h"
#include "core/tracked_event_feature.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
PointToSurfelCorr::PointToSurfelCorr(double timestamp,
                                     Eigen::Vector3d pInScan,
                                     double weight,
                                     Eigen::Vector4d surfelInW)
    : timestamp(timestamp),
      pInScan(std::move(pInScan)),
      weight(weight),
      surfelInW(std::move(surfelInW)) {}

PointToSurfelCorr::Ptr PointToSurfelCorr::Create(double timestamp,
                                                 const Eigen::Vector3d& pInScan,
                                                 double weight,
                                                 const Eigen::Vector4d& surfelInW) {
    return std::make_shared<PointToSurfelCorr>(timestamp, pInScan, weight, surfelInW);
}

VisualReProjCorr::VisualReProjCorr(double ti,
                                   double tj,
                                   Eigen::Vector2d fi,
                                   Eigen::Vector2d fj,
                                   double li,
                                   double lj,
                                   double weight)
    : ti(ti),
      tj(tj),
      fi(std::move(fi)),
      fj(std::move(fj)),
      li(li),
      lj(lj),
      weight(weight) {}

VisualReProjCorr::Ptr VisualReProjCorr::Create(double ti,
                                               double tj,
                                               const Eigen::Vector2d& fi,
                                               const Eigen::Vector2d& fj,
                                               double li,
                                               double lj,
                                               double weight) {
    return std::make_shared<VisualReProjCorr>(ti, tj, fi, fj, li, lj, weight);
}

VisualReProjCorr::VisualReProjCorr() = default;

VisualReProjCorrSeq::VisualReProjCorrSeq() = default;

OpticalFlowCorr::OpticalFlowCorr(const std::array<double, 3>& timeAry,
                                 const std::array<double, 3>& xTraceAry,
                                 const std::array<double, 3>& yTraceAry,
                                 double depth,
                                 const CameraFramePtr& frame,
                                 double rsExpFactor)
    : timeAry(timeAry),
      xTraceAry(xTraceAry),
      yTraceAry(yTraceAry),
      rdFactorAry(),
      depth(depth),
      invDepth(depth > 1E-3 ? 1.0 / depth : -1.0),
      frame(frame),
      withDepthObservability(false) {
    int imgHeight = frame->GetImage().rows;
    for (int i = 0; i < 3; ++i) {
        rdFactorAry[i] = yTraceAry[i] / (double)imgHeight - rsExpFactor;
    }
}

OpticalFlowCorr::OpticalFlowCorr(const std::array<double, 3>& timeAry,
                                 const std::array<double, 3>& xTraceAry,
                                 const std::array<double, 3>& yTraceAry,
                                 double depth)
    : timeAry(timeAry),
      xTraceAry(xTraceAry),
      yTraceAry(yTraceAry),
      rdFactorAry(),
      depth(depth),
      invDepth(depth > 1E-3 ? 1.0 / depth : -1.0),
      frame(nullptr),
      withDepthObservability(false) {
    for (int i = 0; i < 3; ++i) {
        rdFactorAry[i] = 0.0;
    }
}

OpticalFlowCorr::Ptr OpticalFlowCorr::Create(const std::array<double, 3>& timeAry,
                                             const std::array<double, 3>& xDynamicAry,
                                             const std::array<double, 3>& yDynamicAry,
                                             double depth,
                                             const CameraFramePtr& frame,
                                             double rsExpFactor) {
    return std::make_shared<OpticalFlowCorr>(timeAry, xDynamicAry, yDynamicAry, depth, frame,
                                             rsExpFactor);
}

OpticalFlowCorr::Ptr OpticalFlowCorr::Create(const std::array<double, 3>& timeAry,
                                             const std::array<double, 3>& xDynamicAry,
                                             const std::array<double, 3>& yDynamicAry,
                                             double depth) {
    return std::make_shared<OpticalFlowCorr>(timeAry, xDynamicAry, yDynamicAry, depth);
}

Eigen::Vector2d OpticalFlowCorr::FirPoint() const { return {xTraceAry.at(FIR), yTraceAry.at(FIR)}; }

Eigen::Vector2d OpticalFlowCorr::MidPoint() const { return {xTraceAry.at(MID), yTraceAry.at(MID)}; }

Eigen::Vector2d OpticalFlowCorr::LastPoint() const {
    return {xTraceAry.at(LAST), yTraceAry.at(LAST)};
}

double OpticalFlowCorr::MidReadoutFactor() const { return rdFactorAry[MID]; }

FeatureTrackingCurve::FeatureTrackingCurve(double s_time,
                                           double e_time,
                                           Eigen::Vector3d x_parm,
                                           Eigen::Vector3d y_parm)
    : sTime(s_time),
      eTime(e_time),
      xParm(std::move(x_parm)),
      yParm(std::move(y_parm)) {}

FeatureTrackingCurve::Ptr FeatureTrackingCurve::Create(double s_time,
                                                       double e_time,
                                                       const Eigen::Vector3d& x_parm,
                                                       const Eigen::Vector3d& y_parm) {
    return std::make_shared<FeatureTrackingCurve>(s_time, e_time, x_parm, y_parm);
}

FeatureTrackingCurve::Ptr FeatureTrackingCurve::CreateFrom(
    const std::vector<EventFeature::Ptr>& trackingAry) {
    if (trackingAry.size() < 3) {
        return nullptr;
    }
    const auto size = trackingAry.size();
    std::vector<double> t(size), x(size), y(size);
    double tMin = std::numeric_limits<double>::max();
    double tMax = std::numeric_limits<double>::min();
    for (int i = 0; i != static_cast<int>(size); ++i) {
        const auto& track = trackingAry[i];
        t.at(i) = track->timestamp;
        x.at(i) = track->pos(0);
        y.at(i) = track->pos(1);
        if (tMin > track->timestamp) {
            tMin = track->timestamp;
        }
        if (tMax < track->timestamp) {
            tMax = track->timestamp;
        }
    }
    Eigen::Vector3d xParm = FitQuadraticCurve(t, x);
    Eigen::Vector3d yParm = FitQuadraticCurve(t, y);

    return Create(tMin, tMax, xParm, yParm);
}

std::optional<Eigen::Vector2d> FeatureTrackingCurve::PositionAt(double t) const {
    if (t < sTime || t > eTime) {
        return std::nullopt;
    }
    double x = QuadraticCurveValueAt(t, xParm);
    double y = QuadraticCurveValueAt(t, yParm);
    return Eigen::Vector2d{x, y};
}

std::optional<Eigen::Vector2d> FeatureTrackingCurve::VelocityAt(double t) const {
    if (t < sTime || t > eTime) {
        return std::nullopt;
    }
    double vx = QuadraticCurveVelocityAt(t, xParm);
    double vy = QuadraticCurveVelocityAt(t, yParm);
    return Eigen::Vector2d{vx, vy};
}

std::vector<Eigen::Vector3d> FeatureTrackingCurve::DiscretePositions(double dt) const {
    std::vector<Eigen::Vector3d> positions;
    for (double t = this->sTime; t < this->eTime;) {
        std::optional<Eigen::Vector2d> pos = this->PositionAt(t);
        if (pos == std::nullopt) {
            continue;
        }
        positions.emplace_back(t, (*pos)(0), (*pos)(1));
        t += dt;
    }
    return positions;
}

bool FeatureTrackingCurve::IsTimeInRange(double time) const {
    return time >= sTime && time <= eTime;
}

Eigen::Vector3d FeatureTrackingCurve::FitQuadraticCurve(const std::vector<double>& x,
                                                        const std::vector<double>& y) {
    int n = static_cast<int>(x.size());
    Eigen::MatrixXd A(n, 3);
    Eigen::VectorXd B(n);

    for (int i = 0; i < n; ++i) {
        A(i, 0) = x[i] * x[i];  // x^2
        A(i, 1) = x[i];         // x
        A(i, 2) = 1.0;
        B(i) = y[i];
    }

    // A * p = B，p = (a, b, c)
    Eigen::Vector3d p = A.colPivHouseholderQr().solve(B);
    return p;  // a, b, c
}

OpticalFlowCurveCorr::Ptr OpticalFlowCurveCorr::Create(double midTime,
                                                         double midDepth,
                                                         double reprojTimePadding,
                                                         const FeatureTrackingCurve::Ptr& trace,
                                                         double weight) {
    if (reprojTimePadding < 1E-3) {
        return nullptr;
    }
    double firTime = midTime - reprojTimePadding;
    double lastTime = midTime + reprojTimePadding;
    if (!trace->IsTimeInRange(firTime) || !trace->IsTimeInRange(lastTime)) {
        return nullptr;
    }
    return std::make_shared<OpticalFlowCurveCorr>(midTime, midDepth, 1.0 / midDepth, firTime,
                                                   lastTime, trace, weight);
}

OpticalFlowCurveCorr::OpticalFlowCurveCorr(double mid_time,
                                             double mid_depth,
                                             double mid_inv_depth,
                                             double fir_time,
                                             double last_time,
                                             const FeatureTrackingCurve::Ptr& trace,
                                             double weight)
    : midTime(mid_time),
      midDepth(mid_depth),
      midInvDepth(mid_inv_depth),
      firTime(fir_time),
      lastTime(last_time),
      trace(trace),
      weight(weight) {}
}  // namespace ns_ikalibr
