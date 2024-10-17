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

#include "core/event_trace_sac.h"
#include "spdlog/spdlog.h"
#include "opengv/sac/Ransac.hpp"
#include "core/tracked_event_feature.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

FeatureTrackingTrace::FeatureTrackingTrace(double s_time,
                                           double e_time,
                                           Eigen::Vector3d x_parm,
                                           Eigen::Vector3d y_parm)
    : sTime(s_time),
      eTime(e_time),
      xParm(std::move(x_parm)),
      yParm(std::move(y_parm)) {}

FeatureTrackingTrace::Ptr FeatureTrackingTrace::Create(double s_time,
                                                       double e_time,
                                                       const Eigen::Vector3d& x_parm,
                                                       const Eigen::Vector3d& y_parm) {
    return std::make_shared<FeatureTrackingTrace>(s_time, e_time, x_parm, y_parm);
}

FeatureTrackingTrace::Ptr FeatureTrackingTrace::CreateFrom(
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

std::optional<Eigen::Vector2d> FeatureTrackingTrace::PositionAt(double t) const {
    if (t < sTime || t > eTime) {
        return std::nullopt;
    }
    double x = QuadraticCurveValueAt(t, xParm);
    double y = QuadraticCurveValueAt(t, yParm);
    return Eigen::Vector2d{x, y};
}

std::optional<Eigen::Vector2d> FeatureTrackingTrace::VelocityAt(double t) const {
    if (t < sTime || t > eTime) {
        return std::nullopt;
    }
    double vx = QuadraticCurveVelocityAt(t, xParm);
    double vy = QuadraticCurveVelocityAt(t, yParm);
    return Eigen::Vector2d{vx, vy};
}

std::vector<Eigen::Vector3d> FeatureTrackingTrace::DiscretePositions(double dt) const {
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

double FeatureTrackingTrace::QuadraticCurveValueAt(double x, const Eigen::Vector3d& params) {
    double a = params[0];
    double b = params[1];
    double c = params[2];
    return a * x * x + b * x + c;
}

double FeatureTrackingTrace::QuadraticCurveVelocityAt(double x, const Eigen::Vector3d& params) {
    double a = params[0];
    double b = params[1];
    return 2.0 * a * x + b;
}

Eigen::Vector3d FeatureTrackingTrace::FitQuadraticCurve(const std::vector<double>& x,
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

EventTrackingTraceSacProblem::EventTrackingTraceSacProblem(
    const std::vector<EventFeature::Ptr>& trackingAry, bool randomSeed)
    : opengv::sac::SampleConsensusProblem<model_t>(randomSeed),
      _trackingAry(trackingAry) {
    setUniformIndices(static_cast<int>(_trackingAry.size()));
}

std::pair<FeatureTrackingTrace::Ptr, std::vector<EventFeature::Ptr>>
EventTrackingTraceSacProblem::EventTrackingTraceSac(
    const std::vector<EventFeature::Ptr>& trackingAry, double thd) {
    opengv::sac::Ransac<EventTrackingTraceSacProblem> ransac;
    std::shared_ptr<EventTrackingTraceSacProblem> probPtr(
        new EventTrackingTraceSacProblem(trackingAry));
    ransac.sac_model_ = probPtr;
    ransac.threshold_ = thd;
    ransac.max_iterations_ = 20;
    if (ransac.computeModel()) {
        FeatureTrackingTrace trace;
        probPtr->optimizeModelCoefficients(ransac.inliers_, ransac.model_coefficients_, trace);
        std::vector<EventFeature::Ptr> goods(ransac.inliers_.size());
        for (int i = 0; i != static_cast<int>(ransac.inliers_.size()); ++i) {
            goods.at(i) = trackingAry.at(ransac.inliers_.at(i));
        }

        // spdlog::info(
        //     "total: {}, inliner: {}, rate: {:.3f}", trackingAry.size(), ransac.inliers_.size(),
        //     static_cast<double>(ransac.inliers_.size()) /
        //     static_cast<double>(trackingAry.size()));
        return {std::make_shared<FeatureTrackingTrace>(trace), goods};
    } else {
        spdlog::warn("compute event trace using RANSAC failed!!!");
        return {nullptr, {}};
    }
}

bool EventTrackingTraceSacProblem::computeModelCoefficients(const std::vector<int>& indices,
                                                            model_t& outModel) const {
    std::vector<EventFeature::Ptr> selected(indices.size());
    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        selected.at(i) = _trackingAry.at(indices.at(i));
    }
    FeatureTrackingTrace::Ptr trace = FeatureTrackingTrace::CreateFrom(selected);

    if (trace == nullptr) {
        return false;
    } else {
        outModel = *trace;
        return true;
    }
}

void EventTrackingTraceSacProblem::getSelectedDistancesToModel(const model_t& model,
                                                               const std::vector<int>& indices,
                                                               std::vector<double>& scores) const {
    scores.resize(indices.size());
    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        auto track = _trackingAry.at(indices.at(i));
        std::optional<Eigen::Vector2d> pos = model.PositionAt(track->timestamp);
        double score = std::numeric_limits<double>::max();
        if (pos != std::nullopt) {
            score = (*pos - track->pos).norm();
        }
        scores.at(i) = score;
    }
}

void EventTrackingTraceSacProblem::optimizeModelCoefficients(const std::vector<int>& inliers,
                                                             const model_t& model,
                                                             model_t& optimized_model) {
    computeModelCoefficients(inliers, optimized_model);
}

int EventTrackingTraceSacProblem::getSampleSize() const { return 3; }

}  // namespace ns_ikalibr
