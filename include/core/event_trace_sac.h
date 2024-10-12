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

#ifndef IKALIBR_EVENT_TRACE_SAC_H
#define IKALIBR_EVENT_TRACE_SAC_H

#include "core/haste_data_io.h"
#include "opengv/sac/SampleConsensusProblem.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct EventTrackingTrace {
public:
    using Ptr = std::shared_ptr<EventTrackingTrace>;

public:
    double sTime{}, eTime{};
    Eigen::Vector3d xParm, yParm;

public:
    EventTrackingTrace(double s_time,
                       double e_time,
                       Eigen::Vector3d x_parm,
                       Eigen::Vector3d y_parm);

    EventTrackingTrace() = default;

    static Ptr Create(double s_time,
                      double e_time,
                      const Eigen::Vector3d &x_parm,
                      const Eigen::Vector3d &y_parm);

    static Ptr CreateFrom(const std::vector<HASTEFeature::Ptr> &trackingAry);

    [[nodiscard]] std::optional<Eigen::Vector2d> PositionAt(double t) const;

    [[nodiscard]] std::vector<Eigen::Vector3d> DiscretePositions(double dt = 0.005) const;

protected:
    static double QuadraticCurveValueAt(double x, const Eigen::Vector3d &params);

    static Eigen::Vector3d FitQuadraticCurve(const std::vector<double> &x,
                                             const std::vector<double> &y);
};

class EventTrackingTraceSacProblem
    : public opengv::sac::SampleConsensusProblem<EventTrackingTrace> {
public:
    /** The model we are trying to fit */
    typedef EventTrackingTrace model_t;

public:
    /**
     * \brief Constructor.
     */
    explicit EventTrackingTraceSacProblem(const std::vector<HASTEFeature::Ptr> &trackingAry,
                                          bool randomSeed = true);

    static std::pair<EventTrackingTrace::Ptr, std::vector<HASTEFeature::Ptr>> EventTrackingTraceSac(
        const std::vector<HASTEFeature::Ptr> &trackingAry, double thd);

    /**
     * Destructor.
     */
    ~EventTrackingTraceSacProblem() override = default;

    /**
     * \brief See parent-class.
     */
    bool computeModelCoefficients(const std::vector<int> &indices,
                                  model_t &outModel) const override;

    /**
     * \brief See parent-class.
     */
    void getSelectedDistancesToModel(const model_t &model,
                                     const std::vector<int> &indices,
                                     std::vector<double> &scores) const override;

    /**
     * \brief See parent-class.
     */
    void optimizeModelCoefficients(const std::vector<int> &inliers,
                                   const model_t &model,
                                   model_t &optimized_model) override;

    /**
     * \brief See parent-class.
     */
    [[nodiscard]] int getSampleSize() const override;

protected:
    static double QuadraticCurve(double x, double a, double b, double c);

    static double DistanceSquared(double x0, double y0, double x1, double a, double b, double c);

    static double Gradient(double x0, double y0, double x1, double a, double b, double c);

protected:
    /** The adapter holding all input data */
    std::vector<HASTEFeature::Ptr> _trackingAry;
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_EVENT_TRACE_SAC_H
