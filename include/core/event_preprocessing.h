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

#ifndef EVENT_PREPROCESSING_H
#define EVENT_PREPROCESSING_H

#include "util/utils.h"
#include "sensor/event.h"
#include "opencv2/imgproc.hpp"
#include "core/visual_distortion.h"
#include "opengv/sac/SampleConsensusProblem.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_veta {
struct PinholeIntrinsic;
using PinholeIntrinsicPtr = std::shared_ptr<PinholeIntrinsic>;
}  // namespace ns_veta

namespace ns_ikalibr {

struct VisualUndistortionMap;
using VisualUndistortionMapPtr = std::shared_ptr<VisualUndistortionMap>;

class ActiveEventSurface {
public:
    using Ptr = std::shared_ptr<ActiveEventSurface>;

private:
    const double FILTER_THD;
    ns_veta::PinholeIntrinsicPtr _intri;
    VisualUndistortionMapPtr _undistoMap;

    Eigen::MatrixXd _sae[2];        // save sae
    Eigen::MatrixXd _saeLatest[2];  // save previous sae
    double _timeLatest;

    cv::Mat _eventImgMat;

public:
    explicit ActiveEventSurface(const ns_veta::PinholeIntrinsicPtr &intri, double filterThd = 0.01);

    static Ptr Create(const ns_veta::PinholeIntrinsicPtr &intri, double filterThd = 0.01);

    void GrabEvent(const Event::Ptr &event, bool drawEventMat = false);

    void GrabEvent(const EventArray::Ptr &events, bool drawEventMat = false);

    [[nodiscard]] cv::Mat GetEventImgMat(bool resetMat, bool undistoMat = false);

    cv::Mat TimeSurface(bool ignorePolarity = false,
                        bool undistoMat = false,
                        int medianBlurKernelSize = 0,
                        double decaySec = 0.02);

    std::pair<cv::Mat, cv::Mat> RawTimeSurface(bool ignorePolarity = false,
                                               bool undistoMat = false);

    [[nodiscard]] double GetTimeLatest() const;
};

struct NormFlow;
using NormFlowPtr = std::shared_ptr<NormFlow>;

class EventNormFlow {
public:
    using Ptr = std::shared_ptr<EventNormFlow>;

    struct NormFlowPack {
        using Ptr = std::shared_ptr<NormFlowPack>;
        using NormFlowContainer = std::map<NormFlowPtr, std::vector<std::tuple<int, int, double>>>;

    public:
        // norm flow, [x, y, timestamp]
        std::map<NormFlowPtr, std::vector<std::tuple<int, int, double>>> nfs;
        cv::Mat rawTimeSurfaceMap;  // ex, ey, et
        cv::Mat polarityMap;        // the polarity map
        double timestamp;
        // for visualization
        cv::Mat nfSeedsImg;
        cv::Mat nfsImg;
        cv::Mat tsImg;

    public:
        std::list<Event::Ptr> ActiveEvents(double dt = 0.02) const;

        std::list<Event::Ptr> NormFlowInlierEvents() const;

        cv::Mat Visualization(double dt = 0.02) const;

        cv::Mat InliersOccupyMat() const;

        cv::Mat NormFlowInlierEventMat() const;

        cv::Mat AccumulativeEventMat(double dt = 0.02) const;

        std::list<NormFlowPtr> TemporallySortedNormFlows() const;

        std::list<NormFlowPtr> NormFlows() const;

        int Rows() const;

        int Cols() const;
    };

private:
    ActiveEventSurface::Ptr _sea;

public:
    explicit EventNormFlow(const ActiveEventSurface::Ptr &sea)
        : _sea(sea) {}

    NormFlowPack::Ptr ExtractNormFlows(double decaySec = 0.02,
                                       int winSize = 2,
                                       int neighborDist = 2,
                                       double goodRatioThd = 0.9,
                                       bool undistorted = true,
                                       double timeDistEventToPlaneThd = 2E-3,
                                       int ransacMaxIter = 3) const;

protected:
    static std::vector<std::tuple<double, double, double>> Centralization(
        const std::vector<std::tuple<int, int, double>> &inRangeData);
};

class EventLocalPlaneSacProblem : public opengv::sac::SampleConsensusProblem<Eigen::Vector3d> {
public:
    typedef Eigen::Vector3d model_t;

public:
    explicit EventLocalPlaneSacProblem(const std::vector<std::tuple<double, double, double>> &data,
                                       bool randomSeed = true)
        : opengv::sac::SampleConsensusProblem<model_t>(randomSeed),
          _data(data) {
        setUniformIndices(static_cast<int>(_data.size()));
    };

    ~EventLocalPlaneSacProblem() override = default;

    bool computeModelCoefficients(const std::vector<int> &indices,
                                  model_t &outModel) const override;

    void getSelectedDistancesToModel(const model_t &model,
                                     const std::vector<int> &indices,
                                     std::vector<double> &scores) const override;

    void optimizeModelCoefficients(const std::vector<int> &inliers,
                                   const model_t &model,
                                   model_t &optimized_model) override;

    [[nodiscard]] int getSampleSize() const override;

    static double PointToPlaneDistance(double x, double y, double t, double A, double B, double C);

protected:
    /** The adapter holding all input data */
    const std::vector<std::tuple<double, double, double>> &_data;
};

struct EventLine {
    using Ptr = std::shared_ptr<EventLine>;

public:
    double timestamp;
    double activity;
    // norm (n: cos, sin) and the distance from the line to the origin (rho)
    Eigen::Vector3d param;

public:
    EventLine(const NormFlowPtr &nf);

    static Ptr Create(const NormFlowPtr &nf);

    double PointToLine(const Eigen::Vector2d &p) const;

    double DirectionDifferenceCos(const Eigen::Vector2d &nfDir) const;

    void Normalize();
};

class EventLineTracking {
public:
    using Ptr = std::shared_ptr<EventLineTracking>;

    struct LineParamUpdate {
    public:
        using Ptr = std::shared_ptr<LineParamUpdate>;

        double x;
        double y;
        double xx;
        double yy;
        double xy;

        LineParamUpdate(double x, double y);

        static Ptr Create(double x, double y);

        void Update(double x, double y, double delta);

        std::array<Eigen::Vector3d, 2> ObtainLine(double omega) const;
    };

public:
    constexpr static double DEG2RAD = M_PI / 180.0;
    constexpr static double P2L_DISTANCE_THD = 2.0 /*pixel*/;
    constexpr static double L2L_DISTANCE_THD = 5.0 /*pixel*/;
    constexpr static double P2L_ORIENTATION_THD = std::cos(1.0 /*degree*/ * DEG2RAD);
    constexpr static double L2L_ORIENTATION_THD = std::cos(3.0 /*degree*/ * DEG2RAD);
    constexpr static int MAX_LINE_COUNT = 20;
    constexpr static double VISIBLE_LINE_ACTIVITY_THD = 3.0;

private:
    std::set<EventLine::Ptr> visibleLines;
    std::map<EventLine::Ptr, cv::Scalar> visibleLineColors;

public:
    EventLineTracking() {}

    static Ptr Create() { return std::make_shared<EventLineTracking>(); }

    void TrackingUsingNormFlow(const EventNormFlow::NormFlowPack::Ptr &nfPack);

protected:
    static void DrawLine(cv::Mat &m,
                         const EventLine::Ptr &el,
                         const cv::Scalar &color = cv::Scalar(0, 255, 0));

    static void DrawLine(cv::Mat &m,
                         const double a,
                         const double b,
                         const double c,
                         const cv::Scalar &color = cv::Scalar(0, 255, 0));

    static Eigen::Vector3d EstimateLine(const EventLine::Ptr &ol,
                                        const std::list<NormFlowPtr> &data,
                                        double lambda);

    static EventLine::Ptr FindPossibleSameLine(const std::set<EventLine::Ptr> &lines,
                                               const EventLine::Ptr &ml);

    static std::list<std::pair<EventLine::Ptr, EventLine::Ptr>> FilterSameLines(
        std::set<EventLine::Ptr> &lines);
};

}  // namespace ns_ikalibr

#endif  // EVENT_PREPROCESSING_H
