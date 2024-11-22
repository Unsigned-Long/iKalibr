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

class EventNormFlow {
public:
    using Ptr = std::shared_ptr<EventNormFlow>;
    struct NormFlow {
    public:
        using Ptr = std::shared_ptr<NormFlow>;

        double timestamp;
        Eigen::Vector2i p;
        Eigen::Vector2d nf;

        NormFlow(double timestamp, const Eigen::Vector2i &p, const Eigen::Vector2d &nf);

        static Ptr Create(double timestamp, const Eigen::Vector2i &p, const Eigen::Vector2d &nf);
    };

    struct NormFlowPack {
        std::list<NormFlow::Ptr> nfs;
        cv::Mat rawTimeSurfaceMap;  // ex, ey, et
        cv::Mat inliersOccupy;      // is inlier in norm flow estimation
        cv::Mat polarityMap;        // the polarity map
        // for visualization
        cv::Mat diagram;


    };

private:
    ActiveEventSurface::Ptr _sea;

public:
    explicit EventNormFlow(const ActiveEventSurface::Ptr &sea)
        : _sea(sea) {}

    NormFlowPack ExtractNormFlows(int winSize = 2,
                                  int neighborDist = 2,
                                  double goodRatioThd = 0.9,
                                  double timeDistEventToPlaneThd = 2E-3,
                                  int ransacMaxIter = 3) const;
};

class EventLocalPlaneSacProblem : public opengv::sac::SampleConsensusProblem<Eigen::Vector3d> {
public:
    typedef Eigen::Vector3d model_t;

public:
    explicit EventLocalPlaneSacProblem(const std::vector<std::tuple<int, int, double>> &data,
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
    const std::vector<std::tuple<int, int, double>> &_data;
};

}  // namespace ns_ikalibr

#endif  // EVENT_PREPROCESSING_H
