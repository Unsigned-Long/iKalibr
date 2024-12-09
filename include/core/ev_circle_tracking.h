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

#ifndef EV_CIRCLE_TRACKING_H
#define EV_CIRCLE_TRACKING_H

#include "util/utils.h"
#include "core/event_preprocessing.h"

#include <ostream>

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct NormFlow;
using NormFlowPtr = std::shared_ptr<NormFlow>;
struct Viewer;
using ViewerPtr = std::shared_ptr<Viewer>;

class EventCircleTracking {
public:
    using Ptr = std::shared_ptr<EventCircleTracking>;

    constexpr static double DEG2RAD = M_PI / 180.0;

    enum class CircleClusterType : int { CHASE = 0, RUN = 1, OTHER = 2 };

    struct CircleClusterInfo {
        using Ptr = std::shared_ptr<CircleClusterInfo>;

        std::list<NormFlowPtr> nfCluster;
        bool polarity;
        Eigen::Vector2d center;
        Eigen::Vector2d dir;

        CircleClusterInfo(const std::list<NormFlowPtr>& nf_cluster,
                          bool polarity,
                          const Eigen::Vector2d& center,
                          const Eigen::Vector2d& dir);

        static Ptr Create(const std::list<NormFlowPtr>& nf_cluster,
                          bool polarity,
                          const Eigen::Vector2d& center,
                          const Eigen::Vector2d& dir);
    };

    struct TimeVaryingCircle {
        using Ptr = std::shared_ptr<TimeVaryingCircle>;
        // center, radius
        using Circle = std::pair<Eigen::Vector2d, double>;

        double st, et;
        Eigen::Vector2d cx;
        Eigen::Vector2d cy;
        Eigen::Vector3d r2;

        TimeVaryingCircle(double st,
                          double et,
                          const Eigen::Vector2d& cx,
                          const Eigen::Vector2d& cy,
                          const Eigen::Vector3d& r2)
            : st(st),
              et(et),
              cx(cx),
              cy(cy),
              r2(r2) {}

        static Ptr Create(double st,
                          double et,
                          const Eigen::Vector2d& cx,
                          const Eigen::Vector2d& cy,
                          const Eigen::Vector3d& r2) {
            return std::make_shared<TimeVaryingCircle>(st, et, cx, cy, r2);
        }

        Eigen::Vector2d PosAt(double t) const {
            Eigen::Vector2d tVec(t, 1.0);
            return {cx.dot(tVec), cy.dot(tVec)};
        }

        std::vector<Eigen::Vector3d> PosVecAt(double dt) const {
            std::list<Eigen::Vector3d> posList;
            double t = st;
            while (t < et) {
                // t, x, y
                Eigen::Vector2d p = PosAt(t);
                posList.push_back({t, p(0), p(1)});
                t += dt;
            }
            return std::vector<Eigen::Vector3d>{posList.cbegin(), posList.cend()};
        }

        double RadiusAt(double t) const {
            Eigen::Vector3d tVec(t * t, t, 1.0);
            return std::sqrt(r2.dot(tVec));
        }

        Circle CircleAt(double t) const { return {PosAt(t), RadiusAt(t)}; }

        double PointToCircleDistance(const Event::Ptr& event) const {
            auto c = CircleAt(event->GetTimestamp());
            return std::abs(c.second - (event->GetPos().cast<double>() - c.first).norm());
        }

        friend std::ostream& operator<<(std::ostream& os, const TimeVaryingCircle& obj) {
            return os << "st: " << obj.st << " et: " << obj.et << " cx: " << obj.cx.transpose()
                      << " cy: " << obj.cy.transpose() << " r2: " << obj.r2.transpose();
        }
    };

protected:
    const double CLUSTER_AREA_THD;
    const double DIR_DIFF_DEG_THD;
    const double POINT_TO_CIRCLE_AVG_THD;

public:
    EventCircleTracking(double CLUSTER_AREA_THD,
                        double DIR_DIFF_DEG_THD,
                        double POINT_TO_CIRCLE_AVG_THD)
        : CLUSTER_AREA_THD(CLUSTER_AREA_THD),
          DIR_DIFF_DEG_THD(DIR_DIFF_DEG_THD),
          POINT_TO_CIRCLE_AVG_THD(POINT_TO_CIRCLE_AVG_THD) {}

    static Ptr Create(double CLUSTER_AREA_THD = 10.0,
                      double DIR_DIFF_DEG_THD = 30.0,
                      double POINT_TO_CIRCLE_AVG_THD = 1.0) {
        return std::make_shared<EventCircleTracking>(CLUSTER_AREA_THD, DIR_DIFF_DEG_THD,
                                                     POINT_TO_CIRCLE_AVG_THD);
    }

    std::pair<double, std::vector<TimeVaryingCircle::Circle>> ExtractCircles(
        const EventNormFlow::NormFlowPack::Ptr& nfPack, const ViewerPtr& viewer = nullptr) const;

    void ExtractCirclesGrid(const EventNormFlow::NormFlowPack::Ptr& nfPack,
                            const ViewerPtr& viewer = nullptr) const;

protected:
    static std::vector<std::pair<EventArray::Ptr, EventArray::Ptr>> ExtractPotentialCircleClusters(
        const EventNormFlow::NormFlowPack::Ptr& nfPack,
        double CLUSTER_AREA_THD,
        double DIR_DIFF_DEG_THD);

    static TimeVaryingCircle::Ptr FitTimeVaryingCircle(const EventArray::Ptr& ary1,
                                                       const EventArray::Ptr& ary2,
                                                       double avgDistThd);

protected:
    static std::vector<std::pair<EventArray::Ptr, EventArray::Ptr>> RawEventsOfCircleClusterPairs(
        const std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr>& pairs,
        const EventNormFlow::NormFlowPack::Ptr& nfPack);

    /**
     * The functions 'SearchMatchesInRunChasePair', 'ReSearchMatchesCirclesOtherPair', and
     * 'ReSearchMatchesOtherOtherPair' are used to search for potential circular matching clusters.
     */
    static std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> SearchMatchesInRunChasePair(
        const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
        double DIR_DIFF_COS_THD);

    static std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> ReSearchMatchesCirclesOtherPair(
        const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
        const std::set<CircleClusterInfo::Ptr>& alreadyMatched,
        double DIR_DIFF_COS_THD);

    static std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> ReSearchMatchesOtherOtherPair(
        const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
        const std::set<CircleClusterInfo::Ptr>& alreadyMatched,
        double DIR_DIFF_COS_THD);

    static double TryMatchRunChaseCircleClusterPair(const CircleClusterInfo::Ptr& rCluster,
                                                    const CircleClusterInfo::Ptr& cCluster,
                                                    double DIR_DIFF_COS_THD);

    static bool ClusterExistsInCurCircle(
        const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
        const CircleClusterInfo::Ptr& p1,
        const CircleClusterInfo::Ptr& p2);

    static void RemovingAmbiguousMatches(
        std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr>& pairs);

    /**
     * The classification of the clusters is determined as either the 'chase' category, the 'run'
     * category, or the 'other' category.
     */
    static std::vector<CircleClusterType> IdentifyCategory(
        const std::vector<std::list<NormFlowPtr>>& clusters,
        const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& cenDirs,
        const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static CircleClusterType IdentifyCategory(
        const std::list<NormFlowPtr>& clusters,
        const std::pair<Eigen::Vector2d, Eigen::Vector2d>& cenDir,
        const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> ComputeCenterDir(
        const std::vector<std::list<NormFlowPtr>>& clusters,
        const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static std::pair<Eigen::Vector2d, Eigen::Vector2d> ComputeCenterDir(
        const std::list<NormFlowPtr>& cluster, const EventNormFlow::NormFlowPack::Ptr& nfPack);

    /**
     * The following functions serve the purpose of 'ClusterNormFlowEvents'.
     */
    static std::pair<std::vector<std::list<NormFlowPtr>>, std::vector<std::list<NormFlowPtr>>>
    ClusterNormFlowEvents(const EventNormFlow::NormFlowPack::Ptr& nfPack, double clusterAreaThd);

    static void InterruptionInTimeDomain(cv::Mat& pMat, const cv::Mat& tMat, double thd);

    static void FilterContoursUsingArea(std::vector<std::vector<cv::Point>>& contours,
                                        double areaThd);

    static std::vector<std::vector<cv::Point>> FindContours(const cv::Mat& binaryImg);

protected:
    /**
     * The following functions are used for visualization.
     */
    static void DrawCircleCluster(
        cv::Mat& mat,
        const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
        const EventNormFlow::NormFlowPack::Ptr& nfPack,
        double scale);

    static void DrawCircleClusterPair(
        cv::Mat& mat,
        const std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr>& pair,
        const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static void DrawCircleCluster(cv::Mat& mat,
                                  const std::vector<CircleClusterInfo::Ptr>& clusters,
                                  CircleClusterType type,
                                  const EventNormFlow::NormFlowPack::Ptr& nfPack,
                                  double scale);

    static void DrawCluster(cv::Mat& mat,
                            const std::vector<std::list<NormFlowPtr>>& clusters,
                            const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static void DrawCluster(cv::Mat& mat,
                            const std::list<NormFlowPtr>& clusters,
                            const EventNormFlow::NormFlowPack::Ptr& nfPack,
                            std::optional<ns_viewer::Colour> color = std::nullopt);

    static void DrawContours(cv::Mat& mat,
                             const std::vector<std::vector<cv::Point>>& contours,
                             const cv::Vec3b& color = {255, 255, 255});
};
}  // namespace ns_ikalibr

#endif  // EV_CIRCLE_TRACKING_H
