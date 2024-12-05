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

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct NormFlow;
using NormFlowPtr = std::shared_ptr<NormFlow>;

class EventCircleTracking {
public:
    using Ptr = std::shared_ptr<EventCircleTracking>;

    constexpr static int CLUSTER_AREA_THD = 10;
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

public:
    EventCircleTracking() = default;

    static Ptr Create() { return std::make_shared<EventCircleTracking>(); }

    void ExtractCircles(const EventNormFlow::NormFlowPack::Ptr& nfPack);

protected:
    static std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> SearchMatchesInRunChasePair(
        const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
        double DIR_DIFF_COS_THD);

    static std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> ReSearchMatchesCirclesOtherPair(
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
    ;

    static void RemovingAmbiguousMatches(
        std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr>& pairs);

    template <typename Type1, typename Type2>
    static void RemoveClusterTypes(std::vector<CircleClusterType>& pClusterType,
                                   std::vector<Type1>& seq1,
                                   std::vector<Type2>& seq2,
                                   CircleClusterType typeToRemove) {
        assert(pClusterType.size() == seq.size());
        const auto size = pClusterType.size();

        std::vector<CircleClusterType> newPClusterType;
        newPClusterType.reserve(size);
        std::vector<Type1> newSeq1;
        newSeq1.reserve(size);
        std::vector<Type2> newSeq2;
        newSeq2.reserve(size);

        for (size_t i = 0; i < pClusterType.size(); ++i) {
            if (pClusterType[i] != typeToRemove) {
                newPClusterType.push_back(pClusterType[i]);
                newSeq1.push_back(seq1[i]);
                newSeq2.push_back(seq2[i]);
            }
        }

        pClusterType = std::move(newPClusterType);
        seq1 = std::move(newSeq1);
        seq2 = std::move(newSeq2);
    }

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

    static std::pair<std::vector<std::list<NormFlowPtr>>, std::vector<std::list<NormFlowPtr>>>
    ClusterNormFlowEvents(const EventNormFlow::NormFlowPack::Ptr& nfPack, double clusterAreaThd);

    static void InterruptionInTimeDomain(cv::Mat& pMat, const cv::Mat& tMat, double thd);

    static void FilterContoursUsingArea(std::vector<std::vector<cv::Point>>& contours, int areaThd);

    static std::vector<std::vector<cv::Point>> FindContours(const cv::Mat& binaryImg);

protected:
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
