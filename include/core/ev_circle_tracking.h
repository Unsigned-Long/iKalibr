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

public:
    EventCircleTracking() = default;

    static Ptr Create() { return std::make_shared<EventCircleTracking>(); }

    void ExtractCircles(const EventNormFlow::NormFlowPack::Ptr& nfPack);

protected:
    static std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> ComputeCenterDir(
        const std::vector<std::list<NormFlowPtr>>& clusters,
        const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static std::pair<Eigen::Vector2d, Eigen::Vector2d> ComputeCenterDir(
        const std::list<NormFlowPtr>& cluster, const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static std::pair<std::vector<std::list<NormFlowPtr>>, std::vector<std::list<NormFlowPtr>>>
    ClusterNormFlowEvents(const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static void FilterContoursUsingArea(std::vector<std::vector<cv::Point>>& contours, int areaThd);

    static std::vector<std::vector<cv::Point>> FindContours(const cv::Mat& binaryImg);

protected:
    static void DrawCenterDir(
        cv::Mat& mat,
        const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& cenDirVec,
        double scale);

    static void DrawCluster(cv::Mat& mat,
                            const std::vector<std::list<NormFlowPtr>>& clusters,
                            const EventNormFlow::NormFlowPack::Ptr& nfPack);

    static void DrawContours(cv::Mat& mat,
                             const std::vector<std::vector<cv::Point>>& contours,
                             const cv::Vec3b& color = {255, 255, 255});
};
}  // namespace ns_ikalibr

#endif  // EV_CIRCLE_TRACKING_H
