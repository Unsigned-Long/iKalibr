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

#ifndef IKALIBR_ROT_ONLY_VO_H
#define IKALIBR_ROT_ONLY_VO_H

#include "utility"
#include "util/utils.h"
#include "opencv2/core.hpp"
#include "opengv/types.hpp"
#include "veta/camera/pinhole.h"
#include "core/feature_tracking.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
class CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;

class RotOnlyVisualOdometer {
public:
    using Ptr = std::shared_ptr<RotOnlyVisualOdometer>;

    // landmark id, track lists [camera frame, feature point]
    using FeatTrackingInfo =
        std::map<ns_veta::IndexT, std::list<std::pair<CameraFramePtr, Feature>>>;

private:
    FeatureTracking::Ptr _featTracking;
    ns_veta::PinholeIntrinsic::Ptr _intri;
    FeatureTracking::TrackedFeaturePack::Ptr _trackFeatLast;

    // landmark id, track lists [camera frame, feature point]
    std::map<ns_veta::IndexT, std::list<std::pair<CameraFramePtr, Feature>>> _lmTrackInfo;
    // feature id, landmark id, only for the last image
    std::map<int, ns_veta::IndexT> _featId2lmIdInLast;

    std::vector<std::pair<double, Sophus::SO3d>> _rotations;

public:
    explicit RotOnlyVisualOdometer(FeatureTracking::Ptr featTracking,
                                   ns_veta::PinholeIntrinsic::Ptr intri);

    static Ptr Create(const FeatureTracking::Ptr &featTracking,
                      const ns_veta::PinholeIntrinsic::Ptr &intri);

    bool GrabFrame(const CameraFramePtr &curFrame,
                   const std::optional<Sophus::SO3d> &SO3_LastToCur = std::nullopt);

    [[nodiscard]] const std::vector<std::pair<double, Sophus::SO3d>> &GetRotations() const;

    virtual ~RotOnlyVisualOdometer();

    [[nodiscard]] const std::map<ns_veta::IndexT, std::list<std::pair<CameraFramePtr, Feature>>> &
    GetLmTrackInfo() const;

    void ShowLmTrackInfo() const;

    void ResetWorkspace();

protected:
    static std::pair<std::vector<int>, std::vector<cv::Point2f>> ExtractFeatMapAsRawFeatVec(
        const FeatureMap &featMap, const std::vector<int> &desiredIds = {});

    static std::pair<std::vector<int>, std::vector<cv::Point2f>> ExtractFeatMapAsUndistoFeatVec(
        const FeatureMap &featMap, const std::vector<int> &desiredIds = {});

    static ns_veta::IndexT GenNewLmId();

    template <class Type>
    static void ReduceVector(std::vector<Type> &v, std::vector<uchar> status) {
        int j = 0;
        for (int i = 0; i < int(v.size()); i++) {
            if (status[i]) {
                v[j++] = v[i];
            }
        }
        v.resize(j);
    }

    [[nodiscard]] std::pair<cv::Mat, std::set<int>> ComputeMaskAndFilterPts(
        const CameraFramePtr &frame,
        const FeatureMap &featMap,
        std::vector<std::pair<int, int>> trackCount) const;

    template <class Type>
    std::vector<Type> FindElements(const std::vector<Type> &vec, const std::vector<int> &idx) {
        std::vector<Type> newVec(idx.size());
        for (int i = 0; i < static_cast<int>(idx.size()); ++i) {
            newVec.at(i) = vec.at(idx.at(i));
        }
        return newVec;
    }

    void ShowCurrentFrame() const;

    static std::vector<uchar> RejectUsingFMat(const std::vector<cv::Point2f> &undistPtsInLast,
                                              const std::vector<cv::Point2f> &undistPtsInCur);

    [[nodiscard]] std::pair<opengv::rotation_t, std::vector<int>> RelRotationRecovery(
        const std::vector<cv::Point2f> &ptsUndisto1,
        const std::vector<cv::Point2f> &ptsUndisto2) const;

    [[nodiscard]] opengv::bearingVectors_t ComputeBeringVec(
        const std::vector<cv::Point2f> &ptsUndist) const;
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_ROT_ONLY_VO_H
