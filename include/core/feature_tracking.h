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

#ifndef FEATURE_TRACKING_H
#define FEATURE_TRACKING_H

#include "utility"
#include "util/utils.h"
#include "veta/camera/pinhole.h"
#include "opencv2/features2d.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
class CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;

struct Feature {
    cv::Point2f raw;
    cv::Point2f undistorted;

    Feature(cv::Point2f raw, cv::Point2f undistorted);
};

using FeatureMap = std::map<int, Feature>;
using FeatureVec = std::vector<Feature>;
using FeatureIdVec = std::vector<int>;
using FeatureMatch = std::map<int, int>;

class FeatureTracking {
public:
    using Ptr = std::shared_ptr<FeatureTracking>;

    struct TrackedFeaturePack {
        using Ptr = std::shared_ptr<TrackedFeaturePack>;

        CameraFramePtr imgLast;
        CameraFramePtr imgCur;
        FeatureMap featLast;
        FeatureMap featCur;
        FeatureMatch featMatchLast2Cur;

        [[nodiscard]] cv::Mat DrawMatches() const;

        [[nodiscard]] cv::Mat DrawMatches(const Ptr& compPack) const;

    protected:
        static void DrawFeatureTracking(const std::map<int, Feature>& ptsInLast,
                                        const std::map<int, Feature>& ptsInCur,
                                        const std::map<int, int>& matches,
                                        cv::Mat& matImg,
                                        const cv::Point2f& bias,
                                        const cv::Scalar& color);
    };

protected:
    const int FEAT_NUM_PER_IMG;
    const int MIN_DIST;

    // visual intrinsics
    ns_veta::PinholeIntrinsic::Ptr _intri;
    // the last image
    CameraFramePtr _imgLast;
    // feature id, raw feature, undistorted feature in the last image
    FeatureMap _featLast;
    // feature id, tracked count
    std::map<int, int> _featTrackCountLast;

public:
    explicit FeatureTracking(int featNumPerImg,
                             int minDist,
                             ns_veta::PinholeIntrinsic::Ptr intri = nullptr);

    virtual ~FeatureTracking() = default;

    TrackedFeaturePack::Ptr GrabImageFrame(
        const CameraFramePtr& imgCur,
        const std::optional<Sophus::SO3d>& SO3_Last2Cur = std::nullopt);

protected:
    /**
     * @param imgCur the current image
     * @return the
     */
    virtual void ExtractFeatures(const CameraFramePtr& imgCur,
                                 const cv::Mat& mask,
                                 int featCountDesired,
                                 std::vector<cv::Point2f>& ptsCurVec,
                                 FeatureIdVec& ptsCurIdVec,
                                 int& ptsIdCounter) = 0;

    /**
     * @param imgLast the last image
     * @param ptsLastVec the features in the last image
     * @param ptsLastIdVec the ids of features in the last image
     * @param imgCur the current image
     * @param ptsCurVec the features in the current image
     * @param status the tracking status
     * @param SO3_Last2Cur the prior rotation
     */
    virtual void GrabNextImageFrame(const CameraFramePtr& imgLast,
                                    const std::vector<cv::Point2f>& ptsLastVec,
                                    const FeatureIdVec& ptsLastIdVec,
                                    const std::optional<Sophus::SO3d>& SO3_Last2Cur,
                                    const CameraFramePtr& imgCur,
                                    std::vector<cv::Point2f>& ptsCurVec,
                                    FeatureIdVec& ptsCurIdVec,
                                    std::vector<uchar>& status,
                                    int& ptsIdCounter) = 0;

    static void ComputeIndexVecOfPoints(const std::vector<cv::Point2f>& ptsVec,
                                        FeatureIdVec& ptsIdVec,
                                        int& ptsIdCounter);

    [[nodiscard]] cv::Point2f UndistortPoint(const cv::Point2f& p) const;

    void ComputePriorPoints(const std::vector<cv::Point2f>& ptsLast,
                            const Sophus::SO3d& SO3_Last2Cur,
                            std::vector<cv::Point2f>& ptsCur) const;

    static bool InImageBorder(const cv::Point2f& pt, const cv::Mat& img, int borderSize);

    [[nodiscard]] std::pair<cv::Mat, FeatureMap> MakeTrackedFeatUniform(
        const FeatureMap& feats,
        const CameraFramePtr& frame,
        const std::vector<int>& featPriority) const;
};

class LKFeatureTracking : public FeatureTracking {
public:
    using Ptr = std::shared_ptr<LKFeatureTracking>;

public:
    explicit LKFeatureTracking(int featNumPerImg,
                               int minDist,
                               const ns_veta::PinholeIntrinsic::Ptr& intri = nullptr);

    static Ptr Create(int featNumPerImg,
                      int minDist,
                      const ns_veta::PinholeIntrinsic::Ptr& intri = nullptr);

protected:
    void ExtractFeatures(const CameraFramePtr& imgCur,
                         const cv::Mat& mask,
                         int featCountDesired,
                         std::vector<cv::Point2f>& ptsCurVec,
                         FeatureIdVec& ptsCurIdVec,
                         int& ptsIdCounter) override;

    void GrabNextImageFrame(const CameraFramePtr& imgLast,
                            const std::vector<cv::Point2f>& ptsLastVec,
                            const FeatureIdVec& ptsLastIdVec,
                            const std::optional<Sophus::SO3d>& SO3_Last2Cur,
                            const CameraFramePtr& imgCur,
                            std::vector<cv::Point2f>& ptsCurVec,
                            FeatureIdVec& ptsCurIdVec,
                            std::vector<uchar>& status,
                            int& ptsIdCounter) override;
};

class DescriptorBasedFeatureTracking : public FeatureTracking {
public:
    using Ptr = std::shared_ptr<DescriptorBasedFeatureTracking>;

protected:
    cv::Ptr<cv::DescriptorMatcher> _matcher;
    std::map<int, cv::KeyPoint> _kptLastMap;
    std::map<int, cv::Mat> _descLastMap;

    static constexpr double NN_MATCH_RATION = 0.8f;

public:
    DescriptorBasedFeatureTracking(
        int featNumPerImg,
        int minDist,
        const ns_veta::PinholeIntrinsic::Ptr& intri = nullptr,
        const cv::Ptr<cv::DescriptorMatcher>& matcher = cv::BFMatcher::create(cv::NORM_HAMMING));

protected:
    void ExtractFeatures(const CameraFramePtr& imgCur,
                         const cv::Mat& mask,
                         int featCountDesired,
                         std::vector<cv::Point2f>& ptsCurVec,
                         FeatureIdVec& ptsCurIdVec,
                         int& ptsIdCounter) override;

    void GrabNextImageFrame(const CameraFramePtr& imgLast,
                            const std::vector<cv::Point2f>& ptsLastVec,
                            const FeatureIdVec& ptsLastIdVec,
                            const std::optional<Sophus::SO3d>& SO3_Last2Cur,
                            const CameraFramePtr& imgCur,
                            std::vector<cv::Point2f>& ptsCurVec,
                            FeatureIdVec& ptsCurIdVec,
                            std::vector<uchar>& status,
                            int& ptsIdCounter) override;

    virtual void DetectAndComputeKeyPoints(const cv::Mat& img,
                                           const cv::Mat& mask,
                                           int featCountDesired,
                                           std::vector<cv::KeyPoint>& kps,
                                           cv::Mat& descriptor) = 0;
};

class ORBFeatureTracking : public DescriptorBasedFeatureTracking {
public:
    using Ptr = std::shared_ptr<ORBFeatureTracking>;

public:
    ORBFeatureTracking(
        int featNumPerImg,
        int minDist,
        const ns_veta::PinholeIntrinsic::Ptr& intri = nullptr,
        const cv::Ptr<cv::DescriptorMatcher>& matcher = cv::BFMatcher::create(cv::NORM_HAMMING));

    static Ptr Create(
        int featNumPerImg,
        int minDist,
        const ns_veta::PinholeIntrinsic::Ptr& intri,
        const cv::Ptr<cv::DescriptorMatcher>& matcher = cv::BFMatcher::create(cv::NORM_HAMMING));

protected:
    void DetectAndComputeKeyPoints(const cv::Mat& img,
                                   const cv::Mat& mask,
                                   int featCountDesired,
                                   std::vector<cv::KeyPoint>& kps,
                                   cv::Mat& descriptor) override;
};

class AKAZEFeatureTracking : public DescriptorBasedFeatureTracking {
public:
    using Ptr = std::shared_ptr<AKAZEFeatureTracking>;

public:
    AKAZEFeatureTracking(
        int featNumPerImg,
        int minDist,
        const ns_veta::PinholeIntrinsic::Ptr& intri = nullptr,
        const cv::Ptr<cv::DescriptorMatcher>& matcher = cv::BFMatcher::create(cv::NORM_HAMMING));

    static Ptr Create(
        int featNumPerImg,
        int minDist,
        const ns_veta::PinholeIntrinsic::Ptr& intri,
        const cv::Ptr<cv::DescriptorMatcher>& matcher = cv::BFMatcher::create(cv::NORM_HAMMING));

protected:
    void DetectAndComputeKeyPoints(const cv::Mat& img,
                                   const cv::Mat& mask,
                                   int featCountDesired,
                                   std::vector<cv::KeyPoint>& kps,
                                   cv::Mat& descriptor) override;
};
}  // namespace ns_ikalibr

#endif  // FEATURE_TRACKING_H
