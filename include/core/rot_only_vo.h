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

#include "util/utils.h"
#include "opencv2/core.hpp"
#include "opengv/types.hpp"
#include "veta/camera/pinhole.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;

struct ORBFeatExtMatConfig {
public:
    int featNumPerImg{};

    ORBFeatExtMatConfig();
};

class RotOnlyVisualOdometer {
public:
    using Ptr = std::shared_ptr<RotOnlyVisualOdometer>;

private:
    const int FEAT_NUM_PER_IMG;
    const int MIN_DIST;
    ns_veta::PinholeIntrinsic::Ptr _intri;

    CameraFramePtr _lastFrame;
    std::vector<cv::Point2f> _ptsInLast;
    std::vector<int> _ptsTrackCount;
    std::vector<cv::Point2f> _ptsUndistInLast;

    std::vector<std::pair<double, Sophus::SO3d>> _rotations;

public:
    explicit RotOnlyVisualOdometer(int featNumPerImg,
                                   int minDist,
                                   ns_veta::PinholeIntrinsic::Ptr intri)
        : FEAT_NUM_PER_IMG(featNumPerImg),
          MIN_DIST(minDist),
          _intri(std::move(intri)),
          _lastFrame(nullptr) {}

    static Ptr Create(int featNumPerImg, int minDist, const ns_veta::PinholeIntrinsic::Ptr &intri);

    bool GrabFrame(const CameraFramePtr &curFrame);

    [[nodiscard]] const std::vector<std::pair<double, Sophus::SO3d>> &GetRotations() const;

    virtual ~RotOnlyVisualOdometer();

protected:
    static bool InImageBorder(const cv::Point2f &pt,
                              const CameraFramePtr &frame,
                              int borderSize = 1);

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

    std::vector<cv::Point2f> UndistortedPoints(const std::vector<cv::Point2f> &pts);

    [[nodiscard]] std::tuple<cv::Mat, std::vector<cv::Point2f>, std::vector<int>, std::vector<int>>
    ComputeMaskAndFilterPts(const CameraFramePtr &frame,
                            const std::vector<cv::Point2f> &pts,
                            const std::vector<int> &trackCount) const;

    template <class Type>
    std::vector<Type> FindElements(const std::vector<Type> &vec, const std::vector<int> &idx) {
        std::vector<Type> newVec(idx.size());
        for (int i = 0; i < static_cast<int>(idx.size()); ++i) {
            newVec.at(i) = vec.at(idx.at(i));
        }
        return newVec;
    }

    void ShowCurrentFrame() const;

    static void ShowFeatureTracking(const CameraFramePtr &lastFrame,
                                    const std::vector<cv::Point2f> &ptsInLast,
                                    const CameraFramePtr &curFrame,
                                    const std::vector<cv::Point2f> &ptsInCur,
                                    const std::vector<int> &inliers,
                                    const std::string &winName);

    static std::vector<uchar> RejectUsingFMat(const std::vector<cv::Point2f> &undistPtsInLast,
                                              const std::vector<cv::Point2f> &undistPtsInCur);

    std::pair<opengv::rotation_t, std::vector<int>> RelRotationRecovery(
        const std::vector<cv::Point2f> &ptsUndisto1, const std::vector<cv::Point2f> &ptsUndisto2);

    opengv::bearingVectors_t ComputeBeringVec(const std::vector<cv::Point2f> &ptsUndist);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_ROT_ONLY_VO_H
