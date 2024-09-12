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

#include "core/feature_tracking.h"
#include "sensor/camera.h"
#include "opencv2/imgproc.hpp"
#include "opencv2/video/tracking.hpp"

namespace ns_ikalibr {

Feature::Feature(cv::Point2f raw, cv::Point2f undistorted)
    : raw(std::move(raw)),
      undistorted(std::move(undistorted)) {}

FeatureTracking::FeatureTracking(int featNumPerImg,
                                 int minDist,
                                 ns_veta::PinholeIntrinsic::Ptr intri)
    : FEAT_NUM_PER_IMG(featNumPerImg),
      MIN_DIST(minDist),
      _intri(std::move(intri)),
      _imgLast(nullptr),
      _featLast() {}

FeatureTracking::TrackedFeaturePack::Ptr FeatureTracking::GrabImageFrame(
    const CameraFrame::Ptr& imgCur, const std::optional<Sophus::SO3d>& SO3_Last2Cur) {
    if (_imgLast == nullptr) {
        // for the first image, we just extract features and store them to '_featLast'
        std::vector<cv::Point2f> ptsLastVec;
        ExtractFeatures(imgCur, ptsLastVec, FEAT_NUM_PER_IMG, cv::Mat());
        int idCount = 0;
        FeatureMap featCur;
        for (const auto& raw : ptsLastVec) {
            featCur.insert({idCount++, Feature(raw, UndistortPoint(raw))});
        }
        TrackedFeaturePack::Ptr pack = std::make_shared<TrackedFeaturePack>();
        pack->imgLast = nullptr;
        pack->featLast = {};
        pack->featMatchLast2Cur = {};
        pack->imgCur = imgCur;
        pack->featCur = featCur;

        _imgLast = imgCur;
        _featLast = featCur;
        for (const auto& [id, feat] : _featLast) {
            _featTrackCountLast.insert({id, 1});
        }
        return pack;
    }
    // last tracking info
    std::vector<cv::Point2f> ptsLastVec;
    ptsLastVec.reserve(_featLast.size());
    FeatureIdVec ptsLastIdVec;
    ptsLastIdVec.reserve(_featLast.size());
    for (const auto& [id, feat] : _featLast) {
        ptsLastVec.push_back(feat.raw);
        ptsLastIdVec.push_back(id);
    }
    // current tracking
    std::vector<cv::Point2f> ptsCurVec;
    std::vector<uchar> status;
    GrabNextImageFrame(_imgLast, ptsLastVec, ptsLastIdVec, imgCur, ptsCurVec, status, SO3_Last2Cur);

    FeatureMap featCur;
    FeatureMatch featMatchLast2Cur;
    int idCount = 0;
    for (int i = 0; i < static_cast<int>(status.size()); ++i) {
        const cv::Point2f& pCur = ptsCurVec.at(i);
        if (status.at(i) && InImageBorder(pCur, imgCur->GetImage(), 5)) {
            auto idCur = idCount++;
            featCur.insert({idCur, Feature(pCur, UndistortPoint(pCur))});
            featMatchLast2Cur.insert({ptsLastIdVec.at(i), idCur});
        }
    }

    // compute the tracking count of the current image, id, count
    std::map<int, int> featTrackCountCur;
    for (const auto& [idLast, idCur] : featMatchLast2Cur) {
        featTrackCountCur.insert({idCur, _featTrackCountLast.at(idLast) + 1});
    }

    // filter tracked features, id, track count
    std::vector<std::pair<int, int>> trackCount;
    trackCount.reserve(featTrackCountCur.size());
    for (const auto& [id, count] : featTrackCountCur) {
        trackCount.emplace_back(id, count);
    }
    // sort features based on the tracking count
    std::sort(trackCount.begin(), trackCount.end(),
              [](const auto& p1, const auto& p2) { return p1.second > p2.second; });

    std::vector<int> featPriority;
    featPriority.reserve(trackCount.size());
    for (const auto& [id, feat] : trackCount) {
        featPriority.push_back(id);
    }
    auto [mask, filteredFeatCur] = MakeTrackedFeatUniform(featCur, imgCur, featPriority);
    // update 'featCur'
    featCur = filteredFeatCur;
    // update 'featTrackCountCur'
    for (auto iter = featTrackCountCur.begin(); iter != featTrackCountCur.end();) {
        auto [idCur, count] = *iter;
        if (featCur.find(idCur) == featCur.end()) {
            // this feature has been filtered
            iter = featTrackCountCur.erase(iter);
        } else {
            ++iter;
        }
    }
    // update 'featMatchLast2Cur'
    for (auto iter = featMatchLast2Cur.begin(); iter != featMatchLast2Cur.end();) {
        auto [idLast, idCur] = *iter;
        if (featCur.find(idCur) == featCur.end()) {
            iter = featMatchLast2Cur.erase(iter);
        } else {
            ++iter;
        }
    }
    // update '_featLast'
    for (auto iter = _featLast.begin(); iter != _featLast.end();) {
        auto [idLast, feat] = *iter;
        if (featMatchLast2Cur.find(idLast) == featMatchLast2Cur.end()) {
            iter = _featLast.erase(iter);
        } else {
            ++iter;
        }
    }

    TrackedFeaturePack::Ptr pack = std::make_shared<TrackedFeaturePack>();
    pack->imgLast = _imgLast;
    pack->featLast = _featLast;
    pack->imgCur = imgCur;
    pack->featCur = featCur;
    pack->featMatchLast2Cur = featMatchLast2Cur;

    // extract new features
    int featNumToExtract = FEAT_NUM_PER_IMG - static_cast<int>(featCur.size());
    if (featNumToExtract > 0) {
        std::vector<cv::Point2f> newPtsCurVec;
        ExtractFeatures(imgCur, newPtsCurVec, featNumToExtract, mask);
        /**
         * append new features to map, attention, we continue to use the 'idCount' above, rather
         * than redefine 'int idCount = 0;'
         */
        for (const auto& raw : newPtsCurVec) {
            auto id = idCount++;
            featCur.insert({id, Feature(raw, UndistortPoint(raw))});
            featTrackCountCur.insert({id, 1});
        }
    }

    _imgLast = imgCur;
    _featLast = featCur;
    _featTrackCountLast = featTrackCountCur;

    return pack;
}

cv::Point2f FeatureTracking::UndistortPoint(const cv::Point2f& p) const {
    ns_veta::Vec2d up = _intri->GetUndistoPixel(ns_veta::Vec2d(p.x, p.y));
    return {static_cast<float>(up(0)), static_cast<float>(up(1))};
}

void FeatureTracking::ComputePriorPoints(const std::vector<cv::Point2f>& ptsLast,
                                         const Sophus::SO3d& SO3_Last2Cur,
                                         std::vector<cv::Point2f>& ptsCur) const {
    ptsCur.clear();
    ptsCur.reserve(ptsLast.size());
    for (const auto& raw : ptsLast) {
        Eigen::Vector2d pCam = _intri->RemoveDisto(_intri->ImgToCam({raw.x, raw.y}));
        Eigen::Vector3d pCamNew = SO3_Last2Cur * Eigen::Vector3d(pCam(0), pCam(1), 1.0);
        Eigen::Vector2d rawNew = _intri->CamToImg(_intri->AddDisto({pCamNew(0), pCamNew(1)}));
        ptsCur.emplace_back(rawNew(0), rawNew(1));
    }
}

bool FeatureTracking::InImageBorder(const cv::Point2f& pt, const cv::Mat& img, int borderSize) {
    int col = img.cols, row = img.rows;
    int imgX = cvRound(pt.x), imgY = cvRound(pt.y);
    return borderSize <= imgX && imgX < col - borderSize && borderSize <= imgY &&
           imgY < row - borderSize;
}

std::pair<cv::Mat, FeatureMap> FeatureTracking::MakeTrackedFeatUniform(
    const FeatureMap& feats,
    const CameraFramePtr& frame,
    const std::vector<int>& featPriority) const {
    int col = frame->GetImage().cols, row = frame->GetImage().rows;
    cv::Mat mask = cv::Mat(row, col, CV_8UC1, cv::Scalar(255));
    std::set<int> filteredFeatId;
    FeatureMap filteredFeatLast;
    for (const auto& featId : featPriority) {
        auto feat = feats.at(featId);
        // if this not has been occupied
        if (mask.at<uchar>(feat.raw) == 255) {
            // draw mask
            cv::circle(mask, feat.raw, MIN_DIST, 0, -1);
            filteredFeatLast.insert({featId, feat});
        }
    }
    return {mask, filteredFeatLast};
}

/**
 * feature tracking based on LK optical flow
 */
LKFeatureTracking::LKFeatureTracking(int featNumPerImg,
                                     int minDist,
                                     const ns_veta::PinholeIntrinsic::Ptr& intri)
    : FeatureTracking(featNumPerImg, minDist, intri) {}

LKFeatureTracking::Ptr LKFeatureTracking::Create(int featNumPerImg,
                                                 int minDist,
                                                 const ns_veta::PinholeIntrinsic::Ptr& intri) {
    return std::make_shared<LKFeatureTracking>(featNumPerImg, minDist, intri);
}

void LKFeatureTracking::ExtractFeatures(const CameraFrame::Ptr& imgCur,
                                        std::vector<cv::Point2f>& featCur,
                                        int featCountDesired,
                                        const cv::Mat& mask) {
    // the mask is empty for the first frame, do not use mask
    cv::goodFeaturesToTrack(imgCur->GetImage(), featCur, featCountDesired, 0.01, MIN_DIST, mask);
}

void LKFeatureTracking::GrabNextImageFrame(const CameraFramePtr& imgLast,
                                           const std::vector<cv::Point2f>& ptsLastVec,
                                           const FeatureIdVec& ptsLastIdVec,
                                           const CameraFramePtr& imgCur,
                                           std::vector<cv::Point2f>& ptsCurVec,
                                           std::vector<uchar>& status,
                                           const std::optional<Sophus::SO3d>& SO3_Last2Cur) {
    if (SO3_Last2Cur != std::nullopt) {
        ComputePriorPoints(ptsLastVec, *SO3_Last2Cur, ptsCurVec);
    } else {
        ptsCurVec = ptsLastVec;
    }
    std::vector<float> errors;
    cv::TermCriteria termCrit =
        cv::TermCriteria(cv::TermCriteria::COUNT | cv::TermCriteria::EPS, 30, 0.01);
    cv::calcOpticalFlowPyrLK(imgLast->GetImage(), imgCur->GetImage(), ptsLastVec, ptsCurVec, status,
                             errors, cv::Size(21, 21), 5, termCrit, cv::OPTFLOW_USE_INITIAL_FLOW);
}

}  // namespace ns_ikalibr