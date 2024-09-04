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

#include "core/rot_only_vo.h"
#include "sensor/camera.h"
#include "calib/calib_param_manager.h"
#include "util/status.hpp"
#include "util/utils_tpl.hpp"

#include "opencv2/highgui.hpp"
#include "opencv2/video/tracking.hpp"
#include "opencv2/calib3d.hpp"

#include "opengv/sac/Ransac.hpp"
#include "opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp"
#include "opengv/relative_pose/CentralRelativeAdapter.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

ORBFeatExtMatConfig::ORBFeatExtMatConfig() = default;

RotOnlyVisualOdometer::Ptr RotOnlyVisualOdometer::Create(
    int featNumPerImg, int minDist, const ns_veta::PinholeIntrinsic::Ptr &intri) {
    return std::make_shared<RotOnlyVisualOdometer>(featNumPerImg, minDist, intri);
}

bool RotOnlyVisualOdometer::GrabFrame(const CameraFrame::Ptr &curFrame) {
#define VISUALIZATION 0

    // first curFrame
    if (_lastFrame == nullptr) {
        _lastFrame = curFrame;
        std::vector<cv::Point2f> goodFeats;
        // the mask is empty for the first frame, do not use mask
        cv::goodFeaturesToTrack(_lastFrame->GetImage(), goodFeats, FEAT_NUM_PER_IMG, 0.01,
                                MIN_DIST);

        // store found good feasters to map
        int idxCounter = 0;
        for (const auto &p : goodFeats) {
            auto newFeatId = idxCounter++;
            auto newLmId = GenNewLmId();

            ns_veta::Vec2d up = _intri->GetUndistoPixel(ns_veta::Vec2d(p.x, p.y));
            auto feat = Feat(p, cv::Point2f(static_cast<float>(up(0)), static_cast<float>(up(1))));

            _ptsInLast.insert({newFeatId, feat});

            _lmTrackInfo.insert({newLmId, {{_lastFrame, feat}}});

            _featId2lmIdInLast.insert({newFeatId, newLmId});
        }

        // rotation
        _rotations.emplace_back(curFrame->GetTimestamp(), Sophus::SO3d());
    } else {
#if VISUALIZATION
        cv::Mat imgLast = _lastFrame->GetColorImage(), imgCur = curFrame->GetColorImage();
        cv::Mat matImg;
        cv::hconcat(imgLast, imgCur, matImg);
        auto bias = cv::Point2f((float)imgLast.cols, 0.0f);
#endif

        // -----------------------------
        // perform optical flow tracking
        // -----------------------------
        auto [idsLast, rawFeatLast] = ExtractFeatMapAsRawFeatVec(_ptsInLast);
        std::vector<cv::Point2f> rawFeatCur;
        std::vector<uchar> status;
        std::vector<float> errors;
        cv::calcOpticalFlowPyrLK(_lastFrame->GetImage(), curFrame->GetImage(), rawFeatLast,
                                 rawFeatCur, status, errors, cv::Size(21, 21), 3);
        for (int i = 0; i < int(rawFeatCur.size()); i++) {
            // if the feature is tracked but not in border, we set its status as 'fail = 0'
            if (status[i] && !InImageBorder(rawFeatCur[i], curFrame, 1)) {
                status[i] = 0;
            }
        }

        // cur feat id, last feat id
        std::map<int, int> trackIdsCur2Last;
        int idxCounter = 0;
        // id, feature
        std::map<int, Feat> ptsInCur;
        for (int i = 0; i < static_cast<int>(status.size()); ++i) {
            // if tracking succeed, store its info
            if (status[i]) {
                int newFeatId = idxCounter++;
                trackIdsCur2Last.insert({newFeatId, idsLast.at(i)});

                cv::Point2f p = rawFeatCur.at(i);
                ns_veta::Vec2d up = _intri->GetUndistoPixel(ns_veta::Vec2d(p.x, p.y));
                auto feat =
                    Feat(p, cv::Point2f(static_cast<float>(up(0)), static_cast<float>(up(1))));
                ptsInCur.insert({newFeatId, feat});
            }
        }

        // -----------------------------
        // find mask and filter features
        // -----------------------------
        // current id, track count
        std::vector<std::pair<int, int>> trackCount;
        trackCount.reserve(trackIdsCur2Last.size());
        for (const auto &[idCur, idLast] : trackIdsCur2Last) {
            auto lmId = _featId2lmIdInLast.at(idLast);
            trackCount.emplace_back(idCur, _lmTrackInfo.at(lmId).size());
        }
        auto [mask, filteredFeatIdCur] = ComputeMaskAndFilterPts(curFrame, ptsInCur, trackCount);

        for (const auto &[id, count] : trackCount) {
            // this feature has been filtered in 'ComputeMaskAndFilterPts'
            if (filteredFeatIdCur.count(id) == 0) {
                ptsInCur.erase(id);
                trackIdsCur2Last.erase(id);
            }
        }

#if VISUALIZATION
        ShowFeatureTracking(_ptsInLast, ptsInCur, trackIdsCur2Last, matImg, bias,
                            "Optical Flow Tracking", {255, 255, 255});
#endif

        // -------------------------------------------
        // outliers rejection using fundamental matrix
        // -------------------------------------------
        auto undistFeatLast =
            ExtractFeatMapAsUndistoFeatVec(_ptsInLast, ExtractValsAsVec(trackIdsCur2Last));
        auto undistFeatCur =
            ExtractFeatMapAsUndistoFeatVec(ptsInCur, ExtractKeysAsVec(trackIdsCur2Last));
        status = RejectUsingFMat(undistFeatLast.second, undistFeatCur.second);

        for (int i = 0; i < static_cast<int>(status.size()); ++i) {
            // is an outlier, erase this match
            if (!status.at(i)) {
                auto id = undistFeatCur.first.at(i);
                ptsInCur.erase(id);
                trackIdsCur2Last.erase(id);
            }
        }
#if VISUALIZATION
        ShowFeatureTracking(_ptsInLast, ptsInCur, trackIdsCur2Last, matImg, bias,
                            "Feature Matching (FMat Rejection)", {0, 0, 255});
#endif

        // ----------------------------------------------------------
        // perform rotation-only estimation (with outliers rejection)
        // ----------------------------------------------------------
        if (trackIdsCur2Last.size() < 10) {
            return false;
        }

        undistFeatLast =
            ExtractFeatMapAsUndistoFeatVec(_ptsInLast, ExtractValsAsVec(trackIdsCur2Last));
        undistFeatCur =
            ExtractFeatMapAsUndistoFeatVec(ptsInCur, ExtractKeysAsVec(trackIdsCur2Last));

        auto res = RelRotationRecovery(undistFeatLast.second, undistFeatCur.second);

        // solving failed
        if (res.second.empty()) {
            return false;
        }

        // organize results
        Eigen::Matrix3d ROT_CurToLast = res.first;
        Sophus::SO3d SO3_CurToLast(Sophus::makeRotationMatrix(ROT_CurToLast));
        Sophus::SO3d SO3_CurToW = _rotations.back().second * SO3_CurToLast;
        _rotations.emplace_back(curFrame->GetTimestamp(), SO3_CurToW);

        // remove outlier
        std::set<int> inlierSet;
        for (const auto &vecIdx : res.second) {
            inlierSet.insert(vecIdx);
        }

        for (int vecIdx = 0; vecIdx < static_cast<int>(undistFeatCur.first.size()); ++vecIdx) {
            if (inlierSet.count(vecIdx) == 0) {
                // this is an outlier
                auto id = undistFeatCur.first.at(vecIdx);
                ptsInCur.erase(id);
                trackIdsCur2Last.erase(id);
            }
        }

#if VISUALIZATION
        ShowFeatureTracking(_ptsInLast, ptsInCur, trackIdsCur2Last, matImg, bias,
                            "Feature Matching (Rot-only RANSAC Rejection)", {0, 255, 0});
#endif

        // store tracking info
        std::map<int, ns_veta::IndexT> featId2lmIdInLastTmp;

        for (const auto &[idCur, idLast] : trackIdsCur2Last) {
            auto lmId = _featId2lmIdInLast.at(idLast);
            _lmTrackInfo.at(lmId).emplace_back(curFrame, ptsInCur.at(idCur));
            featId2lmIdInLastTmp.insert({idCur, lmId});
        }
        _featId2lmIdInLast = featId2lmIdInLastTmp;

        // --------------------------------------
        // perform incremental feature extraction
        // --------------------------------------
        int featNumToExtract = FEAT_NUM_PER_IMG - static_cast<int>(ptsInCur.size());
        if (featNumToExtract > 0) {
            std::vector<cv::Point2f> newGoodFeats;
            cv::goodFeaturesToTrack(curFrame->GetImage(), newGoodFeats, featNumToExtract, 0.01,
                                    MIN_DIST, mask);
            for (const auto &p : newGoodFeats) {
                auto newFeatId = idxCounter++;
                auto newLmId = GenNewLmId();

                ns_veta::Vec2d up = _intri->GetUndistoPixel(ns_veta::Vec2d(p.x, p.y));
                auto feat =
                    Feat(p, cv::Point2f(static_cast<float>(up(0)), static_cast<float>(up(1))));

                ptsInCur.insert({newFeatId, feat});

                _lmTrackInfo.insert({newLmId, {{curFrame, feat}}});

                _featId2lmIdInLast.insert({newFeatId, newLmId});
            }
        }

        _lastFrame = curFrame;
        _ptsInLast = ptsInCur;
    }
    // #if VISUALIZATION
    ShowCurrentFrame();
    cv::waitKey(1);
// #endif
#undef VISUALIZATION
    return true;
}

void RotOnlyVisualOdometer::ResetWorkspace() {
    _ptsInLast.clear();
    _lmTrackInfo.clear();
    _featId2lmIdInLast.clear();
    _lastFrame = nullptr;
    _rotations.clear();
}

std::pair<cv::Mat, std::set<int>> RotOnlyVisualOdometer::ComputeMaskAndFilterPts(
    const CameraFrame::Ptr &frame,
    const std::map<int, Feat> &featMap,
    // current id, track count
    std::vector<std::pair<int, int>> trackCount) const {
    // sort features based on the tracking count
    std::sort(trackCount.begin(), trackCount.end(),
              [](const auto &p1, const auto &p2) { return std::get<1>(p1) > std::get<1>(p2); });

    int col = frame->GetImage().cols, row = frame->GetImage().rows;
    cv::Mat mask = cv::Mat(row, col, CV_8UC1, cv::Scalar(255));

    std::set<int> filteredFeatId;
    for (const auto &[id, count] : trackCount) {
        auto pt = featMap.at(id).raw;
        // if this not has been occupied
        if (mask.at<uchar>(pt) == 255) {
            // draw mask
            cv::circle(mask, pt, MIN_DIST, 0, -1);
            filteredFeatId.insert(id);
        }
    }

    return {mask, filteredFeatId};
}

bool RotOnlyVisualOdometer::InImageBorder(const cv::Point2f &pt,
                                          const CameraFrame::Ptr &frame,
                                          int borderSize) {
    int col = frame->GetImage().cols, row = frame->GetImage().rows;
    int imgX = cvRound(pt.x), imgY = cvRound(pt.y);
    return borderSize <= imgX && imgX < col - borderSize && borderSize <= imgY &&
           imgY < row - borderSize;
}

void RotOnlyVisualOdometer::ShowCurrentFrame() const {
    cv::Mat img = _lastFrame->GetImage().clone();
    cv::cvtColor(img, img, cv::COLOR_GRAY2BGR);

    for (const auto &[id, feat] : _ptsInLast) {
        const auto &pt = feat.raw;
        const auto &upt = feat.undistorted;

        // connect between point and its undistorted one
        DrawLineOnCVMat(img, pt, upt, cv::Scalar(255, 255, 255));
        // undistorted point
        DrawKeypointOnCVMat(img, upt, false, cv::Scalar(255, 255, 255));
        // raw point
        DrawKeypointOnCVMat(img, pt);

        // text: track count
        int count = static_cast<int>(_lmTrackInfo.at(_featId2lmIdInLast.at(id)).size());
        PutTextOnCVMat(img, std::to_string(count), pt);
    }
    // auto [K, D] = CalibParamManager().INTRI.ObtainKDMatForUndisto(_intri);
    // cv::Mat undistImg;
    // cv::undistort(_lastFrame->GetImage(), undistImg, K, D);

    const static std::string winName = "Feature Tracking";
    cv::imshow(winName, img);
    // static int count = 0;
    // cv::imwrite(Configor::DataStream::DebugPath + "/tracking" + std::to_string(count++) + ".png",
    // img); const static std::string unWinName = "Undisto Image"; cv::namedWindow(unWinName,
    // cv::WindowFlags::WINDOW_NORMAL); cv::imshow(unWinName, undistImg);
}

void RotOnlyVisualOdometer::ShowFeatureTracking(const std::map<int, Feat> &ptsInLast,
                                                const std::map<int, Feat> &ptsInCur,
                                                const std::map<int, int> &matches,
                                                cv::Mat matImg,
                                                const cv::Point2f &bias,
                                                const std::string &winName,
                                                const cv::Scalar &color) {
    for (const auto &[idCur, idLast] : matches) {
        cv::Point2f pt1 = ptsInLast.at(idLast).raw;
        cv::Point2f pt2 = ptsInCur.at(idCur).raw + bias;
        cv::drawMarker(matImg, pt1, color, cv::MarkerTypes::MARKER_SQUARE, 10, 1);
        cv::drawMarker(matImg, pt1, color, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
        cv::drawMarker(matImg, pt2, color, cv::MarkerTypes::MARKER_SQUARE, 10, 1);
        cv::drawMarker(matImg, pt2, color, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
        cv::line(matImg, pt1, pt2, color, 1);
    }

    cv::imshow(winName, matImg);
    cv::waitKey();
    // static int count = 0;
    // cv::imwrite(Configor::DataStream::DebugPath + "/matching" + std::to_string(count++) + ".png",
    // matImg);
}

std::vector<uchar> RotOnlyVisualOdometer::RejectUsingFMat(
    const std::vector<cv::Point2f> &undistPtsInLast,
    const std::vector<cv::Point2f> &undistPtsInCur) {
    std::vector<uchar> status;
    cv::findFundamentalMat(undistPtsInLast, undistPtsInCur, cv::FM_RANSAC, 1.0, 0.99, status);
    return status;
}

std::pair<opengv::rotation_t, std::vector<int>> RotOnlyVisualOdometer::RelRotationRecovery(
    const std::vector<cv::Point2f> &ptsUndisto1, const std::vector<cv::Point2f> &ptsUndisto2) {
    assert(ptsUndisto1.size() == ptsUndisto2.size());
    opengv::bearingVectors_t bearingVectors1 = ComputeBeringVec(ptsUndisto1);
    opengv::bearingVectors_t bearingVectors2 = ComputeBeringVec(ptsUndisto2);
    // create a central relative adapter
    opengv::relative_pose::CentralRelativeAdapter adapter(bearingVectors1, bearingVectors2);

    // using all correspondences
    // opengv::rotation_t rotation = opengv::relative_pose::rotationOnly(adapter);

    // Create a RotationOnlySacProblem and Ransac
    opengv::sac::Ransac<opengv::sac_problems::relative_pose::RotationOnlySacProblem> ransac;
    std::shared_ptr<opengv::sac_problems::relative_pose::RotationOnlySacProblem> probPtr(
        new opengv::sac_problems::relative_pose::RotationOnlySacProblem(adapter));
    ransac.sac_model_ = probPtr;
    ransac.threshold_ = _intri->ImagePlaneToCameraPlaneError(1.0);
    ransac.max_iterations_ = 50;

    auto res = ransac.computeModel(0);

    if (!res) {
        return {};
    }

    opengv::rotation_t rotation = ransac.model_coefficients_;

    // spdlog::info("ransac iterations: {}, inliers rate: {}/{}", ransac.iterations_,
    //              ransac.inliers_.size(), ptsUndisto1.size());

    return {rotation, ransac.inliers_};
}

opengv::bearingVectors_t RotOnlyVisualOdometer::ComputeBeringVec(
    const std::vector<cv::Point2f> &ptsUndist) {
    opengv::bearingVectors_t bearingVec(ptsUndist.size());
    for (int i = 0; i < static_cast<int>(ptsUndist.size()); ++i) {
        const auto &p = ptsUndist.at(i);
        ns_veta::Vec2d pCam = _intri->ImgToCam(ns_veta::Vec2d(p.x, p.y));
        bearingVec.at(i) = ns_veta::Vec3d(pCam(0), pCam(1), 1.0).normalized();
    }
    return bearingVec;
}

const std::vector<std::pair<double, Sophus::SO3d>> &RotOnlyVisualOdometer::GetRotations() const {
    return _rotations;
}

RotOnlyVisualOdometer::~RotOnlyVisualOdometer() { cv::destroyAllWindows(); }

ns_veta::IndexT RotOnlyVisualOdometer::GenNewLmId() {
    static ns_veta::IndexT id = 0;
    return ++id;
}

std::pair<std::vector<int>, std::vector<cv::Point2f>>
RotOnlyVisualOdometer::ExtractFeatMapAsRawFeatVec(const std::map<int, Feat> &featMap,
                                                  const std::vector<int> &desiredIds) {
    std::vector<int> ids;
    std::vector<cv::Point2f> feats;
    if (desiredIds.empty()) {
        ids.reserve(featMap.size());
        feats.reserve(featMap.size());
        for (const auto &[id, pt] : featMap) {
            ids.push_back(id);
            // distorted point
            feats.push_back(pt.raw);
        }
    } else {
        ids.reserve(desiredIds.size());
        feats.reserve(desiredIds.size());
        for (const auto &id : desiredIds) {
            if (auto iter = featMap.find(id); iter != featMap.cend()) {
                ids.push_back(id);
                // distorted point
                feats.push_back(iter->second.raw);
            }
        }
    }
    return {ids, feats};
}

std::pair<std::vector<int>, std::vector<cv::Point2f>>
RotOnlyVisualOdometer::ExtractFeatMapAsUndistoFeatVec(const std::map<int, Feat> &featMap,
                                                      const std::vector<int> &desiredIds) {
    std::vector<int> ids;
    std::vector<cv::Point2f> feats;
    if (desiredIds.empty()) {
        ids.reserve(featMap.size());
        feats.reserve(featMap.size());
        for (const auto &[id, pt] : featMap) {
            ids.push_back(id);
            // distorted point
            feats.push_back(pt.undistorted);
        }
    } else {
        ids.reserve(desiredIds.size());
        feats.reserve(desiredIds.size());
        for (const auto &id : desiredIds) {
            if (auto iter = featMap.find(id); iter != featMap.cend()) {
                ids.push_back(id);
                // distorted point
                feats.push_back(iter->second.undistorted);
            }
        }
    }
    return {ids, feats};
}
const std::map<ns_veta::IndexT, std::list<std::pair<CameraFramePtr, RotOnlyVisualOdometer::Feat>>> &
RotOnlyVisualOdometer::GetLmTrackInfo() const {
    return _lmTrackInfo;
}

void RotOnlyVisualOdometer::ShowLmTrackInfo() const {
    for (const auto &[lmId, trackList] : _lmTrackInfo) {
        spdlog::info("landmark id: {}, track count: {}", lmId, trackList.size());
        if (trackList.size() == 1) {
            continue;
        }

        for (const auto &[frame, feat] : trackList) {
            cv::Mat img = frame->GetColorImage().clone();

            const auto &pt = feat.raw;
            const auto &upt = feat.undistorted;

            // connect between point and its undistorted one
            DrawLineOnCVMat(img, pt, upt, cv::Scalar(255, 255, 255));
            // undistorted point
            DrawKeypointOnCVMat(img, upt, false, cv::Scalar(255, 255, 255));
            // raw point
            DrawKeypointOnCVMat(img, pt);

            PutTextOnCVMat(img, std::to_string(trackList.size()), pt);

            cv::imshow(std::to_string(frame->GetId()), img);
        }
        cv::waitKey(0);
    }
}

RotOnlyVisualOdometer::Feat::Feat(cv::Point2f raw, cv::Point2f undistorted)
    : raw(std::move(raw)),
      undistorted(std::move(undistorted)) {}
}  // namespace ns_ikalibr