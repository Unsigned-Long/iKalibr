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

RotOnlyVisualOdometer::RotOnlyVisualOdometer(FeatureTracking::Ptr featTracking,
                                             ns_veta::PinholeIntrinsic::Ptr intri)
    : _featTracking(std::move(featTracking)),
      _intri(std::move(intri)),
      _trackFeatLast(nullptr) {}

RotOnlyVisualOdometer::Ptr RotOnlyVisualOdometer::Create(
    const FeatureTracking::Ptr &featTracking, const ns_veta::PinholeIntrinsic::Ptr &intri) {
    return std::make_shared<RotOnlyVisualOdometer>(featTracking, intri);
}

bool RotOnlyVisualOdometer::GrabFrame(const CameraFrame::Ptr &curFrame,
                                      const std::optional<Sophus::SO3d> &SO3_LastToCur) {
#define VISUALIZATION 0
    // grab the current image
    const auto trackedFeats = _featTracking->GrabImageFrame(curFrame, SO3_LastToCur);

    if (_trackFeatLast == nullptr) {
        // the first tracked frame
        for (const auto &[idCur, feat] : trackedFeats->featCur) {
            auto newLmId = GenNewLmId();
            // create landmarks
            _lmTrackInfo.insert({newLmId, {{curFrame, feat}}});
            // store tracking information, feat id in last image --> landmark id
            _featId2lmIdInLast.insert({idCur, newLmId});
        }
        _rotations.emplace_back(curFrame->GetTimestamp(), Sophus::SO3d());
        _trackFeatLast = trackedFeats;
        return true;
    }
#if VISUALIZATION
    // cv::imshow("Raw Feature Tracking", trackedFeats->DrawMatches());
    // cv::waitKey(0);
    auto matchBackup = std::make_shared<FeatureTracking::TrackedFeaturePack>(*trackedFeats);
#endif
    // -------------------------------------------
    // outliers rejection using fundamental matrix
    // -------------------------------------------
    // spdlog::info("perform outliers rejection using fundamental matrix...");
    if (trackedFeats->featMatchLast2Cur.size() < 15) {
        this->ResetWorkspace();
        return false;
    }

    auto undistFeatLast = ExtractFeatMapAsUndistoFeatVec(
        trackedFeats->featLast, ExtractKeysAsVec(trackedFeats->featMatchLast2Cur));
    auto undistFeatCur = ExtractFeatMapAsUndistoFeatVec(
        trackedFeats->featCur, ExtractValsAsVec(trackedFeats->featMatchLast2Cur));

    auto status = RejectUsingFMat(undistFeatLast.second, undistFeatCur.second);

    for (int i = 0; i < static_cast<int>(status.size()); ++i) {
        // is an outlier, erase this match
        if (!status.at(i)) {
            auto idLast = undistFeatLast.first.at(i);
            auto idCur = undistFeatCur.first.at(i);
            trackedFeats->featLast.erase(idLast);
            trackedFeats->featCur.erase(idCur);
            trackedFeats->featMatchLast2Cur.erase(idLast);
        }
    }

#if VISUALIZATION
    // cv::imshow("FMat Rejection", trackedFeats->DrawMatches(matchBackup));
    // cv::waitKey(0);
    matchBackup = std::make_shared<FeatureTracking::TrackedFeaturePack>(*trackedFeats);
#endif
    // ----------------------------------------------------------
    // perform rotation-only estimation (with outliers rejection)
    // ----------------------------------------------------------
    // spdlog::info("perform rotation-only estimation (with outliers rejection)...");
    if (trackedFeats->featMatchLast2Cur.size() < 10) {
        this->ResetWorkspace();
        return false;
    }

    undistFeatLast = ExtractFeatMapAsUndistoFeatVec(
        trackedFeats->featLast, ExtractKeysAsVec(trackedFeats->featMatchLast2Cur));
    undistFeatCur = ExtractFeatMapAsUndistoFeatVec(
        trackedFeats->featCur, ExtractValsAsVec(trackedFeats->featMatchLast2Cur));

    auto res = RelRotationRecovery(undistFeatLast.second, undistFeatCur.second);

    // solving failed
    if (res.second.empty()) {
        this->ResetWorkspace();
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

    for (int i = 0; i < static_cast<int>(undistFeatCur.first.size()); ++i) {
        if (inlierSet.count(i) == 0) {
            // this is an outlier
            auto idLast = undistFeatLast.first.at(i);
            auto idCur = undistFeatCur.first.at(i);
            trackedFeats->featLast.erase(idLast);
            trackedFeats->featCur.erase(idCur);
            trackedFeats->featMatchLast2Cur.erase(idLast);
        }
    }

#if VISUALIZATION
    cv::imshow("Rot-only RANSAC Rejection", trackedFeats->DrawMatches(matchBackup));
    cv::waitKey(0);
    matchBackup = nullptr;
#endif

    // store tracking info
    // spdlog::info("store tracking information...");
    std::map<int, ns_veta::IndexT> featId2lmIdInLastTmp;
    for (const auto &[idLast, idCur] : trackedFeats->featMatchLast2Cur) {
        if (auto lmIdIter = _featId2lmIdInLast.find(idLast);
            lmIdIter == _featId2lmIdInLast.cend()) {
            // new landmark
            auto newLmId = GenNewLmId();
            _lmTrackInfo.insert({newLmId, {{curFrame, trackedFeats->featCur.at(idCur)}}});
            featId2lmIdInLastTmp.insert({idCur, newLmId});
        } else {
            // old landmark
            auto lmId = lmIdIter->second;
            _lmTrackInfo.at(lmId).emplace_back(curFrame, trackedFeats->featCur.at(idCur));
            featId2lmIdInLastTmp.insert({idCur, lmId});
        }
    }
    _featId2lmIdInLast = featId2lmIdInLastTmp;
    _trackFeatLast = trackedFeats;

    // spdlog::info("show tracked features on the image...");
    ShowCurrentFrame();
    cv::waitKey(1);
#undef VISUALIZATION
    return true;
}

void RotOnlyVisualOdometer::ResetWorkspace() {
    _lmTrackInfo.clear();
    _featId2lmIdInLast.clear();
    _rotations.clear();
}

void RotOnlyVisualOdometer::ShowCurrentFrame() const {
    cv::Mat img = _trackFeatLast->imgCur->GetImage().clone();
    cv::cvtColor(img, img, cv::COLOR_GRAY2BGR);

    for (const auto &[id, feat] : _trackFeatLast->featCur) {
        int count = static_cast<int>(_lmTrackInfo.at(_featId2lmIdInLast.at(id)).size());
        if (count < 2) {
            continue;
        }
        const auto &pt = feat.raw;
        const auto &upt = feat.undistorted;

        // connect between point and its undistorted one
        DrawLineOnCVMat(img, pt, upt, cv::Scalar(255, 255, 255));
        // undistorted point
        DrawKeypointOnCVMat(img, upt, false, cv::Scalar(255, 255, 255));
        // raw point
        DrawKeypointOnCVMat(img, pt);

        // text: track count
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

std::vector<uchar> RotOnlyVisualOdometer::RejectUsingFMat(
    const std::vector<cv::Point2f> &undistPtsInLast,
    const std::vector<cv::Point2f> &undistPtsInCur) {
    std::vector<uchar> status;
    cv::findFundamentalMat(undistPtsInLast, undistPtsInCur, cv::FM_RANSAC, 1.0, 0.99, status);
    return status;
}

std::pair<opengv::rotation_t, std::vector<int>> RotOnlyVisualOdometer::RelRotationRecovery(
    const std::vector<cv::Point2f> &ptsUndisto1,
    const std::vector<cv::Point2f> &ptsUndisto2) const {
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
    const std::vector<cv::Point2f> &ptsUndist) const {
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
RotOnlyVisualOdometer::ExtractFeatMapAsRawFeatVec(const std::map<int, Feature> &featMap,
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
RotOnlyVisualOdometer::ExtractFeatMapAsUndistoFeatVec(const std::map<int, Feature> &featMap,
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

const std::map<ns_veta::IndexT, std::list<std::pair<CameraFramePtr, Feature>>> &
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
}  // namespace ns_ikalibr