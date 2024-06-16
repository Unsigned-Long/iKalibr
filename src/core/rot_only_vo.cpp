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

#include "opencv2/highgui.hpp"
#include "opencv2/video/tracking.hpp"
#include "opencv2/calib3d.hpp"

#include "opengv/sac/Ransac.hpp"
#include "opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp"
#include "opengv/relative_pose/CentralRelativeAdapter.hpp"

_3_

namespace ns_ikalibr {

    ORBFeatExtMatConfig::ORBFeatExtMatConfig() = default;

    RotOnlyVisualOdometer::Ptr
    RotOnlyVisualOdometer::Create(int featNumPerImg, int minDist, const ns_veta::PinholeIntrinsic::Ptr &intri) {
        return std::make_shared<RotOnlyVisualOdometer>(featNumPerImg, minDist, intri);
    }

    bool RotOnlyVisualOdometer::GrabFrame(const CameraFrame::Ptr &curFrame) {
#define VISUALIZATION 0

        // first curFrame
        if (_lastFrame == nullptr) {
            _lastFrame = curFrame;
            // the mask is empty for the first frame do not use mask
            cv::goodFeaturesToTrack(_lastFrame->GetImage(), _ptsInLast, FEAT_NUM_PER_IMG, 0.01, MIN_DIST);
            _ptsUndistInLast = UndistortedPoints(_ptsInLast);
            _ptsTrackCount = std::vector<int>(_ptsInLast.size(), 1);
            // rotation
            _rotations.emplace_back(curFrame->GetTimestamp(), Sophus::SO3d());
        } else {
            // perform optical flow track
            std::vector<cv::Point2f> ptsInCur;
            std::vector<uchar> status;
            std::vector<float> errors;
            cv::calcOpticalFlowPyrLK(
                    _lastFrame->GetImage(), curFrame->GetImage(), _ptsInLast, ptsInCur,
                    status, errors, cv::Size(21, 21), 3
            );
            for (int i = 0; i < int(ptsInCur.size()); i++) {
                // if the feature is tracked but not in border, we set its status as 'fail = 0'
                if (status[i] && !InImageBorder(ptsInCur[i], curFrame, 1)) { status[i] = 0; }
            }
            // if track is failed, remove feature info
            ReduceVector(_ptsInLast, status);
            ReduceVector(_ptsTrackCount, status);
            for (auto &count: _ptsTrackCount) { ++count; }
            ReduceVector(_ptsUndistInLast, status);
            ReduceVector(ptsInCur, status);

            // find mask and filter features
            auto [mask, filteredPts, filteredTrackCount, filteredOldIds] =
                    ComputeMaskAndFilterPts(curFrame, ptsInCur, _ptsTrackCount);

            ptsInCur = filteredPts;
            _ptsTrackCount = filteredTrackCount;
            _ptsInLast = FindElements(_ptsInLast, filteredOldIds);
            _ptsUndistInLast = FindElements(_ptsUndistInLast, filteredOldIds);
            std::vector<cv::Point2f> ptsUndistInCur = UndistortedPoints(ptsInCur);

            // outliers rejection
            status = RejectUsingFMat(_ptsUndistInLast, ptsUndistInCur);
            std::vector<int> inliers;
            inliers.reserve(status.size());
            for (int i = 0; i < static_cast<int>(status.size()); ++i) {
                if (status.at(i)) { inliers.push_back(i); }
            }
#if VISUALIZATION
            ShowFeatureTracking(_lastFrame, _ptsInLast, curFrame, ptsInCur, inliers,
                                "Feature Matching (FMat Rejection)");
#endif
            ReduceVector(_ptsInLast, status);
            ReduceVector(ptsInCur, status);
            ReduceVector(_ptsTrackCount, status);
            ReduceVector(_ptsUndistInLast, status);
            ReduceVector(ptsUndistInCur, status);

            // perform rotation-only estimation
            if (_ptsUndistInLast.size() < 10) {
                // at least two correspondences are required
                _ptsInLast.clear(), _ptsTrackCount.clear(), _ptsUndistInLast.clear(), _lastFrame = nullptr, _rotations.clear();
                return false;
            }
            auto res = RelRotationRecovery(_ptsUndistInLast, ptsUndistInCur);
            // solving failed
            if (res.second.empty()) {
                _ptsInLast.clear(), _ptsTrackCount.clear(), _ptsUndistInLast.clear(), _lastFrame = nullptr, _rotations.clear();
                return false;
            }
            Eigen::Matrix3d ROT_CurToLast = res.first;
            inliers = res.second;
            Sophus::SO3d SO3_CurToLast(Sophus::makeRotationMatrix(ROT_CurToLast));
            Sophus::SO3d SO3_CurToW = _rotations.back().second * SO3_CurToLast;
            _rotations.emplace_back(curFrame->GetTimestamp(), SO3_CurToW);

            // reset outlier track count
            std::set<int> inliersSet(inliers.begin(), inliers.end());
            for (int i = 0; i < static_cast<int>(_ptsTrackCount.size()); ++i) {
                if (inliersSet.find(i) == inliersSet.cend()) { _ptsTrackCount.at(i) = 1; }
            }

#if VISUALIZATION
            ShowFeatureTracking(_lastFrame, _ptsInLast, curFrame, ptsInCur, inliers,
                                "Feature Matching (RANSAC Rejection)");
#endif

            // perform incremental feature extraction
            int featNumToExtract = FEAT_NUM_PER_IMG - static_cast<int>(ptsInCur.size());
            if (featNumToExtract > 0) {
                std::vector<cv::Point2f> newPtsInCur;
                cv::goodFeaturesToTrack(curFrame->GetImage(), newPtsInCur, featNumToExtract, 0.01, MIN_DIST, mask);
                // undisto
                std::vector<cv::Point2f> newPtsUndistInCur = UndistortedPoints(newPtsInCur);
                // append new features
                ptsInCur.reserve(ptsInCur.size() + newPtsInCur.size());
                _ptsTrackCount.reserve(_ptsTrackCount.size() + newPtsInCur.size());
                ptsUndistInCur.reserve(ptsUndistInCur.size() + newPtsUndistInCur.size());
                for (int i = 0; i < static_cast<int>(newPtsInCur.size()); ++i) {
                    ptsInCur.push_back(newPtsInCur.at(i)), _ptsTrackCount.push_back(1);
                    ptsUndistInCur.push_back(newPtsUndistInCur.at(i));
                }
            }
            _lastFrame = curFrame;
            _ptsInLast = ptsInCur;
            _ptsUndistInLast = ptsUndistInCur;

        }
//#if VISUALIZATION
        ShowCurrentFrame();
        cv::waitKey(1);
//#endif
#undef VISUALIZATION
        return true;
    }

    bool RotOnlyVisualOdometer::InImageBorder(const cv::Point2f &pt, const CameraFrame::Ptr &frame, int borderSize) {
        int col = frame->GetImage().cols, row = frame->GetImage().rows;
        int imgX = cvRound(pt.x), imgY = cvRound(pt.y);
        return borderSize <= imgX && imgX < col - borderSize && borderSize <= imgY && imgY < row - borderSize;
    }

    std::vector<cv::Point2f> RotOnlyVisualOdometer::UndistortedPoints(const std::vector<cv::Point2f> &pts) {
        std::vector<cv::Point2f> undistPts(pts.size());
        for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
            const auto &p = pts.at(i);
            ns_veta::Vec2d up = _intri->GetUndistoPixel(ns_veta::Vec2d(p.x, p.y));
            undistPts.at(i).x = static_cast<float>(up(0)), undistPts.at(i).y = static_cast<float>(up(1));
        }
        return undistPts;
    }

    std::tuple<cv::Mat, std::vector<cv::Point2f>, std::vector<int>, std::vector<int>>
    RotOnlyVisualOdometer::ComputeMaskAndFilterPts(const CameraFrame::Ptr &frame,
                                                   const std::vector<cv::Point2f> &pts,
                                                   const std::vector<int> &trackCount) const {
        assert(pts.size() == trackCount.size());
        int size = static_cast<int>(pts.size());

        int col = frame->GetImage().cols, row = frame->GetImage().rows;
        cv::Mat mask = cv::Mat(row, col, CV_8UC1, cv::Scalar(255));

        // point, track count, source id
        std::vector<std::tuple<cv::Point2f, int, int>> ptsCount(pts.size());
        for (int i = 0; i < size; ++i) {
            ptsCount.at(i) = std::make_tuple(pts.at(i), trackCount.at(i), i);
        }
        // sort based on the track count
        std::sort(ptsCount.begin(), ptsCount.end(), [](const auto &p1, const auto &p2) {
            return std::get<1>(p1) > std::get<1>(p2);
        });
        std::vector<cv::Point2f> filteredPts;
        std::vector<int> filteredTrackCount, filteredOldIds;
        filteredPts.reserve(size), filteredTrackCount.reserve(size), filteredOldIds.reserve(size);

        for (auto &[pt, count, id]: ptsCount) {
            if (mask.at<uchar>(pt) == 255) {
                filteredPts.push_back(pt);
                filteredTrackCount.push_back(count);
                filteredOldIds.push_back(id);
                // draw mask
                cv::circle(mask, pt, MIN_DIST, 0, -1);
            }
        }
        return {mask, filteredPts, filteredTrackCount, filteredOldIds};
    }

    void RotOnlyVisualOdometer::ShowCurrentFrame() const {
        cv::Mat img = _lastFrame->GetImage().clone();
        cv::cvtColor(img, img, cv::COLOR_GRAY2BGR);
        int size = static_cast<int>(_ptsTrackCount.size());
        for (int i = 0; i < size; ++i) {
            const auto &pt = _ptsInLast.at(i);
            const auto &upt = _ptsUndistInLast.at(i);
            // square
            cv::drawMarker(img, pt, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 10, 1);
            // key point
            cv::drawMarker(img, pt, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
            cv::drawMarker(img, upt, cv::Scalar(255, 255, 255), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
            cv::line(img, pt, upt, cv::Scalar(255, 255, 255), 1);
            // text: track count
            cv::putText(img, std::to_string(_ptsTrackCount.at(i)), cv::Point2f(pt.x + 10, pt.y),
                        cv::HersheyFonts::FONT_HERSHEY_PLAIN, 1.0, cv::Scalar(255, 0, 0), 2);
        }
        // auto [K, D] = CalibParamManager().INTRI.ObtainKDMatForUndisto(_intri);
        // cv::Mat undistImg;
        // cv::undistort(_lastFrame->GetImage(), undistImg, K, D);

        const static std::string winName = "Feature Tracking";
        cv::imshow(winName, img);
        // static int count = 0;
        // cv::imwrite(Configor::DataStream::DebugPath + "/tracking" + std::to_string(count++) + ".png", img);
        // const static std::string unWinName = "Undisto Image";
        // cv::namedWindow(unWinName, cv::WindowFlags::WINDOW_NORMAL);
        // cv::imshow(unWinName, undistImg);
    }

    void RotOnlyVisualOdometer::ShowFeatureTracking(const CameraFrame::Ptr &lastFrame,
                                                    const std::vector<cv::Point2f> &ptsInLast,
                                                    const CameraFrame::Ptr &curFrame,
                                                    const std::vector<cv::Point2f> &ptsInCur,
                                                    const std::vector<int> &inliers,
                                                    const std::string &winName) {
        assert(ptsInLast.size() == ptsInCur.size());
        std::vector<cv::KeyPoint> kpsInLast(ptsInLast.size()), kpsInCur(ptsInCur.size());
        std::vector<cv::KeyPoint> inKpsInLast(inliers.size()), inKpsInCur(inliers.size());
        std::vector<cv::DMatch> matches(ptsInLast.size()), inMatches(inliers.size());
        for (int i = 0; i < static_cast<int>(ptsInLast.size()); ++i) {
            kpsInLast.at(i).pt = ptsInLast.at(i);
            matches.at(i).queryIdx = matches.at(i).trainIdx = i;
        }
        for (int i = 0; i < static_cast<int>(ptsInCur.size()); ++i) {
            kpsInCur.at(i).pt = ptsInCur.at(i);
        }
        for (int i = 0; i < static_cast<int>(inliers.size()); ++i) {
            int idx = inliers.at(i);
            inKpsInLast.at(i).pt = ptsInLast.at(idx);
            inKpsInCur.at(i).pt = ptsInCur.at(idx);
            inMatches.at(i).queryIdx = inMatches.at(i).trainIdx = i;
        }

        cv::Mat imgLast, imgCur;
        cv::cvtColor(lastFrame->GetImage(), imgLast, cv::COLOR_GRAY2BGR);
        cv::cvtColor(curFrame->GetImage(), imgCur, cv::COLOR_GRAY2BGR);
        cv::Mat matImg;
        cv::hconcat(imgLast, imgCur, matImg);
        auto bias = cv::Point2f((float) imgLast.cols, 0.0f);

        for (const auto &match: matches) {
            cv::Point2f pt1 = kpsInLast[match.queryIdx].pt;
            cv::Point2f pt2 = kpsInCur[match.trainIdx].pt + bias;
            cv::drawMarker(matImg, pt1, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 10, 1);
            cv::drawMarker(matImg, pt1, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
            cv::drawMarker(matImg, pt2, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 10, 1);
            cv::drawMarker(matImg, pt2, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
            cv::line(matImg, pt1, pt2, cv::Scalar(0, 0, 255), 1);
        }

        for (const auto &match: inMatches) {
            cv::Point2f pt1 = inKpsInLast[match.queryIdx].pt;
            cv::Point2f pt2 = inKpsInCur[match.trainIdx].pt + bias;
            cv::drawMarker(matImg, pt1, cv::Scalar(0, 255, 0), cv::MarkerTypes::MARKER_SQUARE, 10, 1);
            cv::drawMarker(matImg, pt1, cv::Scalar(0, 255, 0), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
            cv::drawMarker(matImg, pt2, cv::Scalar(0, 255, 0), cv::MarkerTypes::MARKER_SQUARE, 10, 1);
            cv::drawMarker(matImg, pt2, cv::Scalar(0, 255, 0), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
            cv::line(matImg, pt1, pt2, cv::Scalar(0, 255, 0), 1);
        }

        cv::imshow(winName, matImg);
        // static int count = 0;
        // cv::imwrite(Configor::DataStream::DebugPath + "/matching" + std::to_string(count++) + ".png", matImg);
    }

    std::vector<uchar> RotOnlyVisualOdometer::RejectUsingFMat(const std::vector<cv::Point2f> &undistPtsInLast,
                                                              const std::vector<cv::Point2f> &undistPtsInCur) {
        std::vector<uchar> status;
        cv::findFundamentalMat(undistPtsInLast, undistPtsInCur, cv::FM_RANSAC, 1.0, 0.99, status);
        return status;
    }

    std::pair<opengv::rotation_t, std::vector<int>>
    RotOnlyVisualOdometer::RelRotationRecovery(const std::vector<cv::Point2f> &ptsUndisto1,
                                               const std::vector<cv::Point2f> &ptsUndisto2) {
        opengv::bearingVectors_t bearingVectors1 = ComputeBeringVec(ptsUndisto1);
        opengv::bearingVectors_t bearingVectors2 = ComputeBeringVec(ptsUndisto2);
        //create a central relative adapter
        opengv::relative_pose::CentralRelativeAdapter adapter(bearingVectors1, bearingVectors2);

        // using all correspondences
        // opengv::rotation_t rotation = opengv::relative_pose::rotationOnly(adapter);

        // Create a RotationOnlySacProblem and Ransac
        opengv::sac::Ransac<opengv::sac_problems::relative_pose::RotationOnlySacProblem> ransac;
        std::shared_ptr<opengv::sac_problems::relative_pose::RotationOnlySacProblem> probPtr(
                new opengv::sac_problems::relative_pose::RotationOnlySacProblem(adapter)
        );
        ransac.sac_model_ = probPtr;
        ransac.threshold_ = _intri->ImagePlaneToCameraPlaneError(1.0);
        ransac.max_iterations_ = 50;

        auto res = ransac.computeModel(0);

        if (!res) { return {}; }

        opengv::rotation_t rotation = ransac.model_coefficients_;

        // spdlog::info("ransac iterations: {}, inliers rate: {}/{}", ransac.iterations_,
        //              ransac.inliers_.size(), ptsUndisto1.size());

        return {rotation, ransac.inliers_};
    }

    opengv::bearingVectors_t RotOnlyVisualOdometer::ComputeBeringVec(const std::vector<cv::Point2f> &ptsUndist) {
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

    RotOnlyVisualOdometer::~RotOnlyVisualOdometer() {
        cv::destroyAllWindows();
    }
}