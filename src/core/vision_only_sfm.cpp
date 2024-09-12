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

#include "core/vision_only_sfm.h"
#include "util/status.hpp"
#include "util/utils_tpl.hpp"

#include "calib/calib_param_manager.h"
#include "calib/estimator.h"
#include "sensor/camera.h"
#include "viewer/viewer.h"

#include "tiny-viewer/entity/cube.h"
#include "tiny-viewer/entity/line.h"

#include "boost/geometry.hpp"

#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/calib3d/calib3d.hpp"

#include "opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp"
#include "opengv/sac_problems/relative_pose/TranslationOnlySacProblem.hpp"
#include "opengv/relative_pose/CentralRelativeAdapter.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
#define USE_OPENCV_REL_POSE 1
#define USE_OPENCV_TRIANGULATE 1

VisionOnlySfM::VisionOnlySfM(std::string topic,
                             const std::vector<CameraFrame::Ptr> &frames,
                             CalibParamManager::Ptr parMagr,
                             const So3SplineType &so3Spline,
                             Viewer::Ptr viewer)
    : _topic(std::move(topic)),
      _frames(frames),
      _parMagr(std::move(parMagr)),
      _so3Spline(so3Spline),
      _lmLabeler(0),
      _viewer(std::move(viewer)) {}

VisionOnlySfM::Ptr VisionOnlySfM::Create(const std::string &topic,
                                         const std::vector<CameraFrame::Ptr> &frames,
                                         const CalibParamManager::Ptr &parMagr,
                                         const VisionOnlySfM::So3SplineType &so3Spline,
                                         const Viewer::Ptr &viewer) {
    return std::make_shared<VisionOnlySfM>(topic, frames, parMagr, so3Spline, viewer);
}

bool VisionOnlySfM::PreProcess() {
    _veta = ns_veta::Veta::Create();

    constexpr ns_veta::IndexT intriIdx = 0;

    // -----------------
    // create intrinsics
    // -----------------
    _intri = std::make_shared<ns_veta::PinholeIntrinsic>(*_parMagr->INTRI.Camera.at(_topic));
    _veta->intrinsics.insert({intriIdx, _intri});

    // ----------------------
    // create views and poses
    // ----------------------
    for (const auto &frame : _frames) {
        const ns_veta::IndexT viewId = frame->GetId();
        const ns_veta::IndexT poseId = viewId;
        _frameBackup.insert({viewId, frame});

        auto pose = ComputeCamRotations(frame);
        if (pose) {
            // rough rotation-only pose
            _veta->poses.insert({poseId, *pose});

            // view
            auto view = ns_veta::View::Create(frame->GetTimestamp(), viewId, intriIdx, poseId,
                                              _intri->imgWidth, _intri->imgHeight);
            _veta->views.insert({viewId, view});
        }
    }

    CreateViewCubes();
    _viewer->AddEntity(_viewCubes, Viewer::VIEW_ASSOCIATION);

    // ------------------
    // feature extraction
    // ------------------
    spdlog::info("start extracting features for each image, this would cost some time...");
    std::map<ns_veta::IndexT, FeaturePack> featMap;
#pragma omp parallel for num_threads(omp_get_max_threads()) default(none) shared(featMap, _intri)
    for (int i = 0; i < static_cast<int>(_frames.size()); ++i) {
        // use detector to detect features
        std::vector<cv::KeyPoint> kps;
        cv::Mat descriptor;
        cv::AKAZE::create()->detectAndCompute(_frames.at(i)->GetImage(), cv::noArray(), kps,
                                              descriptor);
        std::vector<ns_veta::IndexT> index(kps.size());
        std::vector<ns_veta::Vec2d> kpsUndisto(kps.size());
        for (int j = 0; j < static_cast<int>(kpsUndisto.size()); ++j) {
            const auto &kp = kps.at(j).pt;
            kpsUndisto.at(j) = _intri->GetUndistoPixel(ns_veta::Vec2d(kp.x, kp.y));
            index.at(j) = j;
        }
#pragma omp critical
        { featMap.insert({_frames.at(i)->GetId(), {index, kps, descriptor, kpsUndisto}}); }
    }
    spdlog::info("feature extraction finished.");

    // ----------------
    // feature matching
    // ----------------
    spdlog::info("start matching exhaustive features, this would cost some time...");
    std::set<IndexPair> hasDone;
    for (int j = 0; j < static_cast<int>(_frames.size()); ++j) {
        const auto &refFrame = _frames.at(j);
        auto covFrames = FindCovisibility(refFrame, hasDone);
        // spdlog::info("covisibility view count of refFrame '{}': '{}'", refFrame->GetId(),
        // covFrames.size());

        const ns_veta::IndexT &refId = refFrame->GetId();
        const auto &refFeat = featMap.at(refId);

        for (const auto &[schFrame, intersection] : covFrames) {
            const ns_veta::IndexT &schId = schFrame->GetId();

            spdlog::info("performing feature matching between frames '{}'-'{}'...", refId, schId);

            // record match info
            hasDone.insert({schId, refId});
            hasDone.insert({refId, schId});

            // extract in-broder key points
            const auto &[polySchInRef, polyRefInSch] = intersection;
            const auto &schFeat = featMap.at(schId);
            auto refFeatInBorder = FindInBorderOnes(refFeat, polySchInRef);
            auto schFeatInBorder = FindInBorderOnes(schFeat, polyRefInSch);

            // matching
            auto [refMatched, schMatched] = MatchFeatures(refFeatInBorder, schFeatInBorder);

            // outlier rejection
            const auto &inlierIdx = RejectOutliers(refMatched, schMatched, _intri);

            // ransac failed
            if (inlierIdx.empty()) {
                std::cout << "\033[1A";
                // spdlog::warn("feature matching failed for co-visible frame '{}' and '{}'", refId,
                // schId);
                continue;
            }

            SfMFeaturePairVec featPairVec(inlierIdx.size());
            for (int i = 0; i < static_cast<int>(inlierIdx.size()); ++i) {
                const int idx = inlierIdx.at(i);
                auto &featPair = featPairVec.at(i);
                featPair.first = refMatched.at(idx);
                featPair.second = schMatched.at(idx);
            }

            _matchRes.insert(
                {IndexPair(refId, schId), {refId, schId, featPairVec, polySchInRef, polyRefInSch}});
            _viewFeatLM.insert({refId, {}});
            _viewFeatLM.insert({schId, {}});

            std::cout << "\033[1A";
#define VISUALIZATION 0
#if VISUALIZATION
            const cv::Scalar red = cv::Scalar(0, 0, 255);
            const cv::Scalar green = cv::Scalar(0, 255, 0);
            auto covImg1 =
                DrawInBorderFeatMatch(refFrame->GetImage(), polySchInRef, refFeat, refFeatInBorder);
            auto covImg2 =
                DrawInBorderFeatMatch(schFrame->GetImage(), polyRefInSch, schFeat, schFeatInBorder);

            cv::Mat covImg;
            cv::hconcat(covImg1, covImg2, covImg);

            assert(refMatched.size() == schMatched.size());
            spdlog::info(" match rate: {}/{}", refMatched.size(),
                         std::get<0>(refFeatInBorder).size());
            spdlog::info("inlier rate: {}/{}", inlierIdx.size(), refMatched.size());

            auto bias = cv::Point2f((float)refFrame->GetImage().cols, 0.0f);
            for (int i = 0; i < static_cast<int>(refMatched.size()); ++i) {
                const auto &pt1 = refMatched.at(i).kp, pt2 = schMatched.at(i).kp + bias;
                cv::line(covImg, pt1, pt2, red, 1);
            }
            for (const auto &i : inlierIdx) {
                const auto &pt1 = refMatched.at(i).kp, pt2 = schMatched.at(i).kp + bias;
                cv::drawMarker(covImg, pt1, green, cv::MarkerTypes::MARKER_SQUARE, 10, 1);
                cv::drawMarker(covImg, pt1, green, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
                cv::drawMarker(covImg, pt2, green, cv::MarkerTypes::MARKER_SQUARE, 10, 1);
                cv::drawMarker(covImg, pt2, green, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
                cv::line(covImg, pt1, pt2, green, 1);
            }

            cv::imshow("Covisibility", covImg);
            // static int count = 0;
            // cv::imwrite(
            //         Configor::DataStream::DebugPath + "/covisibility" + std::to_string(++count) +
            //         ".png", covImg
            // );
            cv::waitKey(0);
#endif
#undef VISUALIZATION
        }

        if (j % 10 == 0) {
            _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION)
                .AddEntityLocal(_viewCubes, Viewer::VIEW_ASSOCIATION);
            DrawMatchesInViewer(ns_viewer::Colour::Red());
        }
    }
    std::cout << std::endl;
    spdlog::info("feature matching finished.");
    hasDone.clear();
    featMap.clear();
    _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION)
        .AddEntityLocal(_viewCubes, Viewer::VIEW_ASSOCIATION);
    DrawMatchesInViewer(ns_viewer::Colour::Green());

    spdlog::info("checking the covisibility graph connection...");
    if (!IsGraphConnect(ExtractKeysAsVec(_matchRes))) {
        spdlog::warn(
            "the covisibility graph is not connected!!! "
            "some frames do not have covisibility with other ones!!!");
    } else {
        spdlog::info("the covisibility graph is connected!");
    }
    return true;
}

bool VisionOnlySfM::StructureFromMotion() {
    // ---------------------
    // structure from motion
    // ---------------------
    spdlog::info(
        "initialize structure using the frame pair with enough covisibility and parallax...");
    auto initViewIdxPair = InitStructure();
    _viewer->AddVeta(_veta, Viewer::VIEW_MAP);
    cv::Mat img = DrawMatchResult(initViewIdxPair);
    cv::imshow("img", img);
    cv::waitKey(0);

    spdlog::info("performing incremental structure from motion...");
    IncrementalSfM(initViewIdxPair);
    return true;
}

std::optional<ns_veta::Posed> VisionOnlySfM::ComputeCamRotations(const CameraFrame::Ptr &frame) {
    const double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(_topic);
    double tByBr = frame->GetTimestamp() + TO_CmToBr;
    if (!_so3Spline.TimeStampInRange(tByBr)) {
        return {};
    } else {
        auto SO3_BrToBr0 = _so3Spline.Evaluate(tByBr);
        const auto &SO3_CmToBr = _parMagr->EXTRI.SO3_CmToBr.at(_topic);
        // the translation is as zero (identity)
        return ns_veta::Posed(SO3_BrToBr0 * SO3_CmToBr);
    }
}

// -------------------------------
// VisionOnlySfM: feature matching
// -------------------------------

std::map<CameraFrame::Ptr, std::pair<polygon_2d, polygon_2d>> VisionOnlySfM::FindCovisibility(
    const CameraFrame::Ptr &refFrame, const std::set<IndexPair> &ignore, double covThd) {
    const auto &view = _veta->views.at(refFrame->GetId());
    const double area = (double)view->imgHeight * (double)view->imgWidth;
    const auto &refSo3 = _veta->poses.at(view->poseId).Rotation();
    const auto &refSo3Inv = refSo3.inverse();

    std::map<CameraFrame::Ptr, std::pair<polygon_2d, polygon_2d>> covFrames;
    for (const auto &schFrame : _frames) {
        // same frame
        if (schFrame == refFrame) {
            continue;
        }

        // has been performed
        if (ignore.find({refFrame->GetId(), schFrame->GetId()}) != ignore.cend()) {
            continue;
        }
        if (ignore.find({schFrame->GetId(), refFrame->GetId()}) != ignore.cend()) {
            continue;
        }

        const auto &schSo3 = _veta->poses.at(_veta->views.at(schFrame->GetId())->poseId).Rotation();
        auto intersection =
            IntersectionArea(refFrame->GetImage(), schFrame->GetImage(), refSo3Inv * schSo3);
        if (!intersection) {
            continue;
        }

        const auto &[sectSchInRef, sectRefInSch] = *intersection;

        double covRate1 = bg::area(sectSchInRef) / area;
        double covRate2 = bg::area(sectRefInSch) / area;
        // spdlog::info("covisibility area between '{}' and '{}': '{:.3f}' | '{:.3f}'",
        //              refFrame->GetId(), schFrame->GetId(), covRate1, covRate2);

        if (covRate1 < covThd || covRate2 < covThd) {
            continue;
        }

        covFrames.insert({schFrame, *intersection});
    }
    return covFrames;
}

std::optional<std::pair<polygon_2d, polygon_2d>> VisionOnlySfM::IntersectionArea(
    const cv::Mat &i1, const cv::Mat &i2, const Sophus::SO3d &SO3_2To1) {
    auto poly2In1 = ProjPolygon(i2, SO3_2To1, _intri);
    auto poly1 = ProjPolygon(i1, Sophus::SO3d(), _intri);
    if (!poly2In1 || !poly1) {
        return {};
    }

    auto poly1In2 = ProjPolygon(i1, SO3_2To1.inverse(), _intri);
    auto poly2 = ProjPolygon(i2, Sophus::SO3d(), _intri);
    if (!poly1In2 || !poly2) {
        return {};
    }

    std::vector<polygon_2d> output1, output2;
    boost::geometry::intersection(*poly1, *poly2In1, output1);
    if (output1.empty()) {
        return {};
    }

    boost::geometry::intersection(*poly2, *poly1In2, output2);
    if (output2.empty()) {
        return {};
    }

    // only one if intersection exists
    polygon_2d sect2In1 = output1[0];
    polygon_2d sect1In2 = output2[0];

    return std::pair<polygon_2d, polygon_2d>{sect2In1, sect1In2};
}

std::optional<polygon_2d> VisionOnlySfM::ProjPolygon(const cv::Mat &img,
                                                     const Sophus::SO3d &so3,
                                                     const ns_veta::PinholeIntrinsic::Ptr &intri) {
    double col = img.cols - 1, row = img.rows - 1;

    auto Corner2To1 = [&intri, &so3](double x2, double y2) -> std::optional<point_2d> {
        ns_veta::Vec2d pCam = intri->ImgToCam(ns_veta::Vec2d(x2, y2));
        ns_veta::Vec3d bearing2 = ns_veta::Vec3d(pCam(0), pCam(1), 1.0).normalized();
        ns_veta::Vec3d bearing1 = so3 * bearing2;
        // for negative and small z-value
        if (bearing1(2) < 0.01) {
            return {};
        }
        ns_veta::Vec2d pImg =
            intri->CamToImg(ns_veta::Vec2d(bearing1(0) / bearing1(2), bearing1(1) / bearing1(2)));
        return point_2d(pImg(0), pImg(1));
    };

    polygon_2d polyProj;

    // polygon should be clockWise
    auto p1 = Corner2To1(0, 0), p2 = Corner2To1(col, 0), p3 = Corner2To1(col, row),
         p4 = Corner2To1(0, row);
    if (!p1 || !p2 || !p3 || !p4) {
        return {};
    }
    bg::append(polyProj, *p1);
    bg::append(polyProj, *p2);
    bg::append(polyProj, *p3);
    bg::append(polyProj, *p4);

    // polygon should be clockWise and should be closed
    boost::geometry::correct(polyProj);

    return polyProj;
}

VisionOnlySfM::FeaturePack VisionOnlySfM::FindInBorderOnes(const FeaturePack &input,
                                                           const polygon_2d &poly) {
    const auto &iIdx = std::get<0>(input);
    const auto &iKps = std::get<1>(input);
    const auto &iDesc = std::get<2>(input);
    const auto &iKpsUndisto = std::get<3>(input);

    assert(iIdx.size() == iKps.size() && iIdx.size() == iDesc.rows &&
           iIdx.size() == iKpsUndisto.size());

    const int size = static_cast<int>(iKps.size());

    FeaturePack output;
    auto &oIdx = std::get<0>(output);
    auto &oKps = std::get<1>(output);
    auto &oDesc = std::get<2>(output);
    auto &oKpsUndisto = std::get<3>(output);

    oIdx.reserve(size), oKps.reserve(size), oKpsUndisto.reserve(size);

    std::vector<cv::Mat> validDesc;
    validDesc.reserve(size);

    for (int i = 0; i < size; ++i) {
        const auto &kp = iKps.at(i).pt;
        // if (bg::within(point_2d(kp.x, kp.y), poly) || bg::distance(point_2d(kp.x, kp.y), poly) <
        // 50) {
        if (bg::distance(point_2d(kp.x, kp.y), poly) < 50) {
            oIdx.push_back(iIdx.at(i));
            oKps.push_back(iKps.at(i));
            validDesc.push_back(iDesc.row(i));
            oKpsUndisto.push_back(iKpsUndisto.at(i));
        }
    }

    oIdx.shrink_to_fit(), oKps.shrink_to_fit(), oKpsUndisto.shrink_to_fit();

    cv::vconcat(validDesc, oDesc);

    return output;
}

std::pair<SfMFeatureVec, SfMFeatureVec> VisionOnlySfM::MatchFeatures(
    const VisionOnlySfM::FeaturePack &feat1, const VisionOnlySfM::FeaturePack &feat2) {
    const auto &idx1 = std::get<0>(feat1), idx2 = std::get<0>(feat2);
    const auto &kps1 = std::get<1>(feat1), kps2 = std::get<1>(feat2);
    const auto &desc1 = std::get<2>(feat1), desc2 = std::get<2>(feat2);
    const auto &kpsUndist1 = std::get<3>(feat1), kpsUndist2 = std::get<3>(feat2);

    cv::BFMatcher matcher(cv::NORM_HAMMING);
    std::vector<std::vector<cv::DMatch>> nnMatches;
    matcher.knnMatch(desc1, desc2, nnMatches, 2);

    SfMFeatureVec matched1, matched2;
    constexpr double nn_match_ratio = 0.8f;
    std::set<int> hasMatched;
    for (auto &match : nnMatches) {
        cv::DMatch first = match[0];
        auto qi = first.queryIdx, ti = first.trainIdx;

        if (hasMatched.find(first.trainIdx) != hasMatched.cend()) {
            continue;
        }

        float dist1 = match[0].distance;
        float dist2 = match[1].distance;
        if (dist1 < nn_match_ratio * dist2) {
            matched1.emplace_back(idx1.at(qi), kps1.at(qi).pt, kpsUndist1.at(qi));
            matched2.emplace_back(idx2.at(ti), kps2.at(ti).pt, kpsUndist2.at(ti));
            hasMatched.insert(ti);
        }
    }
    return {matched1, matched2};
}

std::vector<int> VisionOnlySfM::RejectOutliers(const SfMFeatureVec &feat1,
                                               const SfMFeatureVec &feat2,
                                               const ns_veta::PinholeIntrinsic::Ptr &intri) {
    // matched pair size are not enough
    if (feat1.size() < 15) {
        return {};
    }

    opengv::bearingVectors_t bearingVectors1 = ComputeBeringVec(feat1, intri);
    opengv::bearingVectors_t bearingVectors2 = ComputeBeringVec(feat2, intri);

    std::vector<int> inliers;

    std::vector<int> idxRecorder;
    opengv::rotation_t rot;
    {
        // -------------
        // rotation only
        // -------------
        opengv::relative_pose::CentralRelativeAdapter adapter(bearingVectors1, bearingVectors2);

        // Create a RotationOnlySacProblem and Ransac
        opengv::sac::Ransac<opengv::sac_problems::relative_pose::RotationOnlySacProblem> ransac;
        std::shared_ptr<opengv::sac_problems::relative_pose::RotationOnlySacProblem> probPtr(
            new opengv::sac_problems::relative_pose::RotationOnlySacProblem(adapter));
        ransac.sac_model_ = probPtr;
        ransac.threshold_ = intri->ImagePlaneToCameraPlaneError(1.0);
        ransac.max_iterations_ = 50;

        auto res = ransac.computeModel(0);

        if (!res) {
            return {};
        }

        idxRecorder.resize(ransac.inliers_.size());
        rot = ransac.model_coefficients_;

        for (int i = 0; i < static_cast<int>(ransac.inliers_.size()); ++i) {
            idxRecorder.at(i) = ransac.inliers_.at(i);
        }

        auto FindInliers = [](const opengv::bearingVectors_t &iBeringVec,
                              const std::vector<int> &inliers) {
            opengv::bearingVectors_t oBeringVec(inliers.size());
            for (int i = 0; i < static_cast<int>(inliers.size()); ++i) {
                oBeringVec.at(i) = iBeringVec.at(inliers.at(i));
            }
            return oBeringVec;
        };
        bearingVectors1 = FindInliers(bearingVectors1, ransac.inliers_);
        bearingVectors2 = FindInliers(bearingVectors2, ransac.inliers_);
    }

    {
        // ----------------
        // translation only
        // ----------------
        opengv::relative_pose::CentralRelativeAdapter adapter(bearingVectors1, bearingVectors2,
                                                              rot);

        opengv::sac::Ransac<opengv::sac_problems::relative_pose::TranslationOnlySacProblem> ransac;
        std::shared_ptr<opengv::sac_problems::relative_pose::TranslationOnlySacProblem> probPtr(
            new opengv::sac_problems::relative_pose::TranslationOnlySacProblem(adapter));
        ransac.sac_model_ = probPtr;
        ransac.threshold_ = intri->ImagePlaneToCameraPlaneError(1.0);
        ransac.max_iterations_ = 20;

        auto res = ransac.computeModel(0);
        if (!res) {
            return {};
        }

        inliers.resize(ransac.inliers_.size());
        for (int i = 0; i < static_cast<int>(ransac.inliers_.size()); ++i) {
            inliers.at(i) = idxRecorder.at(ransac.inliers_.at(i));
        }
    }

    if (inliers.size() < 15) {
        return {};
    }

    return inliers;
}

// ------------------------------------
// VisionOnlySfM: structure from motion
// ------------------------------------

std::map<IndexPair, Eigen::Vector3d> VisionOnlySfM::Triangulate(const SfMFeaturePairInfo &featPair,
                                                                const Eigen::Matrix3d &R12,
                                                                const Eigen::Vector3d &t12,
                                                                double reProjThd,
                                                                double parallaxDegThd) {
    Eigen::Matrix3d R21 = R12.transpose();
    Eigen::Vector3d t21 = -R12.transpose() * t12;

    const Eigen::Vector3d &oCam1InCam1 = Eigen::Vector3d::Zero();
    const Eigen::Vector3d &oCam2InCam1 = t12;
    std::map<IndexPair, Eigen::Vector3d> lmsCur;

    const double cosParallaxThd = std::cos(parallaxDegThd / 180.0 * M_PI);

    const auto &featPairVec = featPair.featPairVec;
    const auto &[viewId1, viewId2] = featPair.viewId;

    std::vector<int> idxRecorder;
    idxRecorder.reserve(featPairVec.size());
    for (int i = 0; i < static_cast<int>(featPairVec.size()); ++i) {
        const auto &[feat1, feat2] = featPairVec.at(i);
        // has been triangulated
        if (FindLMByViewFeatId(viewId1, feat1.id) != ns_veta::UndefinedIndexT) {
            continue;
        }
        if (FindLMByViewFeatId(viewId2, feat2.id) != ns_veta::UndefinedIndexT) {
            continue;
        }

        idxRecorder.push_back(i);
    }

    if (idxRecorder.empty()) {
        return {};
    }

#if USE_OPENCV_TRIANGULATE
    // identity
    cv::Mat T1 = (cv::Mat_<float>(3, 4) << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0);
    // from fir to sed
    cv::Mat T2 = (cv::Mat_<float>(3, 4) << R21(0, 0), R21(0, 1), R21(0, 2), t21(0, 0), R21(1, 0),
                  R21(1, 1), R21(1, 2), t21(1, 0), R21(2, 0), R21(2, 1), R21(2, 2), t21(2, 0));
    std::vector<cv::Point2f> pts1(idxRecorder.size()), pts2(idxRecorder.size());
    for (int i = 0; i < static_cast<int>(idxRecorder.size()); ++i) {
        const auto &[feat1, feat2] = featPairVec.at(idxRecorder.at(i));

        ns_veta::Vec2f p1 = _intri->ImgToCam(feat1.kpUndist).cast<float>();
        pts1.at(i).x = p1(0), pts1.at(i).y = p1(1);

        ns_veta::Vec2f p2 = _intri->ImgToCam(feat2.kpUndist).cast<float>();
        pts2.at(i).x = p2(0), pts2.at(i).y = p2(1);
    }
    // points are parameterized in first view
    cv::Mat pts4d;
    cv::triangulatePoints(T1, T2, pts1, pts2, pts4d);
#else
    auto bearingVec1 = ComputeBeringVec<0>(featPairVec, *_intri, idxRecorder);
    auto bearingVec2 = ComputeBeringVec<1>(featPairVec, *_intri, idxRecorder);

    opengv::relative_pose::CentralRelativeAdapter adapter(bearingVec1, bearingVec2, t12, R12);
#endif
    const double reProjThd2 = reProjThd * reProjThd;
    for (int i = 0; i < static_cast<int>(idxRecorder.size()); ++i) {
        const auto &[feat1, feat2] = featPairVec.at(idxRecorder.at(i));

#if USE_OPENCV_TRIANGULATE
        cv::Mat x = pts4d.col(i);
        x /= x.at<float>(3, 0);
        Eigen::Vector3d pInCam1(x.at<float>(0, 0), x.at<float>(1, 0), x.at<float>(2, 0));
#else
        Eigen::Vector3d pInCam1 = opengv::triangulation::triangulate(adapter, i);
#endif

        Eigen::Vector3d pInCam2 = R21 * pInCam1 + t21;

        // check validity
        if (pInCam1(2) < 1E-3 || pInCam2(2) < 1E-3) {
            continue;
        }

        // check parallax
        Eigen::Vector3d normal1 = pInCam1 - oCam1InCam1;
        Eigen::Vector3d normal2 = pInCam1 - oCam2InCam1;
        double cosParallax = normal1.dot(normal2) / (normal1.norm() * normal2.norm());

        // std::cout << "cos parallax: " << cosParallax << std::endl;
        if (cosParallax > cosParallaxThd) {
            continue;
        }

        // check reprojection error in first image
        double invZ1 = 1.0 / pInCam1(2);
        ns_veta::Vec2d ip1 = _intri->CamToImg({pInCam1(0) * invZ1, pInCam1(1) * invZ1});
        double error1 = (feat1.kpUndist - ip1).squaredNorm();
        // std::cout << "error 1: " << error1 << std::endl;
        if (error1 > reProjThd2) {
            continue;
        }

        // check reprojection error in second image
        double invZ2 = 1.0 / pInCam2(2);
        ns_veta::Vec2d ip2 = _intri->CamToImg({pInCam2(0) * invZ2, pInCam2(1) * invZ2});
        double error2 = (feat2.kpUndist - ip2).squaredNorm();
        // std::cout << "error 2: " << error2 << std::endl;
        if (error2 > reProjThd2) {
            continue;
        }

        lmsCur.insert({{feat1.id, feat2.id}, pInCam1});
    }

    spdlog::info("valid landmarks in 'VisionOnlySfM::Triangulate': {}/{}", lmsCur.size(),
                 idxRecorder.size());

    return lmsCur;
}

bool VisionOnlySfM::IsGraphConnect(const std::vector<IndexPair> &graph) {
    std::set<ns_veta::IndexT> s;
    for (const auto &edge : graph) {
        s.insert(edge.first), s.insert(edge.second);
    }
    std::size_t numNodes = s.size();

    std::set<ns_veta::IndexT> visited;
    std::queue<ns_veta::IndexT> bfsQueue;

    // Start BFS from an arbitrary node
    bfsQueue.push(graph[0].first);
    visited.insert(graph[0].first);

    while (!bfsQueue.empty()) {
        ns_veta::IndexT curNode = bfsQueue.front();
        bfsQueue.pop();

        for (const auto &edge : graph) {
            if (edge.first == curNode && visited.find(edge.second) == visited.cend()) {
                bfsQueue.push(edge.second);
                visited.insert(edge.second);
            }
            if (edge.second == curNode && visited.find(edge.first) == visited.cend()) {
                bfsQueue.push(edge.first);
                visited.insert(edge.first);
            }
        }
    }

    return visited.size() == numNodes;
}

IndexPair VisionOnlySfM::InitStructure() {
    auto [maxMatchIdxPair, maxMatchPairInfo] =
        *std::max_element(_matchRes.begin(), _matchRes.end(), [](const auto &p1, const auto &p2) {
            return p1.second.featPairVec.size() < p2.second.featPairVec.size();
        });

    auto maxMatchCount = maxMatchPairInfo.featPairVec.size();
    spdlog::info("max match count pair: '{}'-'{}', max valid count: {}", maxMatchIdxPair.first,
                 maxMatchIdxPair.second, maxMatchCount);
    std::vector<IndexPair> matchPairPool;
    for (const auto &[viewIdxPair, infoPair] : _matchRes) {
        // matches are insufficient
        if ((double)infoPair.featPairVec.size() < 0.5 * (double)maxMatchCount) {
            continue;
        }
        // parallax is probably small
        // if (std::abs((int) viewIdxPair.first - (int) viewIdxPair.second) == 1) { continue; }

        std::set<ns_veta::IndexT> idxSet = {maxMatchIdxPair.first, maxMatchIdxPair.second,
                                            viewIdxPair.first, viewIdxPair.second};
        if (idxSet.size() < 4) {
            matchPairPool.push_back(viewIdxPair);
        }
    }

    // todo: refine this best search
    std::map<IndexPair, Eigen::Vector3d> lmBest;
    IndexPair viewPairBest;
    Sophus::SE3d se3Best;
    for (const auto &viewIdxPair : matchPairPool) {
        const auto &pairInfo = _matchRes.at(viewIdxPair);
        auto SE3_2To1 = pairInfo.RecoverRelativePose(_intri);
        if (!SE3_2To1) {
            continue;
        }

        auto lmsCur = Triangulate(_matchRes.at(viewIdxPair), SE3_2To1->rotationMatrix(),
                                  SE3_2To1->translation());
        if (lmsCur.size() > lmBest.size()) {
            lmBest = lmsCur, viewPairBest = viewIdxPair, se3Best = *SE3_2To1;
        }
    }

    spdlog::info("best image pair: {}-{}, match count: {}", viewPairBest.first, viewPairBest.second,
                 lmBest.size());

    // the first view is the world frame
    InsertTriangulateLM(_matchRes.at(viewPairBest), lmBest, ns_veta::Posed());

    // clear rotation-only rough poses and add new initial pose
    _veta->poses.clear();
    _veta->poses.insert({viewPairBest.first, ns_veta::Posed()});
    _veta->poses.insert(
        {viewPairBest.second, ns_veta::Posed(se3Best.so3(), se3Best.translation())});

    return viewPairBest;
}

void VisionOnlySfM::IncrementalSfM(const IndexPair &sViewIdxPair) {
    std::set<ns_veta::IndexT> visited;
    std::queue<ns_veta::IndexT> bfsQueue;

    // Start BFS from an arbitrary node
    bfsQueue.push(sViewIdxPair.first);
    visited.insert(sViewIdxPair.first);
    std::set<IndexPair> viewPairsHaveDone;
    // such pair has been in the 'InitStructure'
    viewPairsHaveDone.insert(sViewIdxPair);
    int closureCount = 0;

    while (!bfsQueue.empty()) {
        ns_veta::IndexT curNode = bfsQueue.front();
        bfsQueue.pop();

        for (auto &[edge, pairInfo] : _matchRes) {
            // no connection
            if (edge.first != curNode && edge.second != curNode) {
                continue;
            }
            // such pair has been performed
            if (viewPairsHaveDone.find(edge) != viewPairsHaveDone.cend()) {
                continue;
            }

            // connect but not been performed
            if (edge.first == curNode && visited.find(edge.second) == visited.cend()) {
                // perform absolute pose solving for 'edge.second'
                if (!SolveAbsolutePose<0, 1>(pairInfo)) {
                    continue;
                }

                auto T_FirToW = _veta->poses.at(edge.first);
                auto T_SedToW = _veta->poses.at(edge.second);
                // from 1 to 0 (from second to first)
                auto T_SedToFir = T_FirToW.Inverse() * T_SedToW;

                // landmarks are parameterized in first view
                auto lms =
                    Triangulate(pairInfo, T_SedToFir.Rotation().matrix(), T_SedToFir.Translation());
                InsertTriangulateLM(pairInfo, lms, T_FirToW);

                bfsQueue.push(edge.second);
                visited.insert(edge.second);

                // cv::Mat img = DrawMatchResult(edge);
                // cv::imshow("img", img);
            } else if (edge.second == curNode && visited.find(edge.first) == visited.cend()) {
                // perform absolute pose solving for 'edge.first'
                if (!SolveAbsolutePose<1, 0>(pairInfo)) {
                    continue;
                }

                auto T_FirToW = _veta->poses.at(edge.first);
                auto T_SedToW = _veta->poses.at(edge.second);
                // from 1 to 0 (from second to first)
                auto T_SedToFir = T_FirToW.Inverse() * T_SedToW;

                // landmarks are parameterized in first view
                auto lms =
                    Triangulate(pairInfo, T_SedToFir.Rotation().matrix(), T_SedToFir.Translation());
                InsertTriangulateLM(pairInfo, lms, T_FirToW);

                bfsQueue.push(edge.first);
                visited.insert(edge.first);

                // cv::Mat img = DrawMatchResult(edge);
                // cv::imshow("img", img);
            } else if (visited.find(edge.first) != visited.cend() &&
                       visited.find(edge.second) != visited.cend()) {
                ++closureCount;
                // closure optimization
                if (closureCount % 100 == 0) {
                    BatchOptimization();
                }

                auto T_FirToW = _veta->poses.at(edge.first);
                auto T_SedToW = _veta->poses.at(edge.second);
                // from 1 to 0 (from second to first)
                auto T_SedToFir = T_FirToW.Inverse() * T_SedToW;

                // landmarks are parameterized in first view
                auto lms =
                    Triangulate(pairInfo, T_SedToFir.Rotation().matrix(), T_SedToFir.Translation());
                InsertTriangulateLM(pairInfo, lms, T_FirToW);
            }

            viewPairsHaveDone.insert(edge);
            if (viewPairsHaveDone.size() % 10 == 0) {
                _viewer->ClearViewer(Viewer::VIEW_MAP).AddVeta(_veta, Viewer::VIEW_MAP);
                _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION)
                    .AddEntityLocal(_viewCubes, Viewer::VIEW_ASSOCIATION);
                DrawMatchesInViewer(ns_viewer::Colour::Green(), viewPairsHaveDone,
                                    ns_viewer::Colour::Black());
                // cv::waitKey(0);
            }
        }
    }
    _viewer->ClearViewer(Viewer::VIEW_MAP).AddVeta(_veta, Viewer::VIEW_MAP);
    _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION)
        .AddEntityLocal(_viewCubes, Viewer::VIEW_ASSOCIATION);
    DrawMatchesInViewer(ns_viewer::Colour::Green(), viewPairsHaveDone, ns_viewer::Colour::Black());
}

void VisionOnlySfM::InsertTriangulateLM(const SfMFeaturePairInfo &featPair,
                                        const std::map<IndexPair, Eigen::Vector3d> &lms,
                                        const ns_veta::Posed &curLMToW) {
    const auto &[viewId1, viewId2] = featPair.viewId;
    for (auto &[feat1, feat2] : featPair.featPairVec) {
        auto iter = lms.find({feat1.id, feat2.id});
        if (iter == lms.cend()) {
            continue;
        }

        // add landmarks
        ns_veta::Landmark lm(curLMToW(iter->second), {});
        lm.obs.insert({viewId1, {feat1.kpUndist, feat1.id}});
        lm.obs.insert({viewId2, {feat2.kpUndist, feat2.id}});

        auto lmId = GenNewLMId();
        if (_veta->structure.insert({lmId, lm}).second) {
            if (!InsertViewFeatLM(viewId1, feat1.id, lmId)) {
                spdlog::warn(
                    "[VisionOnlySfM::InitStructure] view feat has connected to landmark "
                    "but not inserted to 'ViewFeatLM'");
            }
            if (!InsertViewFeatLM(viewId2, feat2.id, lmId)) {
                spdlog::warn(
                    "[VisionOnlySfM::InitStructure] view feat has connected to landmark "
                    "but not inserted to 'ViewFeatLM'");
            }
        }
    }
}

// -----------------
// drawing functions
// -----------------

polygon_2d VisionOnlySfM::BufferPolygon(const polygon_2d &polygon, double bufferDistance) {
    const int pointsPerCircle = 36;
    bg::strategy::buffer::distance_symmetric<double> distanceStrategy(bufferDistance);
    bg::strategy::buffer::join_round joinStrategy(pointsPerCircle);
    bg::strategy::buffer::end_round endStrategy(pointsPerCircle);
    bg::strategy::buffer::point_circle circleStrategy(pointsPerCircle);
    bg::strategy::buffer::side_straight sideStrategy;
    bg::model::multi_polygon<polygon_2d> result;
    bg::buffer(polygon, result, distanceStrategy, sideStrategy, joinStrategy, endStrategy,
               circleStrategy);
    return result.at(0);
}

void VisionOnlySfM::DrawProjPolygon(cv::Mat &img,
                                    const polygon_2d &poly,
                                    const cv::Scalar &fill,
                                    const cv::Scalar &border) {
    std::vector<cv::Point> pts;
    for (const auto &p : poly.outer()) {
        pts.emplace_back((int)p.x(), (int)p.y());
    }

    cv::Mat canvas = cv::Mat::zeros(cv::Size(img.cols, img.rows), CV_8UC3);
    const cv::Point *ppts[1] = {&pts[0]};
    const int n = static_cast<int>(pts.size()), *npt = &n;
    cv::fillPoly(canvas, ppts, npt, 1, fill);

    cv::addWeighted(img, 1.0, canvas, 0.3, 0.0, img);
    cv::polylines(img, pts, true, border, 2, cv::LineTypes::LINE_AA);
}

cv::Mat VisionOnlySfM::DrawMatchResult(const IndexPair &viewIdPair) {
    const auto &featPairInfo = _matchRes.at(viewIdPair);

    const cv::Scalar white = cv::Scalar(255, 255, 255);
    const cv::Scalar black = cv::Scalar(0, 0, 0);
    const cv::Scalar blue = cv::Scalar(255, 0, 0);
    const cv::Scalar green = cv::Scalar(0, 255, 0);
    const cv::Scalar red = cv::Scalar(0, 0, 255);

    cv::Mat img1, img2;
    cv::cvtColor(_frameBackup.at(viewIdPair.first)->GetImage(), img1, cv::COLOR_GRAY2BGR);
    cv::cvtColor(_frameBackup.at(viewIdPair.second)->GetImage(), img2, cv::COLOR_GRAY2BGR);

    auto bufPloy1 = BufferPolygon(featPairInfo.poly.first, 50.0);
    auto bufPloy2 = BufferPolygon(featPairInfo.poly.second, 50.0);

    DrawProjPolygon(img1, bufPloy1, black, blue);
    DrawProjPolygon(img1, featPairInfo.poly.first, white, blue);

    DrawProjPolygon(img2, bufPloy2, black, blue);
    DrawProjPolygon(img2, featPairInfo.poly.second, white, blue);

    cv::Mat covImg;
    cv::hconcat(img1, img2, covImg);

    auto bias = cv::Point2f((float)img2.cols, 0.0f);
    for (const auto &[feat1, feat2] : featPairInfo.featPairVec) {
        // such match is invalid
        cv::Scalar color;
        auto lm1 = FindLMByViewFeatId(viewIdPair.first, feat1.id);
        auto lm2 = FindLMByViewFeatId(viewIdPair.second, feat2.id);
        if (lm1 == ns_veta::UndefinedIndexT || lm2 == ns_veta::UndefinedIndexT || lm1 != lm2) {
            color = red;
        } else {
            color = green;
        }

        const auto &pt1 = feat1.kp, pt2 = feat2.kp + bias;
        cv::drawMarker(covImg, pt1, color, cv::MarkerTypes::MARKER_SQUARE, 10, 1);
        cv::drawMarker(covImg, pt1, color, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
        cv::drawMarker(covImg, pt2, color, cv::MarkerTypes::MARKER_SQUARE, 10, 1);
        cv::drawMarker(covImg, pt2, color, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
        cv::line(covImg, pt1, pt2, color, 1);
    }

    return covImg;
}

cv::Mat VisionOnlySfM::DrawInBorderFeatMatch(const cv::Mat &img,
                                             const polygon_2d &poly,
                                             const FeaturePack &allFeats,
                                             const FeaturePack &ibFeats) {
    constexpr double bufferDistance = 50.0;

    const cv::Scalar white = cv::Scalar(255, 255, 255);
    const cv::Scalar black = cv::Scalar(0, 0, 0);
    const cv::Scalar blue = cv::Scalar(255, 0, 0);

    cv::Mat covImg;
    cv::cvtColor(img, covImg, cv::COLOR_GRAY2BGR);

    DrawProjPolygon(covImg, BufferPolygon(poly, bufferDistance), black, blue);
    DrawProjPolygon(covImg, poly, white, blue);

    for (const auto &kp : std::get<1>(allFeats)) {
        const auto &pt = kp.pt;
        cv::drawMarker(covImg, pt, cv::Scalar(255, 255, 255), cv::MarkerTypes::MARKER_SQUARE, 10,
                       1);
        cv::drawMarker(covImg, pt, cv::Scalar(255, 255, 255), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
    }

    for (const auto &kp : std::get<1>(ibFeats)) {
        const auto &pt = kp.pt;
        cv::drawMarker(covImg, pt, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 10, 1);
        cv::drawMarker(covImg, pt, cv::Scalar(0, 0, 255), cv::MarkerTypes::MARKER_SQUARE, 2, 2);
    }

    return covImg;
}

void VisionOnlySfM::DrawMatchesInViewer(const ns_viewer::Colour &color,
                                        const std::set<IndexPair> &special,
                                        const ns_viewer::Colour &specialColor) const {
    std::vector<ns_viewer::Entity::Ptr> lines;
    lines.reserve(_matchRes.size());
    for (const auto &[viewIdPair, info] : _matchRes) {
        // as visualization is heavy, we do not draw the small connection
        if (info.featPairVec.size() < 600) {
            continue;
        }

        ns_viewer::Colour c;
        if (special.count(viewIdPair) != 0) {
            c = specialColor;
        } else {
            c = color;
        }

        Eigen::Vector3f p1 = _viewCubePoses.at(viewIdPair.first).translation().cast<float>();
        Eigen::Vector3f p2 = _viewCubePoses.at(viewIdPair.second).translation().cast<float>();

        lines.push_back(ns_viewer::Line::Create(p1, p2, c, 1.0f));
    }
    _viewer->AddEntityLocal(lines, Viewer::VIEW_ASSOCIATION);
}

void VisionOnlySfM::CreateViewCubes() {
    // visualization
    auto poses = GenerateUniformPoseOnSphere(static_cast<int>(_veta->views.size()), 5.0);
    std::sort(poses.begin(), poses.end(), [](const Sophus::SE3d &e1, const Sophus::SE3d &e2) {
        double angle1 = std::atan2(e1.translation().y(), e1.translation().x());
        double angle2 = std::atan2(e2.translation().y(), e2.translation().x());
        return angle1 < angle2;
    });
    std::sort(poses.begin(), poses.end(), [](const Sophus::SE3d &e1, const Sophus::SE3d &e2) {
        return e1.translation().z() < e2.translation().z();
    });
    _viewCubes.resize(poses.size());
    auto viewIter = _veta->views.cbegin();
    for (int i = 0; i < static_cast<int>(poses.size()); ++i, ++viewIter) {
        const auto &pose = poses.at(i);
        _viewCubes.at(i) = ns_viewer::Cube::Create(
            ns_viewer::Posed(pose.so3().matrix(), pose.translation()).cast<float>(), true, 0.2f,
            0.2f, 0.2f);
        _viewCubePoses.insert({viewIter->first, pose});
    }
    poses.clear();
}

void VisionOnlySfM::BatchOptimization() {
    auto estimator = Estimator::Create(nullptr, nullptr);

    for (auto &[lmId, lm] : _veta->structure) {
        for (const auto &[viewId, feat] : lm.obs) {
            estimator->AddVisualProjectionFactor(&_veta->poses.at(viewId), &lm.X, _intri, feat.x,
                                                 1.0);
        }
    }
    auto sum = estimator->Solve(Estimator::DefaultSolverOptions(
        Configor::Preference::AvailableThreads(), true, Configor::Preference::UseCudaInSolving));
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
}

const std::map<IndexPair, SfMFeaturePairInfo> &VisionOnlySfM::GetMatchRes() const {
    return _matchRes;
}

std::set<IndexPair> VisionOnlySfM::FindCovisibility(double covThd) {
    _veta = ns_veta::Veta::Create();

    constexpr ns_veta::IndexT intriIdx = 0;

    // -----------------
    // create intrinsics
    // -----------------
    _intri = std::make_shared<ns_veta::PinholeIntrinsic>(*_parMagr->INTRI.Camera.at(_topic));
    _veta->intrinsics.insert({intriIdx, _intri});

    // ----------------------
    // create views and poses
    // ----------------------
    for (const auto &frame : _frames) {
        const ns_veta::IndexT viewId = frame->GetId();
        const ns_veta::IndexT poseId = viewId;
        _frameBackup.insert({viewId, frame});

        auto pose = ComputeCamRotations(frame);
        if (pose) {
            // rough rotation-only pose
            _veta->poses.insert({poseId, *pose});

            // view
            auto view = ns_veta::View::Create(frame->GetTimestamp(), viewId, intriIdx, poseId,
                                              _intri->imgWidth, _intri->imgHeight);
            _veta->views.insert({viewId, view});
        }
    }

    CreateViewCubes();
    _viewer->AddEntity(_viewCubes, Viewer::VIEW_ASSOCIATION);

    std::set<IndexPair> hasDone, covPairs;
    for (int j = 0; j < static_cast<int>(_frames.size()); ++j) {
        const auto &refFrame = _frames.at(j);
        auto covFrames = FindCovisibility(refFrame, hasDone, covThd);
        // spdlog::info("covisibility view count of refFrame '{}': '{}'", refFrame->GetId(),
        // covFrames.size());
        if (j % 10 == 0) {
            _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION)
                .AddEntityLocal(_viewCubes, Viewer::VIEW_ASSOCIATION);
            DrawMatchesInViewer(ns_viewer::Colour::Red());
        }
        for (const auto &[covFrame, _] : covFrames) {
            covPairs.insert({refFrame->GetId(), covFrame->GetId()});
        }
    }
    return covPairs;
}

// ---------------
// SfMFeaturePairInfo
// ---------------

std::optional<Sophus::SE3d> SfMFeaturePairInfo::RecoverRelativePose(
    const ns_veta::PinholeIntrinsic::Ptr &intri) const {
    if (featPairVec.size() < 10) {
        return {};
    }

#if USE_OPENCV_REL_POSE
    std::vector<cv::Point2f> ll, rr;
    for (const auto &[feat1, feat2] : featPairVec) {
        ll.emplace_back(feat1.kpUndist(0), feat1.kpUndist(1));
        rr.emplace_back(feat2.kpUndist(0), feat2.kpUndist(1));
    }

    const cv::Point2d principalPoint(intri->PrincipalPoint()(0), intri->PrincipalPoint()(1));
    const double focalLength = intri->Focal();
    cv::Mat essential_matrix;
    cv::Mat E = findEssentialMat(ll, rr, focalLength, principalPoint, cv::RANSAC, 0.99, 1.0);
    cv::Mat rot, trans;
    cv::recoverPose(E, ll, rr, rot, trans, focalLength, principalPoint);

    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    for (int i = 0; i < 3; i++) {
        t(i) = trans.at<double>(i, 0);
        for (int j = 0; j < 3; j++) {
            R(i, j) = rot.at<double>(i, j);
        }
    }
    t = t / t.norm();

    // opencv computes pose from 1 to 2, thus return its inverse
    return Sophus::SE3d(Sophus::makeRotationMatrix(R), t).inverse();

#else
    auto bearingVec1 = ComputeBeringVec<0>(featPairVec, intri);
    auto bearingVec2 = ComputeBeringVec<1>(featPairVec, intri);

    opengv::relative_pose::CentralRelativeAdapter adapter(bearingVec1, bearingVec2);
    // Create a RelativePoseSac problem and Ransac
    // Set algorithm to NISTER, STEWENIUS, SEVENPT, or EIGHTPT
    opengv::sac::Ransac<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> ransac;
    std::shared_ptr<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> probPtr(
        new opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem(
            adapter,
            opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem::STEWENIUS));
    ransac.sac_model_ = probPtr;
    ransac.threshold_ = intri.ImagePlaneToCameraPlaneError(1.0);
    ransac.max_iterations_ = 50;

    auto res = ransac.computeModel(0);
    if (!res) {
        return {};
    }

    Eigen::Matrix3d rot = ransac.model_coefficients_.block<3, 3>(0, 0);
    Eigen::Vector3d trans = ransac.model_coefficients_.col(3);
    trans /= trans.norm();

    // opengv computes pose from 2 to 1
    return Sophus::SE3d(Sophus::makeRotationMatrix(rot), trans);

#endif
}

#undef USE_OPENCV_TRIANGULATE
#undef USE_OPENCV_REL_POSE
}  // namespace ns_ikalibr