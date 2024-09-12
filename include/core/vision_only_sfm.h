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

#ifndef IKALIBR_VISION_ONLY_SFM_H
#define IKALIBR_VISION_ONLY_SFM_H

#include "config/configor.h"

#include "ctraj/core/spline_bundle.h"
#include "tiny-viewer/entity/entity.h"

#include "veta/veta.h"
#include "veta/camera/pinhole.h"

#include "boost/geometry/geometries/polygon.hpp"
#include "boost/geometry/geometries/point_xy.hpp"

#include "opengv/sac/Ransac.hpp"
#include "opengv/absolute_pose/CentralAbsoluteAdapter.hpp"
#include "opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp"

#include "opencv2/core.hpp"
#include "spdlog/spdlog.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;
struct CalibParamManager;
using CalibParamManagerPtr = std::shared_ptr<CalibParamManager>;
struct Viewer;
using ViewerPtr = std::shared_ptr<Viewer>;

namespace bg = boost::geometry;
typedef boost::geometry::model::d2::point_xy<double, boost::geometry::cs::cartesian> point_2d;
typedef boost::geometry::model::polygon<point_2d> polygon_2d;

struct SfMFeature {
public:
    ns_veta::IndexT id{};
    // just for visualization
    cv::Point2f kp;
    ns_veta::Vec2d kpUndist;

public:
    SfMFeature(ns_veta::IndexT id, cv::Point2f kp, ns_veta::Vec2d kpUndist)
        : id(id),
          kp(std::move(kp)),
          kpUndist(std::move(kpUndist)) {}

    SfMFeature() = default;
};

using SfMFeatureVec = std::vector<SfMFeature>;
using SfMFeaturePair = std::pair<SfMFeature, SfMFeature>;
using SfMFeaturePairVec = std::vector<SfMFeaturePair>;

using IndexPair = std::pair<ns_veta::IndexT, ns_veta::IndexT>;
using PolyPair = std::pair<polygon_2d, polygon_2d>;

struct SfMFeaturePairInfo {
public:
    IndexPair viewId;
    SfMFeaturePairVec featPairVec;
    PolyPair poly;

public:
    SfMFeaturePairInfo(ns_veta::IndexT view1Id,
                       ns_veta::IndexT view2Id,
                       SfMFeaturePairVec featPairVec,
                       const polygon_2d &poly1,
                       const polygon_2d &poly2)
        : viewId(view1Id, view2Id),
          featPairVec(std::move(featPairVec)),
          poly(poly1, poly2) {}

    [[nodiscard]] std::optional<Sophus::SE3d> RecoverRelativePose(
        const ns_veta::PinholeIntrinsic::Ptr &intri) const;
};

inline opengv::bearingVectors_t ComputeBeringVec(const SfMFeatureVec &feat,
                                                 const ns_veta::PinholeIntrinsic::Ptr &intri) {
    opengv::bearingVectors_t bearingVec(feat.size());
    for (int i = 0; i < static_cast<int>(feat.size()); ++i) {
        ns_veta::Vec2d pCam = intri->ImgToCam(feat.at(i).kpUndist);
        bearingVec.at(i) = ns_veta::Vec3d(pCam(0), pCam(1), 1.0).normalized();
    }
    return bearingVec;
}

template <int PairIndex>
opengv::bearingVectors_t ComputeBeringVec(const SfMFeaturePairVec &featPairVec,
                                          const ns_veta::PinholeIntrinsic::Ptr &intri) {
    opengv::bearingVectors_t bearingVec(featPairVec.size());
    for (int i = 0; i < static_cast<int>(featPairVec.size()); ++i) {
        const auto &feat = std::get<PairIndex>(featPairVec.at(i));
        ns_veta::Vec2d pInCam = intri->ImgToCam(feat.kpUndist);
        bearingVec.at(i) = ns_veta::Vec3d(pInCam(0), pInCam(1), 1.0).normalized();
    }
    return bearingVec;
}

template <int PairIndex>
opengv::bearingVectors_t ComputeBeringVec(const SfMFeaturePairVec &featPairVec,
                                          const ns_veta::PinholeIntrinsic::Ptr &intri,
                                          const std::vector<int> &consideredIdx) {
    opengv::bearingVectors_t bearingVec(consideredIdx.size());
    for (int i = 0; i < static_cast<int>(consideredIdx.size()); ++i) {
        auto idx = consideredIdx.at(i);
        const auto &feat = std::get<PairIndex>(featPairVec.at(idx));
        ns_veta::Vec2d pInCam = intri->ImgToCam(feat.kpUndist);
        bearingVec.at(i) = ns_veta::Vec3d(pInCam(0), pInCam(1), 1.0).normalized();
    }
    return bearingVec;
}

class VisionOnlySfM {
public:
    using Ptr = std::shared_ptr<VisionOnlySfM>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;
    using So3SplineType = SplineBundleType::So3SplineType;

    // index, key points, descriptors, undistorted key points
    using FeaturePack = std::tuple<std::vector<ns_veta::IndexT>,
                                   std::vector<cv::KeyPoint>,
                                   cv::Mat,
                                   std::vector<ns_veta::Vec2d>>;

private:
    std::string _topic;
    std::vector<CameraFramePtr> _frames;

    CalibParamManagerPtr _parMagr;
    const So3SplineType &_so3Spline;
    ns_veta::IndexT _lmLabeler;
    ViewerPtr _viewer;

    // generated in process
    ns_veta::Veta::Ptr _veta;
    ns_veta::PinholeIntrinsic::Ptr _intri;
    std::map<ns_veta::IndexT, CameraFramePtr> _frameBackup;
    std::map<IndexPair, SfMFeaturePairInfo> _matchRes;
    std::map<ns_veta::IndexT, std::map<ns_veta::IndexT, ns_veta::IndexT>> _viewFeatLM;

    // for visualization
    std::vector<ns_viewer::Entity::Ptr> _viewCubes;
    std::map<ns_veta::IndexT, Sophus::SE3d> _viewCubePoses;

    // pixel
    constexpr static double REPROJ_THD_PNP = 1.0;
    constexpr static double REPROJ_THD_TRI = 1.0;
    // degree
    constexpr static double PARALLAX_THD_TRI = 2.0;

public:
    VisionOnlySfM(std::string topic,
                  const std::vector<CameraFramePtr> &frames,
                  CalibParamManagerPtr parMagr,
                  const So3SplineType &so3Spline,
                  ViewerPtr viewer);

    static Ptr Create(const std::string &topic,
                      const std::vector<CameraFramePtr> &frames,
                      const CalibParamManagerPtr &parMagr,
                      const So3SplineType &so3Spline,
                      const ViewerPtr &viewer);

    bool PreProcess();

    bool StructureFromMotion();

    [[nodiscard]] const std::map<IndexPair, SfMFeaturePairInfo> &GetMatchRes() const;

    std::set<IndexPair> FindCovisibility(double covThd = 0.2);

protected:
    std::optional<ns_veta::Posed> ComputeCamRotations(const CameraFramePtr &frame);

    std::map<CameraFramePtr, std::pair<polygon_2d, polygon_2d>> FindCovisibility(
        const CameraFramePtr &refFrame, const std::set<IndexPair> &ignore, double covThd = 0.2);

    std::optional<std::pair<polygon_2d, polygon_2d>> IntersectionArea(const cv::Mat &i1,
                                                                      const cv::Mat &i2,
                                                                      const Sophus::SO3d &SO3_2To1);

    IndexPair InitStructure();

    cv::Mat DrawMatchResult(const IndexPair &viewIdPair);

    inline ns_veta::IndexT GenNewLMId() { return ++_lmLabeler; }

    void IncrementalSfM(const IndexPair &sViewIdxPair);

    [[nodiscard]] inline ns_veta::IndexT FindLMByViewFeatId(ns_veta::IndexT viewId,
                                                            ns_veta::IndexT featId) const {
        const auto &feats = _viewFeatLM.at(viewId);
        auto iter = feats.find(featId);
        if (iter == feats.cend()) {
            return ns_veta::UndefinedIndexT;
        } else {
            return iter->second;
        }
    }

    [[nodiscard]] inline bool InsertViewFeatLM(ns_veta::IndexT viewId,
                                               ns_veta::IndexT featId,
                                               ns_veta::IndexT lmId) {
        auto &feats = _viewFeatLM.at(viewId);
        return feats.insert({featId, lmId}).second;
    }

    void CreateViewCubes();

    void DrawMatchesInViewer(
        const ns_viewer::Colour &color = ns_viewer::Colour::Black(),
        const std::set<IndexPair> &special = {},
        const ns_viewer::Colour &specialColor = ns_viewer::Colour::Green()) const;

    std::map<IndexPair, Eigen::Vector3d> Triangulate(const SfMFeaturePairInfo &featPair,
                                                     const Eigen::Matrix3d &R12,
                                                     const Eigen::Vector3d &t12,
                                                     double reProjThd = REPROJ_THD_TRI,
                                                     double parallaxDegThd = PARALLAX_THD_TRI);

    void InsertTriangulateLM(const SfMFeaturePairInfo &featPair,
                             const std::map<IndexPair, Eigen::Vector3d> &lms,
                             const ns_veta::Posed &curLMToW);

    template <int RefViewIndex, int TarViewIndex>
    bool SolveAbsolutePose(SfMFeaturePairInfo &viewPair, double reProjThd = REPROJ_THD_PNP) {
        //            constexpr int RefViewIndex = 0, TarViewIndex = 1;

        assert(RefViewIndex != TarViewIndex);
        assert(RefViewIndex == 0 || RefViewIndex == 1);
        assert(TarViewIndex == 0 || TarViewIndex == 1);

        auto refViewId = std::get<RefViewIndex>(viewPair.viewId);
        auto tarViewId = std::get<TarViewIndex>(viewPair.viewId);
        const auto &featPairVec = viewPair.featPairVec;

        // rough initial pose for target view
        auto T_RefToW = _veta->poses.at(refViewId);

        std::vector<int> idxRecorder;
        idxRecorder.reserve(featPairVec.size());
        for (int i = 0; i < static_cast<int>(featPairVec.size()); ++i) {
            const auto &feats = featPairVec.at(i);
            const auto &refFeat = std::get<RefViewIndex>(feats);
            auto lmId = FindLMByViewFeatId(refViewId, refFeat.id);

            // this feat in the ref view has not been triangulated
            if (lmId == ns_veta::UndefinedIndexT) {
                continue;
            }

            idxRecorder.push_back(i);
        }

        if (idxRecorder.size() < 15) {
            return false;
        }

#define USE_OPENCV_ABS_POSE 0
#if USE_OPENCV_ABS_POSE
        std::vector<cv::Point3f> pts3d(idxRecorder.size());
        std::vector<cv::Point2f> pts2d(idxRecorder.size());
        for (int i = 0; i < static_cast<int>(idxRecorder.size()); ++i) {
            auto idx = idxRecorder.at(i);
            const auto &feats = featPairVec.at(idx);
            const auto &refFeat = std::get<RefViewIndex>(feats);

            auto lmId = FindLMByViewFeatId(refViewId, refFeat.id);

            const Eigen::Vector2d &tarFeatUndist = std::get<TarViewIndex>(feats).kpUndist;
            const Eigen::Vector3d &lmX = _veta->structure.at(lmId).X;

            pts2d.at(i) = cv::Point2f((float)tarFeatUndist(0), (float)tarFeatUndist(1));
            pts3d.at(i) = cv::Point3f((float)lmX(0), (float)lmX(1), (float)lmX(2));
        }
#else
        opengv::bearingVectors_t bearingVec(idxRecorder.size());
        opengv::points_t pointVec(idxRecorder.size());
        for (int i = 0; i < static_cast<int>(idxRecorder.size()); ++i) {
            auto idx = idxRecorder.at(i);
            const auto &feats = featPairVec.at(idx);
            const auto &refFeat = std::get<RefViewIndex>(feats);

            auto lmId = FindLMByViewFeatId(refViewId, refFeat.id);

            const Eigen::Vector2d &tarFeatUndist = std::get<TarViewIndex>(feats).kpUndist;
            const Eigen::Vector3d &lmX = _veta->structure.at(lmId).X;

            ns_veta::Vec2d pInCam = _intri->ImgToCam(tarFeatUndist);
            bearingVec.at(i) = ns_veta::Vec3d(pInCam(0), pInCam(1), 1.0).normalized();
            pointVec.at(i) = lmX;
        }
#endif

#if USE_OPENCV_ABS_POSE
        auto T_WToRef = T_RefToW.Inverse();
        cv::Mat K = (cv::Mat_<double>(3, 3) << _intri->FocalX(), 0, _intri->PrincipalPoint()(0), 0,
                     _intri->FocalY(), _intri->PrincipalPoint()(1), 0, 0, 1);
        // initial rotation from world to target camera
        cv::Mat rVec = [&T_WToRef]() {
            cv::Mat rMat, rVec;
            cv::eigen2cv(T_WToRef.Rotation().matrix(), rMat);
            cv::Rodrigues(rMat, rVec);
            return rVec;
        }();
        // initial translation from world to target camera
        cv::Mat tVec = [&T_WToRef]() {
            cv::Mat t;
            cv::eigen2cv(T_WToRef.Translation(), t);
            return t;
        }();
        // the estimated pose is from world to camera
        std::vector<uchar> cvInliers;
        auto res = cv::solvePnPRansac(pts3d, pts2d, K, cv::Mat(), rVec, tVec, true, 50,
                                      (float)reProjThd, 0.99, cvInliers);
        if (!res) {
            return false;
        }
        cv::Mat RMat;
        cv::Rodrigues(rVec, RMat);

        // from world to camera
        Eigen::Matrix3d rot;
        Eigen::Vector3d pos;
        cv::cv2eigen(RMat, rot);
        cv::cv2eigen(tVec, pos);
        // from camera to world
        pos = -rot.transpose() * pos;
        rot = rot.transpose();

#else
        // the points are in world frame, thus the solved pose are [from camera to world]
        opengv::absolute_pose::CentralAbsoluteAdapter adapter(
            bearingVec, pointVec, T_RefToW.Translation(), T_RefToW.Rotation().matrix());

        // Create an AbsolutePoseSac problem and Ransac
        // The method can be set to KNEIP, GAO or EPNP
        opengv::sac::Ransac<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
        std::shared_ptr<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> probPtr(
            new opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem(
                adapter, opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem::KNEIP));

        ransac.sac_model_ = probPtr;
        ransac.threshold_ = 1.0 - cos(atan(_intri->ImagePlaneToCameraPlaneError(reProjThd)));
        ransac.max_iterations_ = 50;
        auto res = ransac.computeModel(0);
        if (!res) {
            return false;
        }

        Eigen::Matrix3d rot = ransac.model_coefficients_.block<3, 3>(0, 0);
        Eigen::Vector3d pos = ransac.model_coefficients_.col(3);
#endif

        Sophus::SE3d SE3_TarToW(Sophus::makeRotationMatrix(rot), pos);

        _veta->poses.insert(
            {tarViewId, ns_veta::Posed(SE3_TarToW.so3(), SE3_TarToW.translation())});

        std::set<int> inliers;
#if USE_OPENCV_ABS_POSE
        for (int i = 0; i < static_cast<int>(cvInliers.size()); ++i) {
            if (cvInliers.at(i)) {
                inliers.insert(idxRecorder.at(i));
            }
        }
#else
        for (int inlierIdx : ransac.inliers_) {
            inliers.insert(idxRecorder.at(inlierIdx));
        }
#endif
        spdlog::info("reference view: {}, target view: {}, inlier size: {}/{}", refViewId,
                     tarViewId, inliers.size(), idxRecorder.size());

        // for these involved matches
        for (const auto &idx : inliers) {
            const auto &feats = featPairVec.at(idx);
            const auto &refFeat = std::get<RefViewIndex>(feats);
            const auto &tarFeat = std::get<TarViewIndex>(feats);

            // such landmark must exist
            auto lm = FindLMByViewFeatId(refViewId, refFeat.id);
            if (_veta->structure.at(lm)
                    .obs.insert({tarViewId, {tarFeat.kpUndist, tarFeat.id}})
                    .second) {
                if (!InsertViewFeatLM(tarViewId, tarFeat.id, lm)) {
                    spdlog::warn(
                        "[VisionOnlySfM::SolveAbsolutePose] view feat has connected to landmark "
                        "but not inserted to 'ViewFeatLM'");
                }
            }
        }

        // target to reference
        return true;
#undef USE_OPENCV_ABS_POSE
    }

    void BatchOptimization();

protected:
    static polygon_2d BufferPolygon(const polygon_2d &polygon, double bufferDistance);

    static std::optional<polygon_2d> ProjPolygon(const cv::Mat &img,
                                                 const Sophus::SO3d &so3,
                                                 const ns_veta::PinholeIntrinsic::Ptr &intri);

    static void DrawProjPolygon(cv::Mat &img,
                                const polygon_2d &poly,
                                const cv::Scalar &fill,
                                const cv::Scalar &border);

    static FeaturePack FindInBorderOnes(const FeaturePack &input, const polygon_2d &poly);

    static cv::Mat DrawInBorderFeatMatch(const cv::Mat &img,
                                         const polygon_2d &poly,
                                         const FeaturePack &allFeats,
                                         const FeaturePack &ibFeats);

    static std::pair<SfMFeatureVec, SfMFeatureVec> MatchFeatures(const FeaturePack &feat1,
                                                                 const FeaturePack &feat2);

    static std::vector<int> RejectOutliers(const SfMFeatureVec &feat1,
                                           const SfMFeatureVec &feat2,
                                           const ns_veta::PinholeIntrinsic::Ptr &intri);

    static bool IsGraphConnect(const std::vector<IndexPair> &graph);
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_VISION_ONLY_SFM_H
