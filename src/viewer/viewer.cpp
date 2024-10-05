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

#include "viewer/viewer.h"
#include "util/status.hpp"
#include "core/pts_association.h"
#include "calib/calib_param_manager.h"
#include "sensor/rgbd.h"
#include "factor/data_correspondence.h"
#include "tiny-viewer/object/landmark.h"
#include "tiny-viewer/object/camera.h"
#include "tiny-viewer/core/pose.hpp"
#include "tiny-viewer/entity/arrow.h"
#include "tiny-viewer/entity/line.h"
#include "tiny-viewer/object/surfel.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
// static members
const std::string Viewer::VIEW_SENSORS = "VIEW_SENSORS";
const std::string Viewer::VIEW_SPLINE = "VIEW_SPLINE";
const std::string Viewer::VIEW_MAP = "VIEW_MAP";
const std::string Viewer::VIEW_ASSOCIATION = "VIEW_ASSOCIATION";

ns_ikalibr::Viewer::Viewer(CalibParamManager::Ptr parMagr, SplineBundleType::Ptr splines)
    : Parent(GenViewerConfigor()),
      _parMagr(std::move(parMagr)),
      _splines(std::move(splines)) {
    // create containers for entities
    _entities.insert({VIEW_SENSORS, {}});
    _entities.insert({VIEW_SPLINE, {}});
    _entities.insert({VIEW_MAP, {}});
    _entities.insert({VIEW_ASSOCIATION, {}});

    // run
    this->RunInMultiThread();
}

std::shared_ptr<Viewer> ns_ikalibr::Viewer::Create(const CalibParamManager::Ptr &parMagr,
                                                   const SplineBundleType::Ptr &splines) {
    return std::make_shared<Viewer>(parMagr, splines);
}

Viewer &Viewer::UpdateSensorViewer() {
    ClearViewer(VIEW_SENSORS);
    _entities.at(VIEW_SENSORS) = _parMagr->VisualizationSensors(*this, VIEW_SENSORS);
    return *this;
}

Viewer &Viewer::UpdateSplineViewer(double dt) {
    ClearViewer(VIEW_SPLINE);

    // spline poses
    std::vector<ns_viewer::Entity::Ptr> entities;
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    double minTime = std::max(so3Spline.MinTime(), scaleSpline.MinTime());
    double maxTime = std::min(so3Spline.MaxTime(), scaleSpline.MaxTime());
    for (double t = minTime; t < maxTime;) {
        if (!_splines->TimeInRange(t, so3Spline) || !_splines->TimeInRange(t, scaleSpline)) {
            t += dt;
            continue;
        }

        Sophus::SO3d so3 = so3Spline.Evaluate(t);
        Eigen::Vector3d linScale =
            scaleSpline.Evaluate(t) * Configor::Preference::SplineScaleInViewer;
        // coordinate
        entities.push_back(ns_viewer::Coordinate::Create(
            ns_viewer::Posed(so3.matrix(), linScale).cast<float>(),
            static_cast<float>(Configor::Preference::CoordSScaleInViewer)));
        t += dt;
    }
    const auto &knots = scaleSpline.GetKnots();
    for (const auto &k : knots) {
        entities.push_back(
            ns_viewer::Landmark::Create(k.cast<float>() * Configor::Preference::SplineScaleInViewer,
                                        0.1f, ns_viewer::Colour::Black()));
    }
    for (int i = 0; i < static_cast<int>(knots.size()) - 1; ++i) {
        const int j = i + 1;
        const Eigen::Vector3d ki = knots.at(i) * Configor::Preference::SplineScaleInViewer;
        const Eigen::Vector3d kj = knots.at(j) * Configor::Preference::SplineScaleInViewer;
        entities.push_back(ns_viewer::Line::Create(ki.cast<float>(), kj.cast<float>(), 0.1f,
                                                   ns_viewer::Colour::Black()));
    }
    // gravity
    entities.push_back(Gravity());

    _entities.at(VIEW_SPLINE) = this->AddEntity(entities, VIEW_SPLINE);

    return *this;
}

Viewer &Viewer::AddCloud(const IKalibrPointCloud::Ptr &cloud,
                         const std::string &view,
                         const ns_viewer::Colour &color,
                         float size) {
    AddEntityLocal({ns_viewer::Cloud<IKalibrPoint>::Create(cloud, color, size), Gravity()}, view);
    return *this;
}

Viewer &Viewer::AddStarMarkCloud(const IKalibrPointCloud::Ptr &cloud,
                                 const std::string &view,
                                 float size) {
    PosPointCloud::Ptr posCloud(new PosPointCloud);
    pcl::copyPointCloud(*cloud, *posCloud);
    AddEntityLocal({ns_viewer::Cloud<ns_viewer::Landmark>::Create(posCloud, size), Gravity()},
                   view);
    return *this;
}

Viewer &Viewer::AddCloud(const IKalibrPointCloud::Ptr &cloud, const std::string &view, float size) {
    AddEntityLocal({ns_viewer::Cloud<IKalibrPoint>::Create(cloud, size)}, view);
    return *this;
}

Viewer &Viewer::AddAlignedCloud(const IKalibrPointCloud::Ptr &cloud,
                                const std::string &view,
                                const Eigen::Vector3f &dir,
                                float size) {
    AddEntityLocal({ns_viewer::AlignedCloud<IKalibrPoint>::Create(cloud, dir, size), Gravity()},
                   view);
    return *this;
}

ns_viewer::MultiViewerConfigor Viewer::GenViewerConfigor() {
    ns_viewer::MultiViewerConfigor viewConfig(
        {VIEW_SENSORS, VIEW_SPLINE, VIEW_MAP, VIEW_ASSOCIATION}, "iKalibr");

    viewConfig.grid.at(VIEW_SENSORS).showGrid = false;
    viewConfig.grid.at(VIEW_MAP).showGrid = false;
    viewConfig.grid.at(VIEW_ASSOCIATION).showGrid = false;
    viewConfig.grid.at(VIEW_SENSORS).showIdentityCoord = false;

    viewConfig.WithScreenShotSaveDir(Configor::DataStream::OutputPath);

    viewConfig.window.width = 640 * 2;
    viewConfig.window.height = 480 * 2;

    viewConfig.camera.at(VIEW_SENSORS).initPos = {0.7f, 0.7f, 0.7f};
    viewConfig.camera.at(VIEW_MAP).far = 1000.0;
    viewConfig.camera.at(VIEW_ASSOCIATION).far = 1000.0;

    viewConfig.callBacks.insert({'a', [this]() { ZoomOutSplineCallBack(); }});
    viewConfig.callBacks.insert({'d', [this]() { ZoomInSplineCallBack(); }});
    viewConfig.callBacks.insert({'s', [this]() { ZoomOutCoordCallBack(); }});
    viewConfig.callBacks.insert({'w', [this]() { ZoomInCoordCallBack(); }});

    return viewConfig;
}

Viewer &Viewer::ClearViewer(const std::string &view) {
    this->RemoveEntity(_entities.at(view), view);
    _entities.at(view).clear();
    return *this;
}

void Viewer::ZoomInSplineCallBack() {
    Configor::Preference::SplineScaleInViewer += 0.1;
    UpdateSplineViewer();
}

void Viewer::ZoomOutSplineCallBack() {
    if (Configor::Preference::SplineScaleInViewer > 0.2) {
        Configor::Preference::SplineScaleInViewer -= 0.1;
    }
    UpdateSplineViewer();
}

void Viewer::ZoomInCoordCallBack() {
    Configor::Preference::CoordSScaleInViewer += 0.1;
    UpdateSplineViewer();
}

void Viewer::ZoomOutCoordCallBack() {
    if (Configor::Preference::CoordSScaleInViewer > 0.2) {
        Configor::Preference::CoordSScaleInViewer -= 0.1;
    }
    UpdateSplineViewer();
}

Viewer &Viewer::AddSurfelMap(const ufo::map::SurfelMap &smp,
                             const PointToSurfelCondition &condition,
                             const std::string &view) {
    namespace ufopred = ufo::map::predicate;
    std::vector<ns_viewer::Entity::Ptr> entities;

    auto pred = ufopred::HasSurfel() && ufopred::DepthMin(condition.queryDepthMin) &&
                ufopred::DepthMax(condition.queryDepthMax) &&
                ufopred::NumSurfelPointsMin(condition.surfelPointMin) &&
                ufopred::SurfelPlanarityMin(condition.planarityMin);

    for (const auto &node : smp.query(pred)) {
        // create entities
        auto cen = smp.getNodeCenter(node);
        auto pose =
            ns_viewer::Posef(Eigen::Matrix3f::Identity(), Eigen::Vector3f(cen.x, cen.y, cen.z));
        auto min = smp.getNodeMin(node), max = smp.getNodeMax(node);
        auto cube = ns_viewer::Cube::Create(pose, true, max.x - min.x, max.y - min.y, max.z - min.z,
                                            ns_viewer::Colour::Black().WithAlpha(0.4f));

        entities.push_back(cube);
    }

    entities.push_back(Gravity());

    AddEntityLocal(entities, view);
    return *this;
}

Viewer &Viewer::AddPointToSurfel(
    const ufo::map::SurfelMap &smp,
    const std::map<std::string, std::vector<PointToSurfelCorrPtr>> &corrs,
    const std::string &view) {
    std::map<ufo::map::Node, std::vector<PointToSurfelCorr::Ptr>> nodes;
    for (const auto &[topic, corrVec] : corrs) {
        for (const auto &corr : corrVec) {
            nodes[corr->node].push_back(corr);
        }
    }

    std::vector<ns_viewer::Entity::Ptr> entities;
    for (const auto &[node, nodeCorrs] : nodes) {
        auto color = ns_viewer::Entity::GetUniqueColour();

        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->reserve(nodeCorrs.size());
        for (const auto &corr : nodeCorrs) {
            pcl::PointXYZ p;
            p.x = static_cast<float>(corr->pInMap(0));
            p.y = static_cast<float>(corr->pInMap(1));
            p.z = static_cast<float>(corr->pInMap(2));
            cloud->push_back(p);
        }
        entities.push_back(ns_viewer::Cloud<pcl::PointXYZ>::Create(cloud, DefaultPointSize, color));

        // create surfels
        auto cen = smp.getNodeCenter(node);
        auto pose =
            ns_viewer::Posef(Eigen::Matrix3f::Identity(), Eigen::Vector3f(cen.x, cen.y, cen.z));
        auto min = smp.getNodeMin(node), max = smp.getNodeMax(node);
        auto cube = ns_viewer::Cube::Create(pose, true, max.x - min.x, max.y - min.y, max.z - min.z,
                                            ns_viewer::Colour::Black().WithAlpha(0.2f));
        auto s = ns_viewer::Surfel::Create(nodeCorrs.front()->surfelInW.cast<float>(), *cube, false,
                                           true, color.WithAlpha(0.2f));
        entities.push_back(s);
    }

    entities.push_back(Gravity());

    AddEntityLocal(entities, view);
    return *this;
}

Viewer &Viewer::PopBackEntity(const std::string &view) {
    auto &curEntities = _entities.at(view);
    if (!curEntities.empty()) {
        this->RemoveEntity(curEntities.back(), view);
        curEntities.pop_back();
    }
    return *this;
}

ns_viewer::Entity::Ptr Viewer::Gravity() const {
    return ns_viewer::Arrow::Create(
        _parMagr->GRAVITY.normalized().cast<float>() * Configor::Preference::SplineScaleInViewer,
        Eigen::Vector3f::Zero(), ns_viewer::Colour::Blue());
}

Viewer &Viewer::AddVeta(const ns_veta::Veta::Ptr &veta,
                        const std::string &view,
                        const std::optional<ns_viewer::Colour> &camColor,
                        const std::optional<ns_viewer::Colour> &lmColor) {
    std::vector<ns_viewer::Entity::Ptr> entities;
    if (camColor != std::nullopt) {
        for (const auto &[viewId, se3] : veta->poses) {
            ns_viewer::Posed pose(se3.Rotation().matrix(), se3.Translation());
            entities.push_back(ns_viewer::Camera::Create(pose.cast<float>(), 0.04, *camColor));
        }
    }

    if (!veta->structure.empty()) {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        cloud->reserve(veta->structure.size());
        for (const auto &[lmId, lm] : veta->structure) {
            pcl::PointXYZRGB p;
            p.x = (float)lm.X(0), p.y = (float)lm.X(1), p.z = (float)lm.X(2);
            if (lmColor != std::nullopt) {
                p.r = static_cast<uchar>(lmColor->r * 255.0f);
                p.g = static_cast<uchar>(lmColor->g * 255.0f);
                p.b = static_cast<uchar>(lmColor->b * 255.0f);
            } else {
                p.r = lm.color(0), p.g = lm.color(1), p.b = lm.color(2);
            }
            cloud->push_back(p);
        }
        entities.push_back(ns_viewer::Cloud<pcl::PointXYZRGB>::Create(cloud, 8.0f));
    }
    AddEntityLocal(entities, view);
    return *this;
}

Viewer &Viewer::AddEntityLocal(const std::vector<ns_viewer::Entity::Ptr> &entities,
                               const std::string &view) {
    auto ids = this->AddEntity(entities, view);
    _entities.at(view).insert(_entities.at(view).end(), ids.cbegin(), ids.cend());
    return *this;
}

void Viewer::SetNewSpline(const SplineBundleType::Ptr &splines) { _splines = splines; }

Viewer &Viewer::FillEmptyViews(const std::string &objPath) {
    std::array<std::map<std::string, bool>, 6> occupy;
    std::array<bool, 6> senIntegrated{};
    // imus
    occupy[0] = {
        {VIEW_SENSORS, true}, {VIEW_SPLINE, true}, {VIEW_MAP, false}, {VIEW_ASSOCIATION, false}};
    senIntegrated[0] = true;

    // radars
    occupy[1] = {
        {VIEW_SENSORS, true}, {VIEW_SPLINE, true}, {VIEW_MAP, false}, {VIEW_ASSOCIATION, false}};
    senIntegrated[1] = Configor::IsRadarIntegrated();

    // pos cameras
    occupy[2] = {
        {VIEW_SENSORS, true}, {VIEW_SPLINE, true}, {VIEW_MAP, true}, {VIEW_ASSOCIATION, true}};
    senIntegrated[2] = Configor::IsPosCameraIntegrated();

    // vel cameras
    occupy[3] = {
        {VIEW_SENSORS, true}, {VIEW_SPLINE, true}, {VIEW_MAP, false}, {VIEW_ASSOCIATION, false}};
    senIntegrated[3] = Configor::IsVelCameraIntegrated();

    // lidars
    occupy[4] = {
        {VIEW_SENSORS, true}, {VIEW_SPLINE, true}, {VIEW_MAP, true}, {VIEW_ASSOCIATION, true}};
    senIntegrated[4] = Configor::IsLiDARIntegrated();

    // rgbds
    occupy[5] = {
        {VIEW_SENSORS, true}, {VIEW_SPLINE, true}, {VIEW_MAP, true}, {VIEW_ASSOCIATION, false}};
    senIntegrated[5] = Configor::IsRGBDIntegrated();

    std::map<std::string, bool> viewOccupy;
    for (int i = 0; i < static_cast<int>(occupy.size()); ++i) {
        if (!senIntegrated.at(i)) {
            continue;
        }
        for (const auto &[view, status] : occupy.at(i)) {
            auto iter = viewOccupy.find(view);
            if (iter == viewOccupy.cend()) {
                viewOccupy[view] = status;
            } else {
                iter->second = iter->second || status;
            }
        }
    }

    if (std::filesystem::exists(objPath)) {
        for (const auto &[view, status] : viewOccupy) {
            if (!status) {
                this->AddObjEntity(objPath, view);
            }
        }
    } else {
        spdlog::warn("can not load models from '{}'!", objPath);
    }

    return *this;
}

Viewer &Viewer::AddRGBDFrame(const RGBDFrame::Ptr &frame,
                             const RGBDIntrinsicsPtr &intri,
                             const std::string &view,
                             bool trueColor,
                             float size) {
    ColorPointCloud ::Ptr cloud(new ColorPointCloud);
    if (trueColor) {
        cloud = frame->CreatePointCloud(intri);
    } else {
        auto cMat = frame->CreateColorDepthMap(intri, false);
        auto dMat = frame->GetDepthImage();
        cloud = RGBDFrame(INVALID_TIME_STAMP, cv::Mat(), cMat, dMat).CreatePointCloud(intri);
    }
    AddEntityLocal({ns_viewer::Cloud<ColorPoint>::Create(cloud, size)}, view);
    return *this;
}
}  // namespace ns_ikalibr