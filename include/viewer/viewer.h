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

#ifndef IKALIBR_VIEWER_H
#define IKALIBR_VIEWER_H

#include "config/configor.h"
#include "tiny-viewer/core/multi_viewer.h"
#include "tiny-viewer/object/aligned_cloud.hpp"
#include "cereal/types/polymorphic.hpp"
#include "ctraj/core/spline_bundle.h"
#include "util/cloud_define.hpp"
#include "veta/veta.h"
#include "ufo/map/surfel_map.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CalibParamManager;
using CalibParamManagerPtr = std::shared_ptr<CalibParamManager>;
struct PointToSurfelCondition;
struct PointToSurfelCorr;
using PointToSurfelCorrPtr = std::shared_ptr<PointToSurfelCorr>;
struct RGBDFrame;
using RGBDFramePtr = std::shared_ptr<RGBDFrame>;
struct RGBDIntrinsics;
using RGBDIntrinsicsPtr = std::shared_ptr<RGBDIntrinsics>;

class Viewer : public ns_viewer::MultiViewer {
public:
    using Ptr = std::shared_ptr<Viewer>;
    using Parent = ns_viewer::MultiViewer;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;

public:
    const static std::string VIEW_SENSORS, VIEW_SPLINE, VIEW_MAP, VIEW_ASSOCIATION;

private:
    CalibParamManagerPtr _parMagr;
    SplineBundleType::Ptr _splines;

    std::map<std::string, std::vector<std::size_t>> _entities;

public:
    explicit Viewer(CalibParamManagerPtr parMagr, SplineBundleType::Ptr splines);

    static Ptr Create(const CalibParamManagerPtr &parMagr, const SplineBundleType::Ptr &splines);

    Viewer &FillEmptyViews(const std::string &objPath);

    Viewer &UpdateSensorViewer();

    Viewer &UpdateSplineViewer(double dt = 0.005);

    Viewer &AddAlignedCloud(const IKalibrPointCloud::Ptr &cloud,
                            const std::string &view,
                            const Eigen::Vector3f &dir = {0, 0, 1},
                            float size = DefaultPointSize);

    Viewer &AddCloud(const IKalibrPointCloud::Ptr &cloud,
                     const std::string &view,
                     const ns_viewer::Colour &color,
                     float size = DefaultPointSize);

    Viewer &AddStarMarkCloud(const IKalibrPointCloud::Ptr &cloud,
                             const std::string &view,
                             float size = DefaultLandmarkSize);

    Viewer &AddCloud(const IKalibrPointCloud::Ptr &cloud,
                     const std::string &view,
                     float size = DefaultPointSize);

    Viewer &ClearViewer(const std::string &view);

    Viewer &PopBackEntity(const std::string &view);

    Viewer &AddSurfelMap(const ufo::map::SurfelMap &smp,
                         const PointToSurfelCondition &condition,
                         const std::string &view);

    Viewer &AddPointToSurfel(const ufo::map::SurfelMap &smp,
                             const std::map<std::string, std::vector<PointToSurfelCorrPtr>> &corrs,
                             const std::string &view);

    ns_viewer::Entity::Ptr Gravity() const;

    Viewer &AddVeta(const ns_veta::Veta::Ptr &veta,
                    const std::string &view,
                    const std::optional<ns_viewer::Colour> &camColor = ns_viewer::Colour::Blue(),
                    const std::optional<ns_viewer::Colour> &lmColor = {});

    Viewer &AddEntityLocal(const std::vector<ns_viewer::Entity::Ptr> &entities,
                           const std::string &view);

    void SetNewSpline(const SplineBundleType::Ptr &splines);

    Viewer &AddRGBDFrame(const RGBDFramePtr &frame,
                         const RGBDIntrinsicsPtr &intri,
                         const std::string &view,
                         bool trueColor,
                         float size);

protected:
    ns_viewer::MultiViewerConfigor GenViewerConfigor();

    void ZoomInSplineCallBack();

    void ZoomOutSplineCallBack();

    void ZoomInCoordCallBack();

    void ZoomOutCoordCallBack();
};
}  // namespace ns_ikalibr

namespace ns_viewer {
template <>
struct Cloud<IKalibrPoint> : public Entity {
public:
    using Ptr = std::shared_ptr<Cloud>;

protected:
    using PointCloud = IKalibrPointCloud;
    using ColorPointCloud = pcl::PointCloud<pcl::PointXYZRGBA>;
    using PointCloudPtr = PointCloud::Ptr;
    using ColorPointCloudPtr = ColorPointCloud::Ptr;

    ColorPointCloudPtr cloud;

    float size{};

public:
    explicit Cloud(const PointCloudPtr &inputCloud,
                   float size = DefaultPointSize,
                   IntensityMode mode = IntensityMode::PCL_VISUALIZER_LUT_HSV)
        : Entity(),
          cloud(new ColorPointCloud),
          size(size) {
        pcl::visualization::PointCloudColorHandlerGenericField<IKalibrPoint> colorHandler(
            inputCloud, "intensity");
        auto colors = colorHandler.getColor();
        double minmax[2];
        colors->GetRange(minmax);
        vtkSmartPointer<vtkLookupTable> table = GetColormapLUT(mode, minmax);
        // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
        // double rgb[3];
        cloud->resize(inputCloud->size());
        for (int i = 0; i < static_cast<int>(inputCloud->size()); ++i) {
            auto &ip = inputCloud->at(i);
            auto &op = cloud->at(i);
            op.x = ip.x;
            op.y = ip.y;
            op.z = ip.z;
            // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
            // table->GetColor(ip.intensity, rgb);
            // op.r = static_cast<std::uint8_t>(rgb[0] * 255.0f);
            // op.g = static_cast<std::uint8_t>(rgb[1] * 255.0f);
            // op.b = static_cast<std::uint8_t>(rgb[2] * 255.0f);
            op.a = static_cast<std::uint8_t>(255.0f);
        }
    }

    explicit Cloud(const PointCloudPtr &inputCloud,
                   const Colour &color,
                   float size = DefaultPointSize)
        : Entity(),
          cloud(new ColorPointCloud),
          size(size) {
        cloud->resize(inputCloud->size());
        for (int i = 0; i < static_cast<int>(inputCloud->size()); ++i) {
            auto &ip = inputCloud->at(i);
            auto &op = cloud->at(i);
            op.x = ip.x;
            op.y = ip.y;
            op.z = ip.z;
            op.r = static_cast<std::uint8_t>(color.r * 255.0f);
            op.g = static_cast<std::uint8_t>(color.g * 255.0f);
            op.b = static_cast<std::uint8_t>(color.b * 255.0f);
            op.a = static_cast<std::uint8_t>(color.a * 255.0f);
        }
    }

    static Ptr Create(const PointCloudPtr &cloud, float size = DefaultPointSize) {
        return std::make_shared<Cloud>(cloud, size);
    }

    static Ptr Create(const PointCloudPtr &cloud,
                      const Colour &color,
                      float size = DefaultPointSize) {
        return std::make_shared<Cloud>(cloud, color, size);
    }

    ~Cloud() override = default;

    void Draw() const override {
        glPointSize(size);
        glBegin(GL_POINTS);
        for (const auto &p : cloud->points) {
            glColor4f(ExpandPCLColor(p));
            glVertex3f(ExpandPCLPointXYZ(p));
        }
        glEnd();
    }

    Cloud()
        : cloud(new ColorPointCloud) {}

    [[nodiscard]] const ColorPointCloudPtr &GetCloud() const { return cloud; }

public:
    template <class Archive>
    void serialize(Archive &archive) {
        Entity::serialize(archive);
        archive(cereal::make_nvp("data", *cloud), CEREAL_NVP(size));
    }
};
}  // namespace ns_viewer
CEREAL_REGISTER_TYPE_WITH_NAME(ns_viewer::Cloud<IKalibrPoint>, "Cloud::IKalibrPoint")
CEREAL_REGISTER_POLYMORPHIC_RELATION(ns_viewer::Entity, ns_viewer::Cloud<IKalibrPoint>)

CEREAL_REGISTER_TYPE_WITH_NAME(ns_viewer::AlignedCloud<IKalibrPoint>, "AlignedCloud::IKalibrPoint")
CEREAL_REGISTER_POLYMORPHIC_RELATION(ns_viewer::Entity, ns_viewer::AlignedCloud<IKalibrPoint>)

#endif  // IKALIBR_VIEWER_H
