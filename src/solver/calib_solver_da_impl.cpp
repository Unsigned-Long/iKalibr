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

#include "calib/calib_data_manager.h"
#include "calib/calib_param_manager.h"
#include "core/optical_flow_trace.h"
#include "core/pts_association.h"
#include "core/scan_undistortion.h"
#include "core/visual_reproj_association.h"
#include "factor/data_correspondence.h"
#include "pcl/common/transforms.h"
#include "pcl/filters/random_sample.h"
#include "pcl/filters/voxel_grid.h"
#include "solver/calib_solver.h"
#include "spdlog/spdlog.h"
#include "util/cloud_define.hpp"
#include "util/tqdm.h"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

std::tuple<IKalibrPointCloud::Ptr, std::map<std::string, std::vector<LiDARFrame::Ptr>>>
CalibSolver::BuildGlobalMapOfLiDAR() const {
    if (!Configor::IsLiDARIntegrated()) {
        return {};
    }

    // -------------------------------------------
    // undisto frames and build map in world frame
    // -------------------------------------------
    auto undistHelper = ScanUndistortion::Create(_splines, _parMagr);

    std::map<std::string, std::vector<LiDARFrame::Ptr>> undistFrames;
    IKalibrPointCloud::Ptr mapCloud(new IKalibrPointCloud);
    for (const auto &[topic, data] : _dataMagr->GetLiDARMeasurements()) {
        spdlog::info("undistort scans for lidar '{}'...", topic);
        undistFrames[topic] =
            undistHelper->UndistortToRef(data, topic, ScanUndistortion::Option::ALL);

        spdlog::info("marge scans from lidar '{}' to map...", topic);
        for (auto &frame : undistFrames.at(topic)) {
            if (frame == nullptr) {
                continue;
            }
            *mapCloud += *frame->GetScan();
        }
    }

    IKalibrPointCloud::Ptr newMapCloud(new IKalibrPointCloud);
    std::vector<int> index;
    pcl::removeNaNFromPointCloud(*mapCloud, *newMapCloud, index);

    return {newMapCloud, undistFrames};
}

IKalibrPointCloud::Ptr CalibSolver::BuildGlobalMapOfRadar() const {
    if (!Configor::IsRadarIntegrated() || GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
        return {};
    }

    IKalibrPointCloud::Ptr radarCloud(new IKalibrPointCloud);

    for (const auto &[topic, data] : _dataMagr->GetRadarMeasurements()) {
        const double TO_RjToBr = _parMagr->TEMPORAL.TO_RjToBr.at(topic);
        IKalibrPointCloud::Ptr curRadarCloud(new IKalibrPointCloud);
        curRadarCloud->reserve(data.size() * data.front()->GetTargets().size());
        for (const auto &ary : data) {
            for (const auto &frame : ary->GetTargets()) {
                auto SE3_CurRjToW = CurRjToW(frame->GetTimestamp(), topic);
                if (SE3_CurRjToW == std::nullopt) {
                    continue;
                }
                Eigen::Vector3d p = *SE3_CurRjToW * frame->GetTargetXYZ();
                IKalibrPoint p2;
                p2.timestamp = frame->GetTimestamp() + TO_RjToBr;
                p2.x = static_cast<float>(p(0));
                p2.y = static_cast<float>(p(1));
                p2.z = static_cast<float>(p(2));
                curRadarCloud->push_back(p2);
            }
        }
        *radarCloud += *curRadarCloud;
    }

    // ---------------------------
    // down sample the radar cloud
    // ---------------------------
    pcl::VoxelGrid<IKalibrPoint> filter;
    filter.setInputCloud(radarCloud);
    auto size = static_cast<float>(Configor::Prior::MapDownSample * 2.0);
    filter.setLeafSize(size, size, size);

    IKalibrPointCloud::Ptr radarCloudSampled(new IKalibrPointCloud);
    filter.filter(*radarCloudSampled);
    _viewer->AddStarMarkCloud(radarCloudSampled, Viewer::VIEW_MAP);

    return radarCloud;
}

ColorPointCloud::Ptr CalibSolver::BuildGlobalColorMapOfRGBD(const std::string &topic) const {
    if (!Configor::IsRGBDIntegrated() || GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
        return {};
    }

    spdlog::info("build global map for rgbd '{}'...", topic);

    auto intri = _parMagr->INTRI.RGBD.at(topic);

    ColorPointCloud::Ptr map(new ColorPointCloud);

    auto bar = std::make_shared<tqdm>();
    const auto &frames = _dataMagr->GetRGBDMeasurements(topic);
    for (int i = 0; i < static_cast<int>(frames.size()); ++i) {
        bar->progress(i, static_cast<int>(frames.size()));
        const auto &frame = frames.at(i);

        // transformation
        auto SE3_CurDnToW = CurDnToW(frame->GetTimestamp(), topic);
        if (SE3_CurDnToW == std::nullopt) {
            continue;
        }

        // save points to 'cloud'
        ColorPointCloud::Ptr cloud = frame->CreatePointCloud(intri, 0.1f, 8.0f);
        if (cloud == nullptr) {
            continue;
        }

        // down sample the map cloud
        ColorPointCloud::Ptr cloudDownSampled(new ColorPointCloud);

        pcl::RandomSample<ColorPoint> filter;
        filter.setInputCloud(cloud);
        filter.setSample(1E4);
        filter.filter(*cloudDownSampled);

        // transform cloud to map coordinate frame
        ColorPointCloud::Ptr cloudTransformed(new ColorPointCloud);
        pcl::transformPointCloud(*cloudDownSampled, *cloudTransformed,
                                 SE3_CurDnToW->matrix().cast<float>());

        *map += *cloudTransformed;
    }
    bar->finish();

    return map;
}

std::tuple<IKalibrPointCloud::Ptr,
           // scans in global frame
           std::map<std::string, std::vector<IKalibrPointCloud::Ptr>>,
           // scans in local frame
           std::map<std::string, std::vector<IKalibrPointCloud::Ptr>>>
CalibSolver::BuildGlobalMapOfRGBD() const {
    if (!Configor::IsRGBDIntegrated() || GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
        return {};
    }
    IKalibrPointCloud::Ptr globalMap(new IKalibrPointCloud);
    std::map<std::string, std::vector<IKalibrPointCloud::Ptr>> scanInGFrame;
    std::map<std::string, std::vector<IKalibrPointCloud::Ptr>> scanInLFrame;

    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        spdlog::info("build global map for rgbd '{}'...", topic);

        auto intri = _parMagr->INTRI.RGBD.at(topic);
        const double rsExpFactor =
            CameraModel::RSCameraExposureFactor(EnumCast::stringToEnum<CameraModelType>(
                Configor::DataStream::RGBDTopics.at(topic).Type));
        const double readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);

        const auto &frames = _dataMagr->GetRGBDMeasurements(topic);
        IKalibrPointCloud::Ptr map(new IKalibrPointCloud);
        auto &curScanInGFrame = scanInGFrame[topic];
        curScanInGFrame.reserve(frames.size());
        auto &curScanInLFrame = scanInLFrame[topic];
        curScanInLFrame.reserve(frames.size());

        auto bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(frames.size()); ++i) {
            bar->progress(i, static_cast<int>(frames.size()));
            const auto &frame = frames.at(i);

            // transformation
            auto SE3_CurDnToW = CurDnToW(frame->GetTimestamp(), topic);
            if (SE3_CurDnToW == std::nullopt) {
                continue;
            }

            // save points to 'cloud'
            IKalibrPointCloud::Ptr cloud =
                frame->CreatePointCloud(rsExpFactor, readout, intri, 0.1f, 8.0f);
            if (cloud == nullptr) {
                continue;
            }

            // down sample the map cloud
            IKalibrPointCloud::Ptr cloudDownSampled(new IKalibrPointCloud);

            pcl::RandomSample<IKalibrPoint> filter;
            filter.setInputCloud(cloud);
            filter.setSample(1E4);
            filter.filter(*cloudDownSampled);

            // transform cloud to map coordinate frame
            IKalibrPointCloud::Ptr cloudTransformed(new IKalibrPointCloud);
            pcl::transformPointCloud(*cloudDownSampled, *cloudTransformed,
                                     SE3_CurDnToW->matrix().cast<float>());

            curScanInLFrame.push_back(cloudDownSampled);
            curScanInGFrame.push_back(cloudTransformed);
            *map += *cloudTransformed;
        }
        bar->finish();
        *globalMap += *map;
    }

    return {globalMap, scanInGFrame, scanInLFrame};
}

std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> CalibSolver::DataAssociationForLiDARs(
    const IKalibrPointCloud::Ptr &map,
    const std::map<std::string, std::vector<LiDARFrame::Ptr>> &undistFrames,
    int ptsCountInEachScan) const {
    if (!Configor::IsLiDARIntegrated()) {
        return {};
    }

    // ---------------------------------
    // Step 1: down sample the map cloud
    // ---------------------------------
    pcl::VoxelGrid<IKalibrPoint> filter;
    filter.setInputCloud(map);
    auto size = static_cast<float>(Configor::Prior::MapDownSample);
    filter.setLeafSize(size, size, size);

    IKalibrPointCloud::Ptr mapDownSampled(new IKalibrPointCloud);
    filter.filter(*mapDownSampled);

    _viewer->AddAlignedCloud(mapDownSampled, Viewer::VIEW_MAP, -_parMagr->GRAVITY.cast<float>(),
                             2.0f);

    // ------------------------------------------------
    // Step 2: perform data association for each frames
    // ------------------------------------------------
    auto associator = PointToSurfelAssociator::Create(
        // we use the dense map to create data associator for high-perform point-to-surfel search
        map, Configor::Prior::LiDARDataAssociate::MapResolution,
        Configor::Prior::LiDARDataAssociate::MapDepthLevels);
    auto condition = PointToSurfelCondition();
    _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
    _viewer->AddSurfelMap(associator->GetSurfelMap(), condition, Viewer::VIEW_ASSOCIATION);
    _viewer->AddCloud(mapDownSampled, Viewer::VIEW_ASSOCIATION,
                      ns_viewer::Colour::Black().WithAlpha(0.2f), 2.0f);

    // deconstruction
    mapDownSampled.reset();

    std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> pointToSurfel;

    std::size_t count = 0;
    std::shared_ptr<tqdm> bar;
    for (const auto &[topic, framesInMap] : undistFrames) {
        const auto &rawFrames = _dataMagr->GetLiDARMeasurements(topic);
        spdlog::info("perform point to surfel association for lidar '{}'...", topic);

        // for each scan, we keep 'ptsCountInEachScan' point to surfel corrs
        pointToSurfel[topic] = {};
        auto &curPointToSurfel = pointToSurfel.at(topic);
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(framesInMap.size()); ++i) {
            bar->progress(i, static_cast<int>(framesInMap.size()));

            if (framesInMap.at(i) == nullptr || rawFrames.at(i) == nullptr) {
                continue;
            }

            auto ptsVec = associator->Association(framesInMap.at(i)->GetScan(),
                                                  rawFrames.at(i)->GetScan(), condition);

            curPointToSurfel.insert(curPointToSurfel.end(), ptsVec.cbegin(), ptsVec.cend());
        }
        bar->finish();

        // downsample
        int expectCount = ptsCountInEachScan * static_cast<int>(rawFrames.size());
        if (static_cast<int>(curPointToSurfel.size()) > expectCount) {
            std::map<ufo::map::Node, std::vector<PointToSurfelCorr::Ptr>> nodes;
            for (const auto &corr : curPointToSurfel) {
                nodes[corr->node].push_back(corr);
            }
            std::size_t numEachNode = expectCount / nodes.size() + 1;
            curPointToSurfel.clear();

            // uniform sampling
            std::default_random_engine engine(
                std::chrono::steady_clock::now().time_since_epoch().count());
            for (const auto &[node, corrs] : nodes) {
                auto newCorrs =
                    SamplingWoutReplace2(engine, corrs, std::min(corrs.size(), numEachNode));
                curPointToSurfel.insert(curPointToSurfel.end(), newCorrs.cbegin(), newCorrs.cend());
            }
        }
        count += curPointToSurfel.size();
    }
    spdlog::info("total point to surfel count for LiDARs: {}", count);
    _viewer->AddPointToSurfel(associator->GetSurfelMap(), pointToSurfel, Viewer::VIEW_ASSOCIATION);

    return pointToSurfel;
}

std::map<std::string, std::vector<PointToSurfelCorrPtr>> CalibSolver::DataAssociationForRGBDs(
    const IKalibrPointCloud::Ptr &map,
    const std::map<std::string, std::vector<IKalibrPointCloud::Ptr>> &scanInGFrame,
    const std::map<std::string, std::vector<IKalibrPointCloud::Ptr>> &scanInLFrame,
    int ptsCountInEachScan) const {
    if (!Configor::IsRGBDIntegrated() || GetScaleType() != TimeDeriv::LIN_POS_SPLINE) {
        return {};
    }

    // ---------------------------------
    // Step 1: down sample the map cloud
    // ---------------------------------
    pcl::VoxelGrid<IKalibrPoint> filter;
    filter.setInputCloud(map);
    auto size = static_cast<float>(Configor::Prior::MapDownSample);
    filter.setLeafSize(size, size, size);

    IKalibrPointCloud::Ptr mapDownSampled(new IKalibrPointCloud);
    filter.filter(*mapDownSampled);

    _viewer->AddAlignedCloud(mapDownSampled, Viewer::VIEW_MAP, -_parMagr->GRAVITY.cast<float>(),
                             2.0f);

    // ------------------------------------------------
    // Step 2: perform data association for each frames
    // ------------------------------------------------
    auto associator = PointToSurfelAssociator::Create(
        // we use the dense map to create data associator for high-perform point-to-surfel search
        map, Configor::Prior::LiDARDataAssociate::MapResolution,
        Configor::Prior::LiDARDataAssociate::MapDepthLevels);
    auto condition = PointToSurfelCondition(Configor::Prior::LiDARDataAssociate::PointToSurfelMax,
                                            // for rgbds we relax the conditions
                                            Configor::Prior::LiDARDataAssociate::QueryDepthMin + 1,
                                            // for rgbds we relax the conditions
                                            Configor::Prior::LiDARDataAssociate::QueryDepthMax + 1,
                                            Configor::Prior::LiDARDataAssociate::SurfelPointMin,
                                            Configor::Prior::LiDARDataAssociate::PlanarityMin);
    _viewer->ClearViewer(Viewer::VIEW_ASSOCIATION);
    _viewer->AddSurfelMap(associator->GetSurfelMap(), condition, Viewer::VIEW_ASSOCIATION);
    _viewer->AddCloud(mapDownSampled, Viewer::VIEW_ASSOCIATION,
                      ns_viewer::Colour::Black().WithAlpha(0.2f), 2.0f);

    // deconstruction
    mapDownSampled.reset();

    std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> pointToSurfel;

    std::size_t count = 0;
    std::shared_ptr<tqdm> bar;
    for (const auto &[topic, framesInMap] : scanInGFrame) {
        spdlog::info("perform point to surfel association for rgbd '{}'...", topic);

        // for each scan, we keep 'ptsCountInEachScan' point to surfel corrs
        const auto &rawFrames = scanInLFrame.at(topic);
        auto &curPointToSurfel = pointToSurfel[topic];
        bar = std::make_shared<tqdm>();
        for (int i = 0; i < static_cast<int>(framesInMap.size()); ++i) {
            bar->progress(i, static_cast<int>(framesInMap.size()));

            if (framesInMap.at(i) == nullptr || rawFrames.at(i) == nullptr) {
                continue;
            }

            auto ptsVec = associator->Association(framesInMap.at(i), rawFrames.at(i), condition);

            curPointToSurfel.insert(curPointToSurfel.end(), ptsVec.cbegin(), ptsVec.cend());
        }
        bar->finish();

        // downsample
        int expectCount = ptsCountInEachScan * static_cast<int>(rawFrames.size());
        if (static_cast<int>(curPointToSurfel.size()) > expectCount) {
            std::map<ufo::map::Node, std::vector<PointToSurfelCorr::Ptr>> nodes;
            for (const auto &corr : curPointToSurfel) {
                nodes[corr->node].push_back(corr);
            }
            std::size_t numEachNode = expectCount / nodes.size() + 1;
            curPointToSurfel.clear();

            // uniform sampling
            std::default_random_engine engine(
                std::chrono::steady_clock::now().time_since_epoch().count());
            for (const auto &[node, corrs] : nodes) {
                auto newCorrs =
                    SamplingWoutReplace2(engine, corrs, std::min(corrs.size(), numEachNode));
                curPointToSurfel.insert(curPointToSurfel.end(), newCorrs.cbegin(), newCorrs.cend());
            }
        }
        count += curPointToSurfel.size();
    }
    spdlog::info("total point to surfel count for RGBDs: {}", count);
    _viewer->AddPointToSurfel(associator->GetSurfelMap(), pointToSurfel, Viewer::VIEW_ASSOCIATION);

    return pointToSurfel;
}

std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>>
CalibSolver::DataAssociationForPosCameras() const {
    if (!Configor::IsPosCameraIntegrated()) {
        return {};
    }

    std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>> corrs;
    for (const auto &[topic, sfmData] : _dataMagr->GetSfMData()) {
        spdlog::info("performing visual reprojection data association for camera '{}'...", topic);
        corrs[topic] =
            VisualReProjAssociator::Create(EnumCast::stringToEnum<CameraModelType>(
                                               Configor::DataStream::CameraTopics.at(topic).Type))
                ->Association(*sfmData, _parMagr->INTRI.Camera.at(topic));
        _viewer->AddVeta(sfmData, Viewer::VIEW_MAP);
        spdlog::info("visual reprojection sequences for '{}': {}", topic, corrs.at(topic).size());
    }
    return corrs;
}

std::vector<OpticalFlowTripleTrace::Ptr> CalibSolver::CreateOpticalFlowTrace(
    const std::list<RotOnlyVisualOdometer::FeatTrackingInfo> &trackInfoList, int trackThd) {
    std::vector<OpticalFlowTripleTrace::Ptr> dynamics;
    for (const auto &trackInfo : trackInfoList) {
        for (const auto &[id, info] : trackInfo) {
            // for each tracked feature
            if (static_cast<int>(info.size()) < std::max(trackThd, 3)) {
                continue;
            }
            // store in groups of three
            auto groupSize = info.size() - 2;
            auto iter0 = info.cbegin();
            for (int i = 0; i < static_cast<int>(groupSize); ++i, ++iter0) {
                auto iter1 = iter0, iter2 = iter0;
                std::advance(iter1, 1), std::advance(iter2, 2);
                std::array<std::pair<CameraFramePtr, Eigen::Vector2d>, 3> movement;
                movement.at(0) = {
                    iter0->first,
                    Eigen::Vector2d(iter0->second.undistorted.x, iter0->second.undistorted.y),
                };
                movement.at(1) = {
                    iter1->first,
                    Eigen::Vector2d(iter1->second.undistorted.x, iter1->second.undistorted.y),
                };
                movement.at(2) = {
                    iter2->first,
                    Eigen::Vector2d(iter2->second.undistorted.x, iter2->second.undistorted.y),
                };
                dynamics.push_back(OpticalFlowTripleTrace::Create(movement));
            }
        }
    }
    return dynamics;
}

std::map<std::string, std::vector<OpticalFlowCorr::Ptr>> CalibSolver::DataAssociationForRGBDs(
    bool estDepth) {
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

    std::map<std::string, std::vector<OpticalFlowCorr::Ptr>> corrs;

    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        const auto &traceVec = _dataMagr->GetVisualOpticalFlowTrace(topic);
        spdlog::info("perform data association for RGBD camera '{}'...", topic);
        const auto &intri = _parMagr->INTRI.RGBD.at(topic);
        const double fx = intri->intri->FocalX(), fy = intri->intri->FocalY();
        const double cx = intri->intri->PrincipalPoint()(0), cy = intri->intri->PrincipalPoint()(1);

        const auto &rsExposureFactor =
            CameraModel::RSCameraExposureFactor(EnumCast::stringToEnum<CameraModelType>(
                Configor::DataStream::RGBDTopics.at(topic).Type));

        const double readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);
        const double TO_DnToBr = _parMagr->TEMPORAL.TO_DnToBr.at(topic);

        const auto SO3_DnToBr = _parMagr->EXTRI.SO3_DnToBr.at(topic);
        Sophus::SO3d SO3_BrToDn = SO3_DnToBr.inverse();
        Eigen::Vector3d POS_DnInBr = _parMagr->EXTRI.POS_DnInBr.at(topic);

        auto &curCorrs = corrs[topic];
        curCorrs.reserve(traceVec.size());

        int estDepthCount = 0;

        for (const auto &trace : traceVec) {
            auto corr = trace->CreateOpticalFlowCorr(rsExposureFactor, intri);

            double timeByBr = corr->MidPointTime(readout) + TO_DnToBr;
            if (!so3Spline.TimeStampInRange(timeByBr) || !scaleSpline.TimeStampInRange(timeByBr)) {
                continue;
            }

            Eigen::Vector3d LIN_VEL_BrToBr0InBr0;
            switch (GetScaleType()) {
                case TimeDeriv::LIN_ACCE_SPLINE:
                    // this would not happen
                    continue;
                case TimeDeriv::LIN_VEL_SPLINE:
                    LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<0>(timeByBr);
                    break;
                case TimeDeriv::LIN_POS_SPLINE:
                    LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<1>(timeByBr);
                    break;
            }

            Sophus::SO3d SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
            Eigen::Vector3d ANG_VEL_BrToBr0InBr = so3Spline.VelocityBody(timeByBr);
            Eigen::Vector3d ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * ANG_VEL_BrToBr0InBr;
            Eigen::Vector3d ANG_VEL_DnToBr0InDn = SO3_BrToDn * ANG_VEL_BrToBr0InBr;

            Eigen::Vector3d LIN_VEL_DnToBr0InBr0 =
                -Sophus::SO3d::hat(SO3_BrToBr0 * POS_DnInBr) * ANG_VEL_BrToBr0InBr0 +
                LIN_VEL_BrToBr0InBr0;

            corr->withDepthObservability =
                // cond. 1: the rgbd camera is moving
                (LIN_VEL_DnToBr0InBr0.norm() > 0.3 /* m/sed */) &&
                // cond. 2: this feature is moving (not necessary but involved here)
                (corr->MidPointVel(readout).norm() > 100.0 /* pixels/sed */);

            // for those tracked features, but without depth information.
            // if the depth is ready to estimate, we assign rough depth for them if they have depth
            // observability. otherwise, if the depth is not ready to be estimated, we do not
            // introduce them to estimator.
            if (auto actualDepth = intri->ActualDepth(corr->depth); actualDepth < 1E-3) {
                if (!estDepth || !corr->withDepthObservability) {
                    // do not  introduce it to estimator.
                    continue;
                }

                Eigen::Vector3d LIN_VEL_DnToBr0InDn =
                    SO3_BrToDn * SO3_BrToBr0.inverse() * LIN_VEL_DnToBr0InBr0;

                Eigen::Matrix<double, 2, 3> subAMat, subBMat;
                OpticalFlowCorr::SubMats<double>(&fx, &fy, &cx, &cy, corr->MidPoint(), &subAMat,
                                                 &subBMat);

                Eigen::Vector2d lVec = subAMat * LIN_VEL_DnToBr0InDn;
                Eigen::Vector2d bMat = corr->MidPointVel(readout) - subBMat * ANG_VEL_DnToBr0InDn;
                Eigen::Vector1d HMat = bMat.transpose() * bMat;
                double estimatedDepth = (HMat.inverse() * bMat.transpose() * lVec)(0, 0);

                // spdlog::info(
                //     "raw depth: {:.3f}, actual depth: {:.3f}, est depth: {:.3f}, new raw depth: "
                //     "{:.3f}",
                //     corr->depth, actualDepth, estDepth, newRawDepth);

                if (estimatedDepth < 1E-3) {
                    // depth is negative, do not  introduce it to estimator.
                    continue;
                }
                corr->depth = estimatedDepth;
                corr->invDepth = 1.0 / corr->depth;
                ++estDepthCount;
            }
            const double pVelNorm = corr->MidPointVel(readout).norm();
            corr->weight = pVelNorm / (pVelNorm + Configor::Prior::LossForOpticalFlowFactor);

            curCorrs.push_back(corr);
        }

        auto dObvCount = std::count_if(curCorrs.begin(), curCorrs.end(), [&](const auto &item) {
            return item->withDepthObservability;
        });
        spdlog::info(
            "total correspondences count for rgbd '{}': {}, with roughly estimated initial depth "
            "count: {}, with depth observability: {}",
            topic, curCorrs.size(), estDepthCount, dObvCount);

        std::default_random_engine engine(
            std::chrono::system_clock::now().time_since_epoch().count());
        constexpr int CorrCountPerFrame = 50;
        auto desiredCount = CorrCountPerFrame * _dataMagr->GetRGBDMeasurements(topic).size();
        if (desiredCount < curCorrs.size()) {
            curCorrs = SamplingWoutReplace2(engine, curCorrs, desiredCount);
            spdlog::info("total correspondences count for rgbd '{}' after down sampled: {}", topic,
                         curCorrs.size());
        }
    }

    // add veta for visualization
    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        const auto &intri = _parMagr->INTRI.RGBD.at(topic);

        auto veta =
            CreateVetaFromOpticalFlow(topic, corrs.at(topic), intri->intri, &CalibSolver::CurDnToW);
        if (veta != nullptr) {
            DownsampleVeta(veta, 10000, Configor::DataStream::RGBDTopics.at(topic).TrackLengthMin);
            // we do not show the pose
            _viewer->AddVeta(veta, Viewer::VIEW_MAP, {}, ns_viewer::Entity::GetUniqueColour());
        }
    }

    return corrs;
}

std::map<std::string, std::vector<OpticalFlowCorr::Ptr>> CalibSolver::DataAssociationForVelCameras()
    const {
    const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);

    std::map<std::string, std::vector<OpticalFlowCorr::Ptr>> corrs;

    for (const auto &[topic, _] : Configor::DataStream::VelCameraTopics()) {
        const auto &traceVec = _dataMagr->GetVisualOpticalFlowTrace(topic);
        spdlog::info("perform data association for camera '{}'...", topic);
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        const double fx = intri->FocalX(), fy = intri->FocalY();
        const double cx = intri->PrincipalPoint()(0), cy = intri->PrincipalPoint()(1);

        const auto &rsExposureFactor =
            CameraModel::RSCameraExposureFactor(EnumCast::stringToEnum<CameraModelType>(
                Configor::DataStream::CameraTopics.at(topic).Type));

        const double readout = _parMagr->TEMPORAL.RS_READOUT.at(topic);
        const double TO_CmToBr = _parMagr->TEMPORAL.TO_CmToBr.at(topic);

        const auto SO3_CmToBr = _parMagr->EXTRI.SO3_CmToBr.at(topic);
        Sophus::SO3d SO3_BrToCm = SO3_CmToBr.inverse();
        Eigen::Vector3d POS_CmInBr = _parMagr->EXTRI.POS_CmInBr.at(topic);

        auto &curCorrs = corrs[topic];
        curCorrs.reserve(traceVec.size());

        for (const auto &dynamic : traceVec) {
            auto corr = dynamic->CreateOpticalFlowCorr(rsExposureFactor);

            double timeByBr = corr->MidPointTime(readout) + TO_CmToBr;
            if (!so3Spline.TimeStampInRange(timeByBr) || !scaleSpline.TimeStampInRange(timeByBr)) {
                continue;
            }

            Eigen::Vector3d LIN_VEL_BrToBr0InBr0;
            switch (GetScaleType()) {
                case TimeDeriv::LIN_ACCE_SPLINE:
                    // this would not happen
                    continue;
                case TimeDeriv::LIN_VEL_SPLINE:
                    LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<0>(timeByBr);
                    break;
                case TimeDeriv::LIN_POS_SPLINE:
                    LIN_VEL_BrToBr0InBr0 = scaleSpline.Evaluate<1>(timeByBr);
                    break;
            }

            Sophus::SO3d SO3_BrToBr0 = so3Spline.Evaluate(timeByBr);
            Eigen::Vector3d ANG_VEL_BrToBr0InBr = so3Spline.VelocityBody(timeByBr);
            Eigen::Vector3d ANG_VEL_BrToBr0InBr0 = SO3_BrToBr0 * ANG_VEL_BrToBr0InBr;
            Eigen::Vector3d ANG_VEL_CmToBr0InCm = SO3_BrToCm * ANG_VEL_BrToBr0InBr;

            Eigen::Vector3d LIN_VEL_CmToBr0InBr0 =
                -Sophus::SO3d::hat(SO3_BrToBr0 * POS_CmInBr) * ANG_VEL_BrToBr0InBr0 +
                LIN_VEL_BrToBr0InBr0;

            corr->withDepthObservability =
                // cond. 1: the rgbd camera is moving
                (LIN_VEL_CmToBr0InBr0.norm() > 0.3 /* m/sed */) &&
                // cond. 2: this feature is moving (not necessary but involved here)
                (corr->MidPointVel(readout).norm() > 100.0 /* pixels/sed */);

            if (!corr->withDepthObservability) {
                // do not  introduce it to estimator.
                continue;
            }

            Eigen::Vector3d LIN_VEL_CmToBr0InCm =
                SO3_BrToCm * SO3_BrToBr0.inverse() * LIN_VEL_CmToBr0InBr0;

            Eigen::Matrix<double, 2, 3> subAMat, subBMat;
            OpticalFlowCorr::SubMats<double>(&fx, &fy, &cx, &cy, corr->MidPoint(), &subAMat,
                                             &subBMat);

            Eigen::Vector2d lVec = subAMat * LIN_VEL_CmToBr0InCm;
            Eigen::Vector2d bMat = corr->MidPointVel(readout) - subBMat * ANG_VEL_CmToBr0InCm;
            Eigen::Vector1d HMat = bMat.transpose() * bMat;
            double estDepth = (HMat.inverse() * bMat.transpose() * lVec)(0, 0);

            if (estDepth < 1E-3) {
                // depth is negative, do not  introduce it to estimator.
                continue;
            }
            corr->depth = estDepth;
            corr->invDepth = 1.0 / corr->depth;
            const double pVelNorm = corr->MidPointVel(readout).norm();
            corr->weight = pVelNorm / (pVelNorm + Configor::Prior::LossForOpticalFlowFactor);

            curCorrs.push_back(corr);
        }

        spdlog::info("total correspondences count for camera '{}': {}", topic, curCorrs.size());

        std::default_random_engine engine(
            std::chrono::system_clock::now().time_since_epoch().count());
        constexpr int CorrCountPerFrame = 50;
        auto desiredCount = CorrCountPerFrame * _dataMagr->GetCameraMeasurements(topic).size();
        if (desiredCount < curCorrs.size()) {
            curCorrs = SamplingWoutReplace2(engine, curCorrs, desiredCount);
            spdlog::info("total correspondences count for camera '{}' after down sampled: {}",
                         topic, curCorrs.size());
        }
    }

    // add veta from pixel dynamics
    for (const auto &[topic, _] : Configor::DataStream::VelCameraTopics()) {
        const auto &intri = _parMagr->INTRI.Camera.at(topic);

        auto veta =
            CreateVetaFromOpticalFlow(topic, corrs.at(topic), intri, &CalibSolver::CurCmToW);
        if (veta != nullptr) {
            DownsampleVeta(veta, 10000,
                           Configor::DataStream::CameraTopics.at(topic).TrackLengthMin);
            // we do not show the pose
            _viewer->AddVeta(veta, Viewer::VIEW_MAP, {}, ns_viewer::Entity::GetUniqueColour());
        }
    }
    return corrs;
}

}  // namespace ns_ikalibr