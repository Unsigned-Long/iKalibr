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

#include "core/lidar_odometer.h"
#include "sensor/lidar.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

LiDAROdometer::LiDAROdometer(float ndtResolution, int threads)
    : _ndtResolution(ndtResolution),
      _threads(threads),
      _map(nullptr),
      _mapTime(0.0),
      _ndt(new pclomp::NormalDistributionsTransform<IKalibrPoint, IKalibrPoint>),
      _initialized(false) {
    // init the ndt omp object
    _ndt->setResolution(ndtResolution);
    _ndt->setNumThreads(threads);
    _ndt->setNeighborhoodSearchMethod(pclomp::DIRECT7);
    _ndt->setTransformationEpsilon(1E-3);
    _ndt->setStepSize(0.01);
    _ndt->setMaximumIterations(50);
}

LiDAROdometer::Ptr LiDAROdometer::Create(float ndtResolution, int threads) {
    return std::make_shared<LiDAROdometer>(ndtResolution, threads);
}

ns_ctraj::Posed LiDAROdometer::FeedFrame(const LiDARFrame::Ptr &frame,
                                         const Eigen::Matrix4d &predCurToLast,
                                         bool updateMap) {
    ns_ctraj::Posed curLtoM;
    if (!_initialized) {
        // identity
        curLtoM = ns_ctraj::Posed(frame->GetTimestamp());

        // create map
        _map = boost::make_shared<IKalibrPointCloud>();
        _mapTime = frame->GetTimestamp();

        // here the pose id identity
        _initialized = true;

    } else {
        // down sample
        IKalibrPointCloud::Ptr filterCloud(new IKalibrPointCloud());
        DownSampleCloud(frame->GetScan(), filterCloud, 0.5);
        _ndt->setInputSource(filterCloud);

        // organize the pred pose from cur frame to map

        Eigen::Matrix4d predCurLtoM = this->_poseSeq.back().se3().matrix() * predCurToLast;
        IKalibrPointCloud::Ptr outputCloud(new IKalibrPointCloud());
        _ndt->align(*outputCloud, predCurLtoM.cast<float>());

        // get pose
        Eigen::Matrix4d pose = _ndt->getFinalTransformation().cast<double>();
        curLtoM = ns_ctraj::Posed::FromT(pose, frame->GetTimestamp());
    }

    if (updateMap && CheckKeyFrame(curLtoM)) {
        UpdateMap(frame, curLtoM);
        _keyFrameIdx.push_back(_frames.size());
    }

    _poseSeq.push_back(curLtoM);
    _frames.push_back(frame);

    return curLtoM;
}

bool LiDAROdometer::CheckKeyFrame(const ns_ctraj::Posed &LtoM) {
    static Eigen::Vector3d lastPos(0.0, 0.0, 0.0);
    static Eigen::Vector3d lastYPR(0.0, 0.0, 0.0);

    Eigen::Vector3d curPos = LtoM.t;
    double posDist = (curPos - lastPos).norm();

    // get current rotMat, ypr
    Eigen::Vector3d curYPR = RotMatToYPR(LtoM.so3.matrix());
    Eigen::Vector3d deltaAngle = curYPR - lastYPR;
    for (int i = 0; i < 3; i++) {
        deltaAngle(i) = NormalizeAngle(deltaAngle(i));
    }
    deltaAngle = deltaAngle.cwiseAbs();

    if (_frames.empty() || posDist > 0.2 || deltaAngle(0) > 5.0 || deltaAngle(1) > 5.0 ||
        deltaAngle(2) > 5.0) {
        // update state
        lastPos = curPos;
        lastYPR = curYPR;
        return true;
    }
    return false;
}

void LiDAROdometer::UpdateMap(const LiDARFrame::Ptr &frame, const ns_ctraj::Posed &LtoM) {
    // update the first map frame using all points after this program is fine
    if (_frames.empty()) {
        // copy the frame point cloud to the map
        *_map += *frame->GetScan();
    } else {
        // down sample
        IKalibrPointCloud::Ptr filteredCloud(new IKalibrPointCloud);
        DownSampleCloud(frame->GetScan(), filteredCloud, _ndtResolution);

        // transform
        IKalibrPointCloud::Ptr transformCloud(new IKalibrPointCloud);
        pcl::transformPointCloud(*filteredCloud, *transformCloud,
                                 LtoM.se3().matrix().cast<float>());

        *_map += *transformCloud;
    }

    // set the target point cloud
    _ndt->setInputTarget(_map);
}

void LiDAROdometer::DownSampleCloud(const IKalibrPointCloud::Ptr &inCloud,
                                    const IKalibrPointCloud::Ptr &outCloud,
                                    float leafSize) {
    pcl::VoxelGrid<IKalibrPoint> filter;
    filter.setInputCloud(inCloud);
    filter.setLeafSize(leafSize, leafSize, leafSize);
    filter.filter(*outCloud);
}

std::size_t LiDAROdometer::KeyFrameSize() const { return _keyFrameIdx.size(); }

std::size_t LiDAROdometer::FrameSize() const { return _frames.size(); }

const std::vector<size_t> &LiDAROdometer::GetKeyFrameIdxVec() const { return _keyFrameIdx; }

const std::vector<ns_ctraj::Posed> &LiDAROdometer::GetOdomPoseVec() const { return _poseSeq; }

const IKalibrPointCloud::Ptr &LiDAROdometer::GetMap() const { return _map; }

const std::vector<LiDARFrame::Ptr> &LiDAROdometer::GetFramesVec() const { return _frames; }

const pclomp::NormalDistributionsTransform<IKalibrPoint, IKalibrPoint>::Ptr &LiDAROdometer::GetNdt()
    const {
    return _ndt;
}

double LiDAROdometer::GetMapTime() const { return _mapTime; }
}  // namespace ns_ikalibr