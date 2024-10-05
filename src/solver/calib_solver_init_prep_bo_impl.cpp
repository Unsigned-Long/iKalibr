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

#include "solver/calib_solver.h"
#include "core/lidar_odometer.h"
#include "calib/calib_data_manager.h"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void CalibSolver::InitPrepBatchOpt() const {
    /**
     * align initialized states to gravity direction
     */
    spdlog::info("aligning all states to gravity direction...");
    AlignStatesToGravity();
    _viewer->UpdateSplineViewer();

    /**
     * transform the veta to world frame if Cameras are integrated
     */
    for (const auto &[camTopic, poseSeq] : _initAsset->sfmPoseSeq) {
        auto SE3_Cm0ToBr0 = this->CurCmToW(poseSeq.front().timeStamp, camTopic);

        if (SE3_Cm0ToBr0 == std::nullopt) {
            throw Status(Status::CRITICAL, "map time of '{}' is out of time range of splines!",
                         camTopic);
        }
        PerformTransformForVeta(_dataMagr->GetSfMData(camTopic),
                                ns_veta::Posed(SE3_Cm0ToBr0->so3(), SE3_Cm0ToBr0->translation()),
                                1.0);
    }

    /**
     * build map and undisto frames if LiDARs are integrated
     */
    spdlog::info("build global map and undisto lidar frames in world...");
    auto &globalMap = _initAsset->globalMap;
    globalMap = IKalibrPointCloud::Ptr(new IKalibrPointCloud);
    auto &undistFramesInMap = _initAsset->undistFramesInMap;
    for (const auto &[lidarTopic, odometer] : _initAsset->lidarOdometers) {
        auto SE3_Lk0ToBr0 = this->CurLkToW(odometer->GetMapTime(), lidarTopic);

        if (SE3_Lk0ToBr0 == std::nullopt) {
            throw Status(Status::CRITICAL, "map time of '{}' is out of time range of splines!",
                         lidarTopic);
        }

        const auto &curUndistFramesInScan = _initAsset->undistFramesInScan.at(lidarTopic);

        undistFramesInMap[lidarTopic] = std::vector<LiDARFrame::Ptr>(curUndistFramesInScan.size());
        auto &curUndistFramesInMap = undistFramesInMap.at(lidarTopic);

        const auto &poseSeq = odometer->GetOdomPoseVec();

        for (int i = 0; i < static_cast<int>(poseSeq.size()); ++i) {
            const auto &undistoScan = curUndistFramesInScan.at(i);
            if (undistoScan == nullptr) {
                curUndistFramesInMap.at(i) = nullptr;
            } else {
                const auto &SE3_ScanToLk0 = poseSeq.at(i);
                Sophus::SE3d SE3_ScanToBr0 = *SE3_Lk0ToBr0 * SE3_ScanToLk0.se3();
                IKalibrPointCloud::Ptr scanInBr0(new IKalibrPointCloud);
                pcl::transformPointCloud(*undistoScan->GetScan(), *scanInBr0,
                                         SE3_ScanToBr0.matrix().cast<float>());
                curUndistFramesInMap.at(i) =
                    LiDARFrame::Create(undistoScan->GetTimestamp(), scanInBr0);

                *globalMap += *scanInBr0;
            }
        }
    }
}
}  // namespace ns_ikalibr