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

#ifndef IKALIBR_LIDAR_ODOMETER_H
#define IKALIBR_LIDAR_ODOMETER_H

#include "util/utils.h"
#include "util/cloud_define.hpp"
#include "ctraj/core/pose.hpp"
#include "pclomp/ndt_omp.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct LiDARFrame;
using LiDARFramePtr = std::shared_ptr<LiDARFrame>;

class LiDAROdometer {
public:
    using Ptr = std::shared_ptr<LiDAROdometer>;

private:
    float _ndtResolution;
    int _threads;

    // the key frame index in the '_frames'
    std::vector<std::size_t> _keyFrameIdx;
    std::vector<LiDARFramePtr> _frames;

    // the global map
    IKalibrPointCloud::Ptr _map;
    double _mapTime;

    // ndt
    pclomp::NormalDistributionsTransform<IKalibrPoint, IKalibrPoint>::Ptr _ndt;

    bool _initialized;

    // pose sequence
    std::vector<ns_ctraj::Posed> _poseSeq;

public:
    LiDAROdometer(float ndtResolution, int threads);

    static LiDAROdometer::Ptr Create(float ndtResolution, int threads);

    ns_ctraj::Posed FeedFrame(const LiDARFramePtr &frame,
                              const Eigen::Matrix4d &predCurToLast = Eigen::Matrix4d::Identity(),
                              bool updateMap = true);

    [[nodiscard]] std::size_t KeyFrameSize() const;

    [[nodiscard]] std::size_t FrameSize() const;

    // getters
    [[nodiscard]] const std::vector<size_t> &GetKeyFrameIdxVec() const;

    [[nodiscard]] const std::vector<ns_ctraj::Posed> &GetOdomPoseVec() const;

    [[nodiscard]] const IKalibrPointCloud::Ptr &GetMap() const;

    [[nodiscard]] const std::vector<LiDARFramePtr> &GetFramesVec() const;

    [[nodiscard]] const pclomp::NormalDistributionsTransform<IKalibrPoint, IKalibrPoint>::Ptr &
    GetNdt() const;

    [[nodiscard]] double GetMapTime() const;

protected:
    bool CheckKeyFrame(const ns_ctraj::Posed &LtoM);

    void UpdateMap(const LiDARFramePtr &frame, const ns_ctraj::Posed &LtoM);

    static void DownSampleCloud(const IKalibrPointCloud::Ptr &inCloud,
                                const IKalibrPointCloud::Ptr &outCloud,
                                float leafSize);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_LIDAR_ODOMETER_H
