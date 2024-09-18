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

#ifndef IKALIBR_VISUAL_ANG_VEL_DRAWER_H
#define IKALIBR_VISUAL_ANG_VEL_DRAWER_H

#include "ctraj/core/spline_bundle.h"
#include "config/configor.h"
#include "veta/veta.h"
#include "veta/camera/pinhole.h"
#include "opencv2/core.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CalibParamManager;
using CalibParamManagerPtr = std::shared_ptr<CalibParamManager>;
class CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;
struct OpticalFlowCorr;
using OpticalFlowCorrPtr = std::shared_ptr<OpticalFlowCorr>;

class VisualAngVelDrawer {
public:
    using Ptr = std::shared_ptr<VisualAngVelDrawer>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;

private:
    std::string _topic;
    ns_veta::Veta::Ptr _veta;
    SplineBundleType::Ptr _splines;

    ns_veta::PinholeIntrinsic::Ptr _intri;

    Sophus::SE3d SE3_CmToBr;
    double TO_CmToBr;

public:
    VisualAngVelDrawer(std::string topic,
                       ns_veta::Veta::Ptr veta,
                       SplineBundleType::Ptr splines,
                       const CalibParamManagerPtr &parMagr);

    static Ptr Create(const std::string &topic,
                      const ns_veta::Veta::Ptr &veta,
                      const SplineBundleType::Ptr &splines,
                      const CalibParamManagerPtr &parMagr);

    cv::Mat CreateAngVelImg(const CameraFramePtr &frame, float scale = 1.0f);
};

class VisualOpticalFlowAngVelDrawer {
public:
    using Ptr = std::shared_ptr<VisualOpticalFlowAngVelDrawer>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;

private:
    // frame id, correspondences
    std::map<ns_veta::IndexT, std::vector<OpticalFlowCorrPtr>> _velCorrs;
    SplineBundleType::Ptr _splines;

    ns_veta::PinholeIntrinsic::Ptr _intri;

    Sophus::SE3d SE3_SenToBr;
    double TO_SenToBr;
    Eigen::Vector3d GRAVITY;

public:
    VisualOpticalFlowAngVelDrawer(const std::vector<OpticalFlowCorrPtr> &corrs,
                                  SplineBundleType::Ptr splines,
                                  ns_veta::PinholeIntrinsic::Ptr intri,
                                  const Sophus::SE3d &SE3_SenToBr,
                                  const double &TO_SenToBr);

    static Ptr CreateDrawerForRGBDs(const std::string &topic,
                                    const std::vector<OpticalFlowCorrPtr> &corrs,
                                    const SplineBundleType::Ptr &splines,
                                    const CalibParamManagerPtr &parMagr);

    static Ptr CreateDrawerForVelCameras(const std::string &topic,
                                         const std::vector<OpticalFlowCorrPtr> &corrs,
                                         const SplineBundleType::Ptr &splines,
                                         const CalibParamManagerPtr &parMagr);

    cv::Mat CreateAngVelImg(const CameraFramePtr &frame, float scale = 0.3f);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_VISUAL_ANG_VEL_DRAWER_H
