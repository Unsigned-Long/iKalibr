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

#ifndef IKALIBR_VISUAL_COLORIZED_CLOUD_MAP_H
#define IKALIBR_VISUAL_COLORIZED_CLOUD_MAP_H

#include "ctraj/core/spline_bundle.h"
#include "config/configor.h"
#include "veta/veta.h"
#include "veta/camera/pinhole.h"
#include "opencv2/core.hpp"
#include "util/cloud_define.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CalibParamManager;
using CalibParamManagerPtr = std::shared_ptr<CalibParamManager>;
struct CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;
struct RGBDFrame;
using RGBDFramePtr = std::shared_ptr<RGBDFrame>;

class ColorizedCloudMap {
public:
    using Ptr = std::shared_ptr<ColorizedCloudMap>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;

private:
    std::string _topic;
    const std::vector<CameraFramePtr> &_frames;
    ns_veta::Veta::Ptr _veta;
    SplineBundleType::Ptr _splines;
    CalibParamManagerPtr _parMagr;

    ns_veta::PinholeIntrinsic::Ptr _intri;

    Sophus::SE3d SE3_CmToBr;
    double TO_CmToBr;

    std::map<ns_veta::IndexT, std::pair<std::optional<Sophus::SE3d>, cv::Mat>> _viewIdToFrame;

public:
    ColorizedCloudMap(const std::vector<CameraFramePtr> &frames,
                      ns_veta::Veta::Ptr veta,
                      SplineBundleType::Ptr splines,
                      ns_veta::PinholeIntrinsic::Ptr intri,
                      const Sophus::SE3d &SE3_SenToBr,
                      const double &TO_SenToBr);

    static Ptr CreateForCameras(const std::string &topic,
                                const std::vector<CameraFramePtr> &frames,
                                const ns_veta::Veta::Ptr &veta,
                                const SplineBundleType::Ptr &splines,
                                const CalibParamManagerPtr &parMagr);

    static Ptr CreateForRGBDs(const std::string &topic,
                              const std::vector<CameraFramePtr> &frames,
                              const ns_veta::Veta::Ptr &veta,
                              const SplineBundleType::Ptr &splines,
                              const CalibParamManagerPtr &parMagr);

    ColorPointCloud::Ptr Colorize(const IKalibrPointCloud::Ptr &cloudMap, int K = 5);

protected:
    std::optional<Sophus::SE3d> CurCmToW(double timeByCm);

    static cv::Mat DrawPoint(const cv::Mat &img,
                             const Eigen::Vector2i &p,
                             const cv::Scalar &color = cv::Scalar(0, 255, 0));
};
}  // namespace ns_ikalibr
#endif  // IKALIBR_VISUAL_COLORIZED_CLOUD_MAP_H
