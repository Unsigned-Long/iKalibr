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

#ifndef LIC_CALIB_SCAN_UNDISTORTION_H
#define LIC_CALIB_SCAN_UNDISTORTION_H

#include "config/configor.h"
#include "ctraj/core/spline_bundle.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct LiDARFrame;
using LiDARFramePtr = std::shared_ptr<LiDARFrame>;

struct CalibParamManager;
using CalibParamManagerPtr = std::shared_ptr<CalibParamManager>;

class ScanUndistortion {
public:
    using Ptr = std::shared_ptr<ScanUndistortion>;
    using SplineBundleType = ns_ctraj::SplineBundle<Configor::Prior::SplineOrder>;

public:
    enum class Option : std::uint32_t {
        /**
         * @brief options
         */
        NONE = 1 << 0,
        UNDIST_SO3 = 1 << 1,
        UNDIST_POS = 1 << 2,
        ALL = UNDIST_SO3 | UNDIST_POS
    };

private:
    const SplineBundleType::So3SplineType &_so3Spline;
    const SplineBundleType::RdSplineType &_posSpline;

    CalibParamManagerPtr _parMagr;

public:
    ScanUndistortion(const SplineBundleType::Ptr &splines, CalibParamManagerPtr calibParamManager);

    static ScanUndistortion::Ptr Create(const SplineBundleType::Ptr &splines,
                                        const CalibParamManagerPtr &calibParamManager);

    std::vector<LiDARFramePtr> UndistortToScan(const std::vector<LiDARFramePtr> &data,
                                               const std::string &topic,
                                               Option option);

    std::vector<LiDARFramePtr> UndistortToRef(const std::vector<LiDARFramePtr> &data,
                                              const std::string &topic,
                                              Option option);

protected:
    std::optional<LiDARFramePtr> UndistortToScan(const LiDARFramePtr &lidarFrame,
                                                 const std::string &topic,
                                                 bool correctPos);

    std::optional<LiDARFramePtr> UndistortToRef(const LiDARFramePtr &lidarFrame,
                                                const std::string &topic,
                                                bool correctPos);
};
}  // namespace ns_ikalibr

#endif  // LIC_CALIB_SCAN_UNDISTORTION_H
