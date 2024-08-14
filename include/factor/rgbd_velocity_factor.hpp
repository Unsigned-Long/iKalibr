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

#ifndef IKALIBR_RGBD_VELOCITY_FACTOR_HPP
#define IKALIBR_RGBD_VELOCITY_FACTOR_HPP

#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/ceres.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct RGBDVelocityCorr {
public:
    using Ptr = std::shared_ptr<RGBDVelocityCorr>;
    static constexpr int MID = 1;

public:
    std::array<double, 3> timeAry;
    std::array<double, 3> xDynamicAry;
    std::array<double, 3> yDynamicAry;
    // row / image height
    std::array<double, 3> rDivHAry;
    double depth;
    double rsExpFactor;

public:
    RGBDVelocityCorr(const std::array<double, 3> &timeAry,
                     const std::array<double, 3> &xDynamicAry,
                     const std::array<double, 3> &yDynamicAry,
                     double depth,
                     int imgHeight,
                     double rsExpFactor)
        : timeAry(timeAry),
          xDynamicAry(xDynamicAry),
          yDynamicAry(yDynamicAry),
          rDivHAry(),
          depth(depth),
          rsExpFactor(rsExpFactor) {
        for (int i = 0; i < 3; ++i) {
            rDivHAry[i] = yDynamicAry[i] / (double)imgHeight;
        }
    }

    static Ptr Create(const std::array<double, 3> &timeAry,
                      const std::array<double, 3> &xDynamicAry,
                      const std::array<double, 3> &yDynamicAry,
                      double depth,
                      int imgHeight,
                      double rsExpFactor) {
        return std::make_shared<RGBDVelocityCorr>(timeAry, xDynamicAry, yDynamicAry, depth,
                                                  imgHeight, rsExpFactor);
    }

    [[nodiscard]] Eigen::Vector2d MidPoint() const {
        return {xDynamicAry.at(MID), yDynamicAry.at(MID)};
    }

    [[nodiscard]] double MidPointTime(double readout) const {
        return timeAry[MID] + (rDivHAry[MID] - rsExpFactor) * readout;
    }

    [[nodiscard]] Eigen::Vector2d MidPointVel(double readout) const {
        std::array<double, 3> newTimeAry{};
        for (int i = 0; i < 3; ++i) {
            newTimeAry[i] = timeAry[i] + (rDivHAry[i] - rsExpFactor) * readout;
        }
        return {ns_ikalibr::LagrangePolynomialTripleMidFOD(newTimeAry, xDynamicAry),
                ns_ikalibr::LagrangePolynomialTripleMidFOD(newTimeAry, yDynamicAry)};
    }
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_RGBD_VELOCITY_FACTOR_HPP
