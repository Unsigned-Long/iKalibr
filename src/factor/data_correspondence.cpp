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

#include "factor/data_correspondence.h"
#include "sensor/camera.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
PointToSurfelCorr::PointToSurfelCorr(double timestamp,
                                     Eigen::Vector3d pInScan,
                                     double weight,
                                     Eigen::Vector4d surfelInW)
    : timestamp(timestamp),
      pInScan(std::move(pInScan)),
      weight(weight),
      surfelInW(std::move(surfelInW)) {}

PointToSurfelCorr::Ptr PointToSurfelCorr::Create(double timestamp,
                                                 const Eigen::Vector3d& pInScan,
                                                 double weight,
                                                 const Eigen::Vector4d& surfelInW) {
    return std::make_shared<PointToSurfelCorr>(timestamp, pInScan, weight, surfelInW);
}

VisualReProjCorr::VisualReProjCorr(double ti,
                                   double tj,
                                   Eigen::Vector2d fi,
                                   Eigen::Vector2d fj,
                                   double li,
                                   double lj,
                                   double weight)
    : ti(ti),
      tj(tj),
      fi(std::move(fi)),
      fj(std::move(fj)),
      li(li),
      lj(lj),
      weight(weight) {}

VisualReProjCorr::Ptr VisualReProjCorr::Create(double ti,
                                               double tj,
                                               const Eigen::Vector2d& fi,
                                               const Eigen::Vector2d& fj,
                                               double li,
                                               double lj,
                                               double weight) {
    return std::make_shared<VisualReProjCorr>(ti, tj, fi, fj, li, lj, weight);
}

VisualReProjCorr::VisualReProjCorr() = default;

VisualReProjCorrSeq::VisualReProjCorrSeq() = default;

OpticalFlowCorr::OpticalFlowCorr(const std::array<double, 3>& timeAry,
                                 const std::array<double, 3>& xTraceAry,
                                 const std::array<double, 3>& yTraceAry,
                                 double depth,
                                 const CameraFramePtr& frame,
                                 double rsExpFactor)
    : timeAry(timeAry),
      xTraceAry(xTraceAry),
      yTraceAry(yTraceAry),
      rdFactorAry(),
      depth(depth),
      invDepth(depth > 1E-3 ? 1.0 / depth : -1.0),
      frame(frame),
      withDepthObservability(false) {
    int imgHeight = frame->GetImage().rows;
    for (int i = 0; i < 3; ++i) {
        rdFactorAry[i] = yTraceAry[i] / (double)imgHeight - rsExpFactor;
    }
}

OpticalFlowCorr::Ptr OpticalFlowCorr::Create(const std::array<double, 3>& timeAry,
                                             const std::array<double, 3>& xDynamicAry,
                                             const std::array<double, 3>& yDynamicAry,
                                             double depth,
                                             const CameraFramePtr& frame,
                                             double rsExpFactor) {
    return std::make_shared<OpticalFlowCorr>(timeAry, xDynamicAry, yDynamicAry, depth, frame,
                                             rsExpFactor);
}

Eigen::Vector2d OpticalFlowCorr::FirPoint() const { return {xTraceAry.at(FIR), yTraceAry.at(FIR)}; }

Eigen::Vector2d OpticalFlowCorr::MidPoint() const { return {xTraceAry.at(MID), yTraceAry.at(MID)}; }

Eigen::Vector2d OpticalFlowCorr::LastPoint() const {
    return {xTraceAry.at(LAST), yTraceAry.at(LAST)};
}

double OpticalFlowCorr::MidReadoutFactor() const { return rdFactorAry[MID]; }
}  // namespace ns_ikalibr
