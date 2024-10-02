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

#include "config/configor.h"
#include "factor/hand_eye_rot_align_factor.hpp"
#include "factor/imu_acce_factor.hpp"
#include "factor/imu_gyro_factor.hpp"
#include "factor/lidar_inertial_align_factor.hpp"
#include "factor/lin_scale_factor.hpp"
#include "factor/point_to_surfel_factor.hpp"
#include "factor/radar_factor.hpp"
#include "factor/radar_inertial_align_factor.hpp"
#include "factor/radar_inertial_rot_align_factor.hpp"
#include "factor/rgbd_inertial_align_factor.hpp"
#include "factor/visual_optical_flow_factor.hpp"
#include "factor/so3_factor.hpp"
#include "factor/visual_inertial_align_factor.hpp"
#include "factor/visual_reproj_factor.hpp"
#include "factor/vel_visual_inertial_align_factor.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
template struct HandEyeRotationAlignFactor<Configor::Prior::SplineOrder>;

template struct IMUAcceFactor<Configor::Prior::SplineOrder, 2>;
template struct IMUAcceFactor<Configor::Prior::SplineOrder, 1>;
template struct IMUAcceFactor<Configor::Prior::SplineOrder, 0>;

template struct IMUGyroFactor<Configor ::Prior::SplineOrder>;

template struct LiDARInertialAlignHelper<Configor::Prior::SplineOrder>;
template struct LiDARInertialAlignFactor<Configor::Prior::SplineOrder>;

template struct LinearScaleDerivFactor<Configor::Prior::SplineOrder, 2>;
template struct LinearScaleDerivFactor<Configor::Prior::SplineOrder, 1>;
template struct LinearScaleDerivFactor<Configor::Prior::SplineOrder, 0>;

template struct PointToSurfelFactor<Configor::Prior::SplineOrder, 2>;
template struct PointToSurfelFactor<Configor::Prior::SplineOrder, 1>;
template struct PointToSurfelFactor<Configor::Prior::SplineOrder, 0>;

template struct RadarFactor<Configor::Prior::SplineOrder, 2>;
template struct RadarFactor<Configor::Prior::SplineOrder, 1>;
template struct RadarFactor<Configor::Prior::SplineOrder, 0>;

template struct RadarInertialAlignHelper<Configor::Prior::SplineOrder>;
template struct RadarInertialAlignFactor<Configor::Prior::SplineOrder>;

template struct RadarInertialRotRoughAlignHelper<Configor::Prior::SplineOrder>;
template struct RadarInertialRotRoughAlignFactor<Configor::Prior::SplineOrder>;

template struct RGBDInertialAlignHelper<Configor::Prior::SplineOrder>;
template struct RGBDInertialAlignFactor<Configor::Prior::SplineOrder>;

template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 2, true>;
template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 2, false>;
template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 1, true>;
template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 1, false>;
template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 0, true>;
template struct VisualOpticalFlowFactor<Configor::Prior::SplineOrder, 0, false>;

template struct SO3Factor<Configor::Prior::SplineOrder>;

template struct VisualInertialAlignHelper<Configor::Prior::SplineOrder>;
template struct VisualInertialAlignFactor<Configor::Prior::SplineOrder>;
template struct VelVisualInertialAlignFactor<Configor::Prior::SplineOrder>;

template struct VisualReProjFactor<Configor::Prior::SplineOrder, 2>;
template struct VisualReProjFactor<Configor::Prior::SplineOrder, 1>;
template struct VisualReProjFactor<Configor::Prior::SplineOrder, 0>;

template struct VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, 2, true>;
template struct VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, 1, true>;
template struct VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, 0, true>;
template struct VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, 2, false>;
template struct VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, 1, false>;
template struct VisualOpticalFlowReProjFactor<Configor::Prior::SplineOrder, 0, false>;

}  // namespace ns_ikalibr
