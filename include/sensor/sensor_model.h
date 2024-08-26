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

#ifndef IKALIBR_SENSOR_MODEL_H
#define IKALIBR_SENSOR_MODEL_H

#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct CameraModel {
    enum class CameraModelType : std::uint32_t {
        /**
         * @brief options
         */
        NONE = 1 << 0,

        GS = 1 << 1,
        RS = 1 << 2,

        FIRST_EXPOSURE = 1 << 3,
        MID_EXPOSURE = 1 << 4,
        LAST_EXPOSURE = 1 << 5,

        SENSOR_IMAGE = 1 << 6,
        SENSOR_IMAGE_COMP = 1 << 7,

        SENSOR_IMAGE_GS = SENSOR_IMAGE | GS,
        SENSOR_IMAGE_COMP_GS = SENSOR_IMAGE_COMP | GS,

        SENSOR_IMAGE_RS_FIRST = SENSOR_IMAGE | RS | FIRST_EXPOSURE,
        SENSOR_IMAGE_RS_MID = SENSOR_IMAGE | RS | MID_EXPOSURE,
        SENSOR_IMAGE_RS_LAST = SENSOR_IMAGE | RS | LAST_EXPOSURE,

        SENSOR_IMAGE_COMP_RS_FIRST = SENSOR_IMAGE_COMP_GS | RS | FIRST_EXPOSURE,
        SENSOR_IMAGE_COMP_RS_MID = SENSOR_IMAGE_COMP_GS | RS | MID_EXPOSURE,
        SENSOR_IMAGE_COMP_RS_LAST = SENSOR_IMAGE_COMP_GS | RS | LAST_EXPOSURE,
    };

    static std::string UnsupportedCameraModelMsg(const std::string &modelStr);

    static double RSCameraExposureFactor(const CameraModelType &type);
};

using CameraModelType = CameraModel::CameraModelType;

struct IMUModel {
    enum class IMUModelType : std::uint32_t {
        /**
         * @brief options
         */
        NONE = 1 << 0,
        SBG_IMU = 1 << 1,
        STD_SENSOR_IMU = 1 << 2,
        G = 1 << 3,
        G_NEG = 1 << 4,
        RAD = 1 << 5,
        DEG = 1 << 6,

        SENSOR_IMU = STD_SENSOR_IMU | RAD,      // rad/s, m/s^2
        SENSOR_IMU_G = SENSOR_IMU | G,          // rad/s, g
        SENSOR_IMU_G_NEG = SENSOR_IMU | G_NEG,  // rad/s, -g

        SENSOR_IMU_DEG = STD_SENSOR_IMU | DEG,          // deg/s, m/s^2
        SENSOR_IMU_DEG_G = SENSOR_IMU_DEG | G,          // deg/s, g
        SENSOR_IMU_DEG_G_NEG = SENSOR_IMU_DEG | G_NEG,  // deg/s, -g
    };

    static std::string UnsupportedIMUModelMsg(const std::string &modelStr);
};

using IMUModelType = IMUModel::IMUModelType;

struct LidarModel {
    enum class LidarModelType {
        VLP_16_PACKET,
        VLP_POINTS,

        OUSTER_POINTS,

        PANDAR_XT_POINTS,

        LIVOX_CUSTOM,
    };

    static std::string UnsupportedLiDARModelMsg(const std::string &modelStr);
};

using LidarModelType = LidarModel::LidarModelType;

struct RadarModel {
    enum class RadarModelType {
        AINSTEIN_RADAR,
        AWR1843BOOST_RAW,
        AWR1843BOOST_CUSTOM,
        POINTCLOUD2_POSV,
        POINTCLOUD2_POSIV,
        POINTCLOUD2_XRIO,
    };

    static std::string UnsupportedRadarModelMsg(const std::string &modelStr);
};

using RadarModelType = RadarModel::RadarModelType;

}  // namespace ns_ikalibr

#endif  // IKALIBR_SENSOR_MODEL_H
