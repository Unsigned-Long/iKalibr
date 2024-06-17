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

#ifndef IKALIBR_CLOUD_DEFINE_HPP
#define IKALIBR_CLOUD_DEFINE_HPP

#include "pcl_ros/point_cloud.h"
#include "pcl/pcl_macros.h"
#include "pcl/point_cloud.h"
#include "pcl/point_types.h"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

// clang-format off

// -------------------------------
// point cloud type for lidar data
// -------------------------------

struct EIGEN_ALIGN16 PointXYZIRT {
    PCL_ADD_POINT4D;              // quad-word XYZ
    float intensity;              // laser intensity reading
    std::uint16_t ring;           // laser ring number
    float time;                   // laser time reading
    PCL_MAKE_ALIGNED_OPERATOR_NEW // ensure proper alignment
};

POINT_CLOUD_REGISTER_POINT_STRUCT(PointXYZIRT,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (float, intensity, intensity)
                                 (std::uint16_t, ring, ring)
                                 (float, time, time))

struct EIGEN_ALIGN16 PointXYZT {
    PCL_ADD_POINT4D;              // quad-word XYZ
    double timestamp;             // laser timestamp
    PCL_MAKE_ALIGNED_OPERATOR_NEW // ensure proper alignment
};

POINT_CLOUD_REGISTER_POINT_STRUCT(PointXYZT,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (double, timestamp, timestamp))

struct EIGEN_ALIGN16 PointXYZIT {
    PCL_ADD_POINT4D;              // quad-word XYZ
    float intensity;              // laser intensity reading
    double timestamp;             // laser timestamp
    PCL_MAKE_ALIGNED_OPERATOR_NEW // ensure proper alignment
};

POINT_CLOUD_REGISTER_POINT_STRUCT(PointXYZIT,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (float, intensity, intensity)
                                 (double, timestamp, timestamp))

// OUSTER lidar
struct EIGEN_ALIGN16 PointXYZIR8Y {
    PCL_ADD_POINT4D;   // quad-word XYZ
    float intensity;   // laser intensity reading
    std::uint8_t ring; // laser ring number
    std::uint32_t t;
    float range;

    PCL_MAKE_ALIGNED_OPERATOR_NEW // ensure proper alignment
};

POINT_CLOUD_REGISTER_POINT_STRUCT(PointXYZIR8Y,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (float, intensity, intensity)
                                 (std::uint8_t, ring, ring)
                                 (std::uint32_t, t, t))

struct EIGEN_ALIGN16 PointXYZITR {
    PCL_ADD_POINT4D;
    float intensity;
    double timestamp;
    uint16_t ring;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

POINT_CLOUD_REGISTER_POINT_STRUCT(PointXYZITR,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (float, intensity, intensity)
                                 (double, timestamp, timestamp)
                                 (std::uint16_t, ring, ring))
// pcl
using PosPoint = pcl::PointXYZ;
using PosPointCloud = pcl::PointCloud<PosPoint>;

using PosIPoint = pcl::PointXYZI;
using PosIPointCloud = pcl::PointCloud<PosIPoint>;

using PosIRTPoint = PointXYZIRT;
using PosIRTPointCloud = pcl::PointCloud<PosIRTPoint>;

using PosTPoint = PointXYZT;
using PosTPointCloud = pcl::PointCloud<PosTPoint>;

using PosITPoint = PointXYZIT;
using PosITPointCloud = pcl::PointCloud<PosITPoint>;

using OusterPoint = PointXYZIR8Y;
using OusterPointCloud = pcl::PointCloud<OusterPoint>;

using ColorPoint = pcl::PointXYZRGBA;
using ColorPointCloud = pcl::PointCloud<ColorPoint>;

using PandarPoint = PointXYZITR;
using PandarPointCloud = pcl::PointCloud<PandarPoint>;

using IKalibrPoint = PointXYZT;
using IKalibrPointCloud = pcl::PointCloud<IKalibrPoint>;

#define IS_POS_NAN(p) (std::isnan((p).x) || std::isnan((p).y) || std::isnan((p).z))

#define SET_POS_NAN(p) {(p).x = (p).y = (p).z = NAN;}

// -------------------------------
// point cloud type for radar data
// -------------------------------

struct EIGEN_ALIGN16 RadarTargetPOSV {
    PCL_ADD_POINT4D;   // quad-word XYZ
    float velocity;    // radial velocity

    PCL_MAKE_ALIGNED_OPERATOR_NEW // ensure proper alignment
};

POINT_CLOUD_REGISTER_POINT_STRUCT(RadarTargetPOSV,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (float, velocity, velocity)
)
using RadarPOSVCloud = pcl::PointCloud<RadarTargetPOSV>;

struct EIGEN_ALIGN16 RadarTargetPOSIV {
    PCL_ADD_POINT4D;   // quad-word XYZ
    float intensity;   // laser intensity reading
    float velocity;    // radial velocity

    PCL_MAKE_ALIGNED_OPERATOR_NEW // ensure proper alignment
};

POINT_CLOUD_REGISTER_POINT_STRUCT(RadarTargetPOSIV,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (float, intensity, intensity)
                                 (float, velocity, velocity)
)
using RadarPOSIVCloud = pcl::PointCloud<RadarTargetPOSIV>;

struct EIGEN_ALIGN16 RadarTargetXRIO {
    PCL_ADD_POINT4D;      // position in [m]
    float snr_db;         // CFAR cell to side noise ratio in [dB]
    float v_doppler_mps;  // Doppler's velocity in [m/s]
    float noise_db;       // CFAR noise level of the side of the detected cell in [dB]
    float range;          // range in [m]
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

POINT_CLOUD_REGISTER_POINT_STRUCT(RadarTargetXRIO,
                                 (float, x, x)
                                 (float, y, y)
                                 (float, z, z)
                                 (float, snr_db, snr_db)
                                 (float, v_doppler_mps, v_doppler_mps)
                                 (float, noise_db, noise_db)
                                 (float, range, range)
)
using RadarXRIOCloud = pcl::PointCloud<RadarTargetXRIO>;

#define SET_POSV_POINT_NAN(p)        \
    {                                \
        (p).velocity = NAN;          \
        (p).x = (p).y = (p).z = NAN; \
    }

#define SET_POSIV_POINT_NAN(p)       \
    {                                \
        (p).velocity = NAN;          \
        (p).x = (p).y = (p).z = NAN; \
        (p).intensity = NAN;         \
    }


#endif //IKALIBR_CLOUD_DEFINE_HPP
