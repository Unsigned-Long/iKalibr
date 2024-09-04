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

#include "sensor/sensor_model.h"
#include "spdlog/spdlog.h"
#include "util/utils_tpl.hpp"

namespace ns_ikalibr {

std::string CameraModel::UnsupportedCameraModelMsg(const std::string &modelStr) {
    return fmt::format(
        "Unsupported camera Type: '{}'. "
        "Currently supported camera types are: \n"
        "1. SENSOR_IMAGE_GS: "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
        "2. SENSOR_IMAGE_RS_FIRST: first-row exposure, "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
        "3. SENSOR_IMAGE_RS_MID: middle-row exposure, "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
        "4. SENSOR_IMAGE_RS_LAST: last-row exposure, "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Image.html\n"
        "5. SENSOR_IMAGE_COMP_GS: "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/CompressedImage.html\n"
        "6. SENSOR_IMAGE_COMP_RS_FIRST: first-row exposure, "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/CompressedImage.html\n"
        "7. SENSOR_IMAGE_COMP_RS_MID: middle-row exposure, "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/CompressedImage.html\n"
        "8. SENSOR_IMAGE_COMP_RS_LAST: last-row exposure, "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/CompressedImage.html\n"
        "...\n"
        "If you need to use other camera types, "
        "please 'Issues' us on the profile of the github repository.",
        modelStr);
}

double CameraModel::RSCameraExposureFactor(const CameraModelType &type) {
    double exposureFactor = 0.0;
    if (IsOptionWith(CameraModelType::RS, type)) {
        // if it is a RS camera
        if (IsOptionWith(CameraModelType::FIRST_EXPOSURE, type)) {
            exposureFactor = 0.0;
            spdlog::info("RS images are stamped by 'FIRST_EXPOSURE' mode, exposureFactor: '{:.2f}'",
                         exposureFactor);
        } else if (IsOptionWith(CameraModelType::MID_EXPOSURE, type)) {
            exposureFactor = 0.5;
            spdlog::info("RS images are stamped by 'MID_EXPOSURE' mode, exposureFactor: '{:.2f}'",
                         exposureFactor);
        } else if (IsOptionWith(CameraModelType::LAST_EXPOSURE, type)) {
            exposureFactor = 1.0;
            spdlog::info("RS images are stamped by 'LAST_EXPOSURE' mode, exposureFactor: '{:.2f}'",
                         exposureFactor);
        }
    }
    return exposureFactor;
}

std::string IMUModel::UnsupportedIMUModelMsg(const std::string &modelStr) {
    return fmt::format(
        "Unsupported IMU Type: '{}'. "
        "Currently supported IMU types are: \n"
        "1.              SBG_IMU: https://github.com/SBG-Systems/sbg_ros_driver.git\n"
        "2.           SENSOR_IMU: gyro unit (rad/s), acce unit (m/s^2), "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Imu.html\n"
        "3.         SENSOR_IMU_G: gyro unit (rad/s), acce unit (G)\n"
        "4.     SENSOR_IMU_G_NEG: gyro unit (rad/s), acce unit (-G)\n"
        "5.       SENSOR_IMU_DEG: gyro unit (deg/s), acce unit (m/s^2), "
        "https://docs.ros.org/en/noetic/api/sensor_msgs/html/msg/Imu.html\n"
        "6.     SENSOR_IMU_DEG_G: gyro unit (deg/s), acce unit (G)\n"
        "7. SENSOR_IMU_DEG_G_NEG: gyro unit (deg/s), acce unit (-G)\n"
        "...\n"
        "If you need to use other IMU types, "
        "please 'Issues' us on the profile of the github repository.",
        modelStr);
}

std::string LidarModel::UnsupportedLiDARModelMsg(const std::string &modelStr) {
    return fmt::format(
        "Unsupported LiDAR Type: '{}'. "
        "Currently supported LiDAR types are: \n"
        "1.        Velodyne LiDARs: VLP_16_PACKET, VLP_POINTS\n"
        "2.          Ouster LiDARs: OUSTER_POINTS\n"
        "3. Hesai Pandar XT LiDARs: PANDAR_XT_POINTS\n"
        "4.           Livox LiDARs: LIVOX_CUSTOM (the official 'xfer_format'=1, "
        "mid-360 and avia is recommend)\n"
        "...\n"
        "If you need to use other IMU types, "
        "please 'Issues' us on the profile of the github repository.",
        modelStr);
}

std::string RadarModel::UnsupportedRadarModelMsg(const std::string &modelStr) {
    return fmt::format(
        "Unsupported Radar Type: '{}'. "
        "Currently supported radar types are: \n"
        "1.      AINSTEIN_RADAR: https://github.com/AinsteinAI/ainstein_radar.git\n"
        "2.    AWR1843BOOST_RAW: https://github.com/Unsigned-Long/ti_mmwave_rospkg.git\n"
        "3. AWR1843BOOST_CUSTOM: https://github.com/Unsigned-Long/ti_mmwave_rospkg.git\n"
        "4.    POINTCLOUD2_POSV: 'sensor_msgs/PointCloud2' with point format: [x, y, z, "
        "velocity]\n"
        "5.   POINTCLOUD2_POSIV: 'sensor_msgs/PointCloud2' with point format: [x, y, z, "
        "intensity, velocity]\n"
        "6.    POINTCLOUD2_XRIO: 'sensor_msgs/PointCloud2' with x-RIO point format (see "
        "https://github.com/christopherdoer/rio.git)\n"
        "...\n"
        "If you need to use other radar types, "
        "please 'Issues' us on the profile of the github repository.",
        modelStr);
}

}  // namespace ns_ikalibr
