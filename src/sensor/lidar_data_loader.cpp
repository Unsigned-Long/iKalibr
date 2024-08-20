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

//
// Created by csl on 10/4/22.
//

#include "sensor/lidar_data_loader.h"
#include "ikalibr/LivoxCustomMsg.h"
#include "sensor_msgs/PointCloud2.h"
#include "spdlog/spdlog.h"
#include "util/status.hpp"
#include "velodyne_msgs/VelodynePacket.h"
#include "velodyne_pointcloud/pointcloudXYZIRT.h"
#include "velodyne_pointcloud/rawdata.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

LiDARDataLoader::Ptr LiDARDataLoader::GetLoader(const std::string &lidarModelStr) {
    // try extract lidar model
    LidarModelType lidarModel;
    try {
        lidarModel = EnumCast::stringToEnum<LidarModelType>(lidarModelStr);
    } catch (...) {
        throw Status(Status::WARNING, LidarModel::UnsupportedLiDARModelMsg(lidarModelStr));
    }
    LiDARDataLoader::Ptr dataLoader;
    switch (lidarModel) {
        case LidarModelType::VLP_16_PACKET:
            dataLoader = Velodyne16::Create(lidarModel);
            break;
        case LidarModelType::VLP_POINTS:
            dataLoader = VelodynePoints::Create(lidarModel);
            break;
        case LidarModelType::OUSTER_POINTS:
            dataLoader = OusterLiDAR::Create(lidarModel);
            break;
        case LidarModelType::PANDAR_XT_POINTS:
            dataLoader = PandarXTLiDAR::Create(lidarModel);
            break;
        case LidarModelType::LIVOX_CUSTOM:
            dataLoader = LivoxLiDAR::Create(lidarModel);
            break;
        default:
            throw Status(Status::WARNING, LidarModel::UnsupportedLiDARModelMsg(lidarModelStr));
    }
    return dataLoader;
}

LidarModelType LiDARDataLoader::GetLiDARModel() const { return _lidarModel; }

// ----------
// Velodyne16
// ----------
Velodyne16::Velodyne16(LidarModelType lidarModel)
    : LiDARDataLoader(lidarModel) {
    SetParameters();
}

Velodyne16::Ptr Velodyne16::Create(LidarModelType lidarModel) {
    return std::make_shared<Velodyne16>(lidarModel);
}

LiDARFrame::Ptr Velodyne16::UnpackScan(const rosbag::MessageInstance &msgInstance) {
    if (_lidarModel == LidarModelType::VLP_16_PACKET) {
        velodyne_msgs::VelodyneScan::ConstPtr scanMsg =
            msgInstance.instantiate<velodyne_msgs::VelodyneScan>();
        CheckMessage<velodyne_msgs::VelodyneScan>(scanMsg);
        return UnpackScan(scanMsg);
    } else {
        throw Status(Status::CRITICAL,
                     "'Velodyne16' data loader only supports the lidar type: 'VLP_16_SIMU' and "
                     "'VLP_16_PACKET'");
    }
}

LiDARFrame::Ptr Velodyne16::UnpackScan(
    const velodyne_msgs::VelodyneScan::ConstPtr &lidarMsg) const {
    LiDARFrame::Ptr output = LiDARFrame::Create(lidarMsg->header.stamp.toSec());

    // point cloud
    auto scan = output->GetScan();
    scan->height = 16;
    scan->width = 24 * (int)lidarMsg->packets.size();
    scan->is_dense = false;
    scan->resize(scan->height * scan->width);
    for (auto &p : scan->points) {
        SET_POS_NAN(p)
    }

    int blockCounter = 0;

    double scanTimestamp = lidarMsg->header.stamp.toSec();

    for (const auto &packet : lidarMsg->packets) {
        float azimuth;
        float azimuthDiff;
        float lastAzimuthDiff = 0;
        float azimuthCorrectedF;
        int azimuthCorrected;
        float x, y, z;

        const auto *raw = (const RAW_PACKET_T *)&packet.data[0];

        for (int block = 0; block < BLOCKS_PER_PACKET; block++, blockCounter++) {
            // Calculate difference between current and next block's azimuth angle.
            azimuth = (float)(raw->blocks[block].rotation);

            if (block < (BLOCKS_PER_PACKET - 1)) {
                azimuthDiff = (float)((ROTATION_MAX_UNITS + raw->blocks[block + 1].rotation -
                                       raw->blocks[block].rotation) %
                                      ROTATION_MAX_UNITS);
                lastAzimuthDiff = azimuthDiff;
            } else {
                azimuthDiff = lastAzimuthDiff;
            }

            for (int firing = 0, k = 0; firing < FIRINGS_PER_BLOCK; firing++) {
                for (int dsr = 0; dsr < SCANS_PER_FIRING; dsr++, k += RAW_SCAN_SIZE) {
                    /** Position Calculation */
                    union TwoBytes tmp {};
                    tmp.bytes[0] = raw->blocks[block].data[k];
                    tmp.bytes[1] = raw->blocks[block].data[k + 1];

                    /** correct for the laser rotation as a function of timing during
                     * the firings **/
                    azimuthCorrectedF =
                        azimuth + (azimuthDiff *
                                   (((float)dsr * DSR_TOFFSET) + ((float)firing * FIRING_TOFFSET)) /
                                   BLOCK_TDURATION);

                    azimuthCorrected = ((int)round(azimuthCorrectedF)) % ROTATION_MAX_UNITS;

                    /*condition added to avoid calculating points which are not
                  in the interesting defined area (minAngle < area < maxAngle)*/
                    if ((azimuthCorrected >= CONFIG.minAngle &&
                         azimuthCorrected <= CONFIG.maxAngle &&
                         CONFIG.minAngle < CONFIG.maxAngle) ||
                        (CONFIG.minAngle > CONFIG.maxAngle &&
                         (azimuthCorrected <= CONFIG.maxAngle ||
                          azimuthCorrected >= CONFIG.minAngle))) {
                        // convert polar coordinates to Euclidean XYZ
                        float distance = (float)tmp.uint * DISTANCE_RESOLUTION;

                        float cos_vert_angle = COS_VERT_ANGLE[dsr];
                        float sin_vert_angle = SIN_VERT_ANGLE[dsr];

                        float cos_rot_angle = COS_ROT_TABLE[azimuthCorrected];
                        float sin_rot_angle = SIN_ROT_TABLE[azimuthCorrected];

                        x = distance * cos_vert_angle * sin_rot_angle;
                        y = distance * cos_vert_angle * cos_rot_angle;
                        z = distance * sin_vert_angle;

                        /** Use standard ROS coordinate system (right-hand rule) */
                        float x_coord = y;
                        float y_coord = -x;
                        float z_coord = z;

                        // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT'
                        // here float intensity = raw->blocks[block].data[k + 2];
                        double point_timestamp =
                            scanTimestamp + GetExactTime(dsr, 2 * blockCounter + firing);

                        IKalibrPoint point_xyz;
                        point_xyz.timestamp = point_timestamp;

                        if (PointInRange(distance)) {
                            point_xyz.x = x_coord;
                            point_xyz.y = y_coord;
                            point_xyz.z = z_coord;
                            // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT'
                            // here point_xyz.intensity = intensity;
                        } else {
                            SET_POS_NAN(point_xyz)
                        }

                        scan->at(2 * blockCounter + firing, SCAN_MAPPING_16[dsr]) = point_xyz;
                    }
                }
            }
        }
    }

    return output;
}

double Velodyne16::GetExactTime(int dsr, int firing) const { return VLP16_TIME_BLOCK[firing][dsr]; }

void Velodyne16::SetParameters() {
    CONFIG.maxRange = 150;
    CONFIG.minRange = 1.0;
    CONFIG.minAngle = 0;
    CONFIG.maxAngle = 36000;
    // Set up cached values for sin and cos of all the possible headings
    for (uint16_t rotIndex = 0; rotIndex < ROTATION_MAX_UNITS; ++rotIndex) {
        auto rotation = (float)(ROTATION_RESOLUTION * (float)rotIndex * M_PI / 180.0);
        COS_ROT_TABLE[rotIndex] = cosf(rotation);
        SIN_ROT_TABLE[rotIndex] = sinf(rotation);
    }

    FIRINGS_PER_BLOCK = 2;
    SCANS_PER_FIRING = 16;
    BLOCK_TDURATION = 110.592f;  // [µs]
    DSR_TOFFSET = 2.304f;        // [µs]
    FIRING_TOFFSET = 55.296f;    // [µs]
    PACKET_TIME = (BLOCKS_PER_PACKET * 2 * FIRING_TOFFSET);

    float vertCorrection[16] = {
        -0.2617993877991494,  0.017453292519943295, -0.22689280275926285,  0.05235987755982989,
        -0.19198621771937624, 0.08726646259971647,  -0.15707963267948966,  0.12217304763960307,
        -0.12217304763960307, 0.15707963267948966,  -0.08726646259971647,  0.19198621771937624,
        -0.05235987755982989, 0.22689280275926285,  -0.017453292519943295, 0.2617993877991494};
    for (int i = 0; i < 16; i++) {
        COS_VERT_ANGLE[i] = std::cos(vertCorrection[i]);
        SIN_VERT_ANGLE[i] = std::sin(vertCorrection[i]);
    }

    SCAN_MAPPING_16[15] = 0;
    SCAN_MAPPING_16[13] = 1;
    SCAN_MAPPING_16[11] = 2;
    SCAN_MAPPING_16[9] = 3;
    SCAN_MAPPING_16[7] = 4;
    SCAN_MAPPING_16[5] = 5;
    SCAN_MAPPING_16[3] = 6;
    SCAN_MAPPING_16[1] = 7;

    SCAN_MAPPING_16[14] = 8;
    SCAN_MAPPING_16[12] = 9;
    SCAN_MAPPING_16[10] = 10;
    SCAN_MAPPING_16[8] = 11;
    SCAN_MAPPING_16[6] = 12;
    SCAN_MAPPING_16[4] = 13;
    SCAN_MAPPING_16[2] = 14;
    SCAN_MAPPING_16[0] = 15;

    for (unsigned int w = 0; w < 1824; w++) {
        for (unsigned int h = 0; h < 16; h++) {
            VLP16_TIME_BLOCK[w][h] = h * 2.304 * 1e-6 + w * 55.296 * 1e-6;  //  16*1824
        }
    }
}

bool Velodyne16::PointInRange(float range) const {
    return (range >= CONFIG.minRange && range <= CONFIG.maxRange);
}

// --------------
// VelodynePoints
// --------------
VelodynePoints::VelodynePoints(LidarModelType lidarModel)
    : LiDARDataLoader(lidarModel) {}

VelodynePoints::Ptr VelodynePoints::Create(LidarModelType lidarModel) {
    return std::make_shared<VelodynePoints>(lidarModel);
}

LiDARFrame::Ptr VelodynePoints::UnpackScan(const rosbag::MessageInstance &msgInstance) {
    sensor_msgs::PointCloud2::ConstPtr lidarMsg =
        msgInstance.instantiate<sensor_msgs::PointCloud2>();

    CheckMessage<sensor_msgs::PointCloud2>(lidarMsg);

    PosIRTPointCloud pcIn;
    pcl::fromROSMsg(*lidarMsg, pcIn);

    double timebase = lidarMsg->header.stamp.toSec();

    IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud());
    cloud->is_dense = false;
    cloud->resize(pcIn.size());

    std::size_t j = 0;
    for (const auto &src : pcIn) {
        if (IS_POS_NAN(src)) {
            continue;
        }

        double depth = sqrt(src.x * src.x + src.y * src.y + src.z * src.z);

        if (depth > 200 || depth < 1) {
            continue;
        } else {
            IKalibrPoint dstPoint;
            dstPoint.x = src.x;
            dstPoint.y = src.y;
            dstPoint.z = src.z;
            dstPoint.timestamp = timebase + src.time;
            // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
            // dstPoint.intensity = src.intensity;
            cloud->at(j++) = dstPoint;
        }
    }
    cloud->resize(j);

    return LiDARFrame::Create(timebase, cloud);
}

// -----------
// OusterLiDAR
// -----------
OusterLiDAR::OusterLiDAR(LidarModelType lidarModel)
    : LiDARDataLoader(lidarModel) {}

OusterLiDAR::Ptr OusterLiDAR::Create(LidarModelType lidarModel) {
    return std::make_shared<OusterLiDAR>(lidarModel);
}

LiDARFrame::Ptr OusterLiDAR::UnpackScan(const rosbag::MessageInstance &msgInstance) {
    sensor_msgs::PointCloud2::ConstPtr lidarMsg =
        msgInstance.instantiate<sensor_msgs::PointCloud2>();

    CheckMessage<sensor_msgs::PointCloud2>(lidarMsg);

    OusterPointCloud pcIn;
    pcl::fromROSMsg(*lidarMsg, pcIn);

    double timebase = lidarMsg->header.stamp.toSec();

    IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud());
    cloud->is_dense = false;
    cloud->resize(pcIn.size());

    std::size_t j = 0;
    for (const auto &src : pcIn) {
        if (IS_POS_NAN(src)) {
            continue;
        }

        double depth = sqrt(src.x * src.x + src.y * src.y + src.z * src.z);

        if (depth > 60 || depth < 1) {
            continue;
        } else {
            IKalibrPoint dstPoint;
            dstPoint.x = src.x;
            dstPoint.y = src.y;
            dstPoint.z = src.z;
            dstPoint.timestamp = timebase + static_cast<double>(src.t) * 1E-9;
            // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
            // dstPoint.intensity = src.intensity;
            cloud->at(j++) = dstPoint;
        }
    }
    cloud->resize(j);

    return LiDARFrame::Create(timebase, cloud);
}

// -------------
// PandarXTLiDAR
// -------------
PandarXTLiDAR::PandarXTLiDAR(LidarModelType lidarModel)
    : LiDARDataLoader(lidarModel) {}

PandarXTLiDAR::Ptr PandarXTLiDAR::Create(LidarModelType lidarModel) {
    return std::make_shared<PandarXTLiDAR>(lidarModel);
}

LiDARFrame::Ptr PandarXTLiDAR::UnpackScan(const rosbag::MessageInstance &msgInstance) {
    sensor_msgs::PointCloud2::ConstPtr lidarMsg =
        msgInstance.instantiate<sensor_msgs::PointCloud2>();

    CheckMessage<sensor_msgs::PointCloud2>(lidarMsg);

    PandarPointCloud pcIn;
    pcl::fromROSMsg(*lidarMsg, pcIn);

    double timebase = lidarMsg->header.stamp.toSec();

    /// point cloud
    IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud);
    cloud->is_dense = false;
    cloud->resize(pcIn.size());

    std::size_t j = 0;
    for (const auto &src : pcIn) {
        if (IS_POS_NAN(src)) {
            continue;
        }

        double depth = sqrt(src.x * src.x + src.y * src.y + src.z * src.z);

        if (depth > 100 || depth < 1) {
            continue;
        } else {
            IKalibrPoint dstPoint;
            dstPoint.x = src.x;
            dstPoint.y = src.y;
            dstPoint.z = src.z;
            dstPoint.timestamp = src.timestamp;
            // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
            // dstPoint.intensity = src.intensity;
            cloud->at(j++) = dstPoint;
        }
    }
    cloud->resize(j);

    return LiDARFrame::Create(timebase, cloud);
}

// ----------
// LivoxLiDAR
// ----------
LivoxLiDAR::LivoxLiDAR(LidarModelType lidarModel)
    : LiDARDataLoader(lidarModel) {}

LivoxLiDAR::Ptr LivoxLiDAR::Create(LidarModelType lidarModel) {
    return std::make_shared<LivoxLiDAR>(lidarModel);
}

LiDARFrame::Ptr LivoxLiDAR::UnpackScan(const rosbag::MessageInstance &msgInstance) {
    ikalibr::LivoxCustomMsg::ConstPtr lidarMsg = msgInstance.instantiate<ikalibr::LivoxCustomMsg>();

    CheckMessage<ikalibr::LivoxCustomMsg>(lidarMsg);

    // from nanosecond to second
    double timebase = lidarMsg->header.stamp.toSec();

    IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud);
    cloud->resize(lidarMsg->point_num);
    std::size_t j = 0;

    for (const auto &src : lidarMsg->points) {
        if (IS_POS_NAN(src)) {
            continue;
        }

        double depth = std::sqrt(src.x * src.x + src.y * src.y + src.z * src.z);

        if (depth > 100 || depth < 1) {
            continue;
        } else {
            IKalibrPoint dstPoint;
            dstPoint.x = src.x;
            dstPoint.y = src.y;
            dstPoint.z = src.z;
            dstPoint.timestamp = timebase + static_cast<double>(src.offset_time) * 1E-9;
            // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
            // dstPoint.intensity = src.reflectivity;
            cloud->at(j++) = dstPoint;
        }
    }
    cloud->resize(j);

    return LiDARFrame::Create(timebase, cloud);
}

}  // namespace ns_ikalibr