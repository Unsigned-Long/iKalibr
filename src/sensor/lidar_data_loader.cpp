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
#include "util/status.hpp"
#include "spdlog/spdlog.h"
#include "ikalibr/LivoxCustomMsg.h"

_3_

namespace ns_ikalibr {

    LiDARDataLoader::Ptr LiDARDataLoader::GetLoader(const std::string &lidarModelStr) {
        // try extract lidar model
        LidarModelType lidarModel;
        try {
            lidarModel = EnumCast::stringToEnum<LidarModelType>(lidarModelStr);
        } catch (...) {
            throw Status(
                    Status::WARNING,
                    "Unsupported LiDAR Type: '{}'. "
                    "Currently supported LiDAR types are: \n"
                    "(1)        Velodyne LiDARs: VLP_16_PACKET, VLP_16_SIMU, VLP_16_POINTS, VLP_32E_POINTS\n"
                    "(2)          Ouster LiDARs: OUSTER_16_POINTS OUSTER_32_POINTS, OUSTER_64_POINTS, OUSTER_128_POINTS\n"
                    "(3) Hesai Pandar XT LiDARs: PANDAR_XT_16, PANDAR_XT_32\n"
                    "(4)           Livox LiDARs: LIVOX_MID_360, LIVOX_AVIA\n"
                    "...\n", lidarModelStr
            );
        }
        LiDARDataLoader::Ptr lidarDataLoader;
        switch (lidarModel) {
            case LidarModelType::VLP_16_PACKET :
                lidarDataLoader = Velodyne16::Create(lidarModel);
                break;
            case LidarModelType::VLP_16_POINTS:
            case LidarModelType::VLP_32E_POINTS:
                lidarDataLoader = VelodynePoints::Create(lidarModel);
                break;
            case LidarModelType::OUSTER_16_POINTS:
            case LidarModelType::OUSTER_32_POINTS:
            case LidarModelType::OUSTER_64_POINTS:
            case LidarModelType::OUSTER_128_POINTS:
                lidarDataLoader = OusterLiDAR::Create(lidarModel);
                break;
            case LidarModelType::PANDAR_XT_16:
            case LidarModelType::PANDAR_XT_32:
                lidarDataLoader = PandarXTLiDAR::Create(lidarModel);
                break;
            case LidarModelType::LIVOX_MID_360:
            case LidarModelType::LIVOX_AVIA:
                lidarDataLoader = LivoxLiDAR::Create(lidarModel);
                break;
        }
        return lidarDataLoader;
    }

    LidarModelType LiDARDataLoader::GetLiDARModel() const {
        return _lidarModel;
    }

    // ----------
    // Velodyne16
    // ----------
    Velodyne16::Velodyne16(LidarModelType lidarModel) : LiDARDataLoader(lidarModel) { SetParameters(); }

    Velodyne16::Ptr Velodyne16::Create(LidarModelType lidarModel) {
        return std::make_shared<Velodyne16>(lidarModel);
    }

    LiDARFrame::Ptr Velodyne16::UnpackScan(const rosbag::MessageInstance &msgInstance) {
        if (_lidarModel == LidarModelType::VLP_16_PACKET) {
            velodyne_msgs::VelodyneScan::ConstPtr scanMsg = msgInstance.instantiate<velodyne_msgs::VelodyneScan>();
            CheckMessage<velodyne_msgs::VelodyneScan>(scanMsg);
            return UnpackScan(scanMsg);
        } else {
            throw Status(
                    Status::CRITICAL,
                    "'Velodyne16' data loader only supports the lidar type: 'VLP_16_SIMU' and 'VLP_16_PACKET'"
            );
        }
    }

    LiDARFrame::Ptr Velodyne16::UnpackScan(const velodyne_msgs::VelodyneScan::ConstPtr &lidarMsg) const {

        LiDARFrame::Ptr output = LiDARFrame::Create(lidarMsg->header.stamp.toSec());

        // point cloud
        auto scan = output->GetScan();
        scan->height = 16;
        scan->width = 24 * (int) lidarMsg->packets.size();
        scan->is_dense = false;
        scan->resize(scan->height * scan->width);
        for (auto &p: scan->points) {SET_POS_NAN(p)}

        int blockCounter = 0;

        double scanTimestamp = lidarMsg->header.stamp.toSec();

        for (const auto &packet: lidarMsg->packets) {
            float azimuth;
            float azimuthDiff;
            float lastAzimuthDiff = 0;
            float azimuthCorrectedF;
            int azimuthCorrected;
            float x, y, z;

            const auto *raw = (const RAW_PACKET_T *) &packet.data[0];

            for (int block = 0; block < BLOCKS_PER_PACKET; block++, blockCounter++) {
                // Calculate difference between current and next block's azimuth angle.
                azimuth = (float) (raw->blocks[block].rotation);

                if (block < (BLOCKS_PER_PACKET - 1)) {
                    azimuthDiff =
                            (float) ((ROTATION_MAX_UNITS + raw->blocks[block + 1].rotation -
                                      raw->blocks[block].rotation) %
                                     ROTATION_MAX_UNITS);
                    lastAzimuthDiff = azimuthDiff;
                } else {
                    azimuthDiff = lastAzimuthDiff;
                }

                for (int firing = 0, k = 0; firing < FIRINGS_PER_BLOCK; firing++) {
                    for (int dsr = 0; dsr < SCANS_PER_FIRING; dsr++, k += RAW_SCAN_SIZE) {
                        /** Position Calculation */
                        union TwoBytes tmp{};
                        tmp.bytes[0] = raw->blocks[block].data[k];
                        tmp.bytes[1] = raw->blocks[block].data[k + 1];

                        /** correct for the laser rotation as a function of timing during
                         * the firings **/
                        azimuthCorrectedF =
                                azimuth + (azimuthDiff * (((float) dsr * DSR_TOFFSET) +
                                                          ((float) firing * FIRING_TOFFSET)) / BLOCK_TDURATION);

                        azimuthCorrected = ((int) round(azimuthCorrectedF)) % ROTATION_MAX_UNITS;

                        /*condition added to avoid calculating points which are not
                      in the interesting defined area (minAngle < area < maxAngle)*/
                        if ((azimuthCorrected >= CONFIG.minAngle &&
                             azimuthCorrected <= CONFIG.maxAngle &&
                             CONFIG.minAngle < CONFIG.maxAngle) ||
                            (CONFIG.minAngle > CONFIG.maxAngle &&
                             (azimuthCorrected <= CONFIG.maxAngle ||
                              azimuthCorrected >= CONFIG.minAngle))) {
                            // convert polar coordinates to Euclidean XYZ
                            float distance = (float) tmp.uint * DISTANCE_RESOLUTION;

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

                            // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
                            // float intensity = raw->blocks[block].data[k + 2];
                            double point_timestamp =
                                    scanTimestamp + GetExactTime(dsr, 2 * blockCounter + firing);

                            // TODO: using macro 'set post point nan'
                            IKalibrPoint point_xyz;
                            point_xyz.timestamp = point_timestamp;

                            if (PointInRange(distance)) {
                                point_xyz.x = x_coord;
                                point_xyz.y = y_coord;
                                point_xyz.z = z_coord;
                                // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
                                // point_xyz.intensity = intensity;
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

    double Velodyne16::GetExactTime(int dsr, int firing) const {
        return VLP16_TIME_BLOCK[firing][dsr];
    }

    void Velodyne16::SetParameters() {
        CONFIG.maxRange = 150;
        CONFIG.minRange = 1.0;
        CONFIG.minAngle = 0;
        CONFIG.maxAngle = 36000;
        // Set up cached values for sin and cos of all the possible headings
        for (uint16_t rotIndex = 0; rotIndex < ROTATION_MAX_UNITS; ++rotIndex) {
            auto rotation = (float) angles::from_degrees(ROTATION_RESOLUTION * (float) rotIndex);
            COS_ROT_TABLE[rotIndex] = cosf(rotation);
            SIN_ROT_TABLE[rotIndex] = sinf(rotation);
        }

        FIRINGS_PER_BLOCK = 2;
        SCANS_PER_FIRING = 16;
        BLOCK_TDURATION = 110.592f;// [µs]
        DSR_TOFFSET = 2.304f;      // [µs]
        FIRING_TOFFSET = 55.296f;  // [µs]
        PACKET_TIME = (BLOCKS_PER_PACKET * 2 * FIRING_TOFFSET);

        float vertCorrection[16] = {
                -0.2617993877991494, 0.017453292519943295, -0.22689280275926285,
                0.05235987755982989, -0.19198621771937624, 0.08726646259971647,
                -0.15707963267948966, 0.12217304763960307, -0.12217304763960307,
                0.15707963267948966, -0.08726646259971647, 0.19198621771937624,
                -0.05235987755982989, 0.22689280275926285, -0.017453292519943295,
                0.2617993877991494};
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
                VLP16_TIME_BLOCK[w][h] = h * 2.304 * 1e-6 + w * 55.296 * 1e-6;//  16*1824
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
            : LiDARDataLoader(lidarModel), FIRST_MSG(true),
              HAS_TIME_FIELD(false), HAS_RING_FIELD(false), ONE_SCAN_ANGLE(360.) {
        if (_lidarModel == LidarModelType::VLP_16_POINTS) {
            NUM_LASERS = 16;
            NUM_FIRING = 1824;// 76*24
        } else if (_lidarModel == LidarModelType::VLP_32E_POINTS) {
            NUM_LASERS = 32;
            NUM_FIRING = 2170;// 180*12
        } else {
            throw Status(
                    Status::CRITICAL,
                    "'VelodynePoints' data loader only supports the lidar type: 'VLP_16_POINTS' and 'VLP_32E_POINTS'"
            );
        }
        HORIZON_RESOLUTION = ONE_SCAN_ANGLE / (double) NUM_FIRING;

        if (_lidarModel == LidarModelType::VLP_32E_POINTS) {
            // laserID(dsr channel)
            // channel 0,1,...,31
            laserID_mapping[31] = 0;
            laserID_mapping[29] = 1;
            laserID_mapping[27] = 2;
            laserID_mapping[25] = 3;
            laserID_mapping[23] = 4;
            laserID_mapping[21] = 5;
            laserID_mapping[19] = 6;
            laserID_mapping[17] = 7;
            laserID_mapping[15] = 8;
            laserID_mapping[13] = 9;
            laserID_mapping[11] = 10;
            laserID_mapping[9] = 11;
            laserID_mapping[7] = 12;
            laserID_mapping[5] = 13;
            laserID_mapping[3] = 14;
            laserID_mapping[1] = 15;

            laserID_mapping[30] = 16;
            laserID_mapping[28] = 17;
            laserID_mapping[26] = 18;
            laserID_mapping[24] = 19;
            laserID_mapping[22] = 20;
            laserID_mapping[20] = 21;
            laserID_mapping[18] = 22;
            laserID_mapping[16] = 23;
            laserID_mapping[14] = 24;
            laserID_mapping[12] = 25;
            laserID_mapping[10] = 26;
            laserID_mapping[8] = 27;
            laserID_mapping[6] = 28;
            laserID_mapping[4] = 29;
            laserID_mapping[2] = 30;
            laserID_mapping[0] = 31;
        }
    }

    VelodynePoints::Ptr VelodynePoints::Create(LidarModelType lidarModel) {
        return std::make_shared<VelodynePoints>(lidarModel);
    }

    bool VelodynePoints::CheckCloudField(const sensor_msgs::PointCloud2::ConstPtr &cloud_msg,
                                         const std::string &desc) {
        bool has_the_field = false;
        for (const auto &field: cloud_msg->fields) {
            if (field.name == desc) {
                has_the_field = true;
                break;
            }
        }

        if (!has_the_field) {
            spdlog::warn("PointCloud2 not has channel '{}', please configure your point cloud data!", desc);
        }
        return has_the_field;
    }

    double VelodynePoints::ClockwiseAngle(double befAngle, double afterAngle) {
        // atan2(py, px)
        // 1quadrant-1quadrant = 45 - 30
        // 1quadrant-4quadrant = 45 - (-30) = 75
        // 1quadrant-3quadrant = 45 - (-120) = 165
        // 1quadrant-2quadrant = 45 - 120 = -75 + 360 = 285
        double dAngle = befAngle - afterAngle;
        if (dAngle < 0) dAngle += 360;

        return dAngle;
    }

    double VelodynePoints::RotationTravelledClockwise(double nowAngle, bool resetCnt) {
        static double startAngleDegree = 0;
        static bool halfRotPointer = false;
        static double rotCircleCnt = 0;

        if (resetCnt) {
            startAngleDegree = nowAngle;
            halfRotPointer = false;
            rotCircleCnt = 0;

            return 0;
        } else {
            double dAngle = ClockwiseAngle(startAngleDegree, nowAngle);
            // half a circle
            if (dAngle > 100 && dAngle < 270) {
                halfRotPointer = true;
            }
            // a circle
            if (halfRotPointer && dAngle < 80) {
                halfRotPointer = false;
                rotCircleCnt += 360;
            }

            return rotCircleCnt + dAngle;// rot_travelled
        }
    }

    bool VelodynePoints::InitScanParam(const sensor_msgs::PointCloud2::ConstPtr &lidarMsg) {
        HAS_TIME_FIELD = CheckCloudField(lidarMsg, "time");
        HAS_RING_FIELD = CheckCloudField(lidarMsg, "ring");

        // if (!HAS_RING_FIELD) return false;

        if (!HAS_TIME_FIELD) {
            spdlog::warn("Input PointCloud2 not has channel [time]. Calculate timestamp "
                         "of each point assuming constant rotation speed");
        }

        PosIPointCloud cloud;
        pcl::fromROSMsg(*lidarMsg, cloud);

        bool is_first_horizon_angle = true;

        double rotation_travelled = 0;
        for (auto const &pointIn: cloud) {
            if (IS_POS_NAN(pointIn))
                continue;
            static double rad2deg = 180.0 / M_PI;
            double horizon_angle = atan2(pointIn.y, pointIn.x) * rad2deg;
            if (is_first_horizon_angle) {
                is_first_horizon_angle = false;
                RotationTravelledClockwise(horizon_angle, true);
            }
            rotation_travelled = RotationTravelledClockwise(horizon_angle);
        }

        ONE_SCAN_ANGLE = round(rotation_travelled / 360.) * 363.;

        HORIZON_RESOLUTION = ONE_SCAN_ANGLE / (double) NUM_FIRING;

        return true;
    }

    LiDARFrame::Ptr VelodynePoints::UnpackScan(const rosbag::MessageInstance &msgInstance) {
        sensor_msgs::PointCloud2::ConstPtr lidarMsg = msgInstance.instantiate<sensor_msgs::PointCloud2>();

        CheckMessage<sensor_msgs::PointCloud2>(lidarMsg);

        if (FIRST_MSG) {
            InitScanParam(lidarMsg);
            FIRST_MSG = false;
        }

        static const double LIDAR_FOV_DOWN = -15.0;
        static const double LIDAR_FOV_RESOLUTION = 2.0;

        PosIRTPointCloud pcIn;
        pcl::fromROSMsg(*lidarMsg, pcIn);

        IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud());
        double timebase = lidarMsg->header.stamp.toSec();

        /// point cloud
        cloud->clear();
        cloud->height = NUM_LASERS;
        cloud->width = NUM_FIRING;
        cloud->is_dense = false;
        cloud->resize(cloud->height * cloud->width);

        int num_points = NUM_LASERS * NUM_FIRING;
        for (int k = 0; k < num_points; k++) {SET_POS_NAN(cloud->points[k])}

        int lastFiring = 0;
        int firingOrderPositive = 0, firingOrderPositiveCnt = 0;
        int firingOrderReverse = 0, firingOrderReverseCnt = 0;

        bool isFirstHorizonAngle = true;
        for (int i = 0; i < static_cast<int>(pcIn.size()); ++i) {
            const PosIRTPoint &pointIn = pcIn.points[i];

            if (IS_POS_NAN(pointIn)) { continue; }

            static double RAD_2_DEG = 180.0 / M_PI;

            double horizonAngle = atan2(pointIn.y, pointIn.x) * RAD_2_DEG;
            if (isFirstHorizonAngle) {
                isFirstHorizonAngle = false;
                RotationTravelledClockwise(horizonAngle, true);
            }
            double rotationTravelled = RotationTravelledClockwise(horizonAngle);

            int firing = (int) std::round(rotationTravelled / HORIZON_RESOLUTION);

            if (i > 0) {
                int df = firing - lastFiring;
                if (df >= 0) {
                    firingOrderPositive += df;
                    firingOrderPositiveCnt++;
                } else {
                    firingOrderReverse += df;
                    firingOrderReverseCnt++;
                }
            }
            lastFiring = firing;

            if (firing < 0 || firing >= NUM_FIRING) continue;

            double dt;
            if (HAS_TIME_FIELD) {
                dt = pointIn.time;
            } else {
                dt = 0.1 * rotationTravelled / ONE_SCAN_ANGLE;
            }

            double depth = sqrt(pointIn.x * pointIn.x + pointIn.y * pointIn.y +
                                pointIn.z * pointIn.z);

            int ring;
            if (HAS_RING_FIELD) {
                ring = pointIn.ring;
            } else {
                double pitch = asin(pointIn.z / depth) * RAD_2_DEG;
                ring = (int) std::round((pitch - LIDAR_FOV_DOWN) / LIDAR_FOV_RESOLUTION);
                if (ring < 0 || ring >= 16) continue;
            }

            IKalibrPoint pointOut;
            pointOut.timestamp = timebase + dt;

            if (depth > 200 || depth < 1) {
                SET_POS_NAN(pointOut)
            } else {
                pointOut.x = pointIn.x;
                pointOut.y = pointIn.y;
                pointOut.z = pointIn.z;
                // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
                // pointOut.intensity = pointIn.intensity;
            }

            cloud->at(firing, ring) = pointOut;
        }

        return LiDARFrame::Create(timebase, cloud);
    }

    // -----------
    // OusterLiDAR
    // -----------
    OusterLiDAR::OusterLiDAR(LidarModelType lidarModel) : LiDARDataLoader(lidarModel), NUM_FIRING(2048) {}

    OusterLiDAR::Ptr OusterLiDAR::Create(LidarModelType lidarModel) {
        return std::make_shared<OusterLiDAR>(lidarModel);
    }

    LiDARFrame::Ptr OusterLiDAR::UnpackScan(const rosbag::MessageInstance &msgInstance) {
        sensor_msgs::PointCloud2::ConstPtr lidarMsg = msgInstance.instantiate<sensor_msgs::PointCloud2>();

        CheckMessage<sensor_msgs::PointCloud2>(lidarMsg);

        OusterPointCloud pcIn;
        pcl::fromROSMsg(*lidarMsg, pcIn);

        int ringNumber;
        switch (_lidarModel) {
            case LidarModelType::OUSTER_16_POINTS:
                ringNumber = 16;
                break;
            case LidarModelType::OUSTER_32_POINTS:
                ringNumber = 32;
                break;
            case LidarModelType::OUSTER_64_POINTS:
                ringNumber = 64;
                break;
            case LidarModelType::OUSTER_128_POINTS:
                ringNumber = 128;
                break;
            default:
                throw Status(
                        Status::CRITICAL,
                        "'OusterLiDAR' data loader only supports the lidar type: "
                        "'OUSTER_16_POINTS', 'OUSTER_32_POINTS', 'OUSTER_64_POINTS' and 'OUSTER_128_POINTS'"
                );
        }

        int ring_step = (int) pcIn.height / ringNumber;

        assert(ring_step >= 1 && "OusterRingNo too large");

        IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud);

        double timebase = lidarMsg->header.stamp.toSec();

        /// point cloud
        cloud->clear();
        cloud->height = ringNumber;
        cloud->width = NUM_FIRING;
        cloud->is_dense = false;
        cloud->resize(cloud->height * cloud->width);

        for (int h = 0; h < ringNumber; h++) {
            int hInSelected = h * ring_step;
            for (int w = 0; w < NUM_FIRING; w++) {
                const auto &src = pcIn.at(w, hInSelected);
                if (IS_POS_NAN(src)) {
                    SET_POS_NAN(cloud->at(w, h))
                    continue;
                }

                double depth = sqrt(src.x * src.x + src.y * src.y + src.z * src.z);

                IKalibrPoint dstPoint;
                dstPoint.timestamp = timebase + (float) src.t * 1E-9f;
                if (depth > 60 || depth < 1) {
                    dstPoint.x = src.x;
                    dstPoint.y = src.y;
                    dstPoint.z = src.z;
                    // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
                    // dstPoint.intensity = src.intensity;
                } else {
                    SET_POS_NAN(dstPoint)
                }

                cloud->at(w, h) = dstPoint;
            }
        }
        return LiDARFrame::Create(timebase, cloud);
    }

    // -------------
    // PandarXTLiDAR
    // -------------
    PandarXTLiDAR::PandarXTLiDAR(LidarModelType lidarModel) : LiDARDataLoader(lidarModel) {}

    PandarXTLiDAR::Ptr PandarXTLiDAR::Create(LidarModelType lidarModel) {
        return std::make_shared<PandarXTLiDAR>(lidarModel);
    }

    LiDARFrame::Ptr PandarXTLiDAR::UnpackScan(const rosbag::MessageInstance &msgInstance) {
        sensor_msgs::PointCloud2::ConstPtr lidarMsg = msgInstance.instantiate<sensor_msgs::PointCloud2>();

        CheckMessage<sensor_msgs::PointCloud2>(lidarMsg);

        PandarPointCloud pcIn;
        pcl::fromROSMsg(*lidarMsg, pcIn);

        double timebase = lidarMsg->header.stamp.toSec();

        /// point cloud
        IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud);
        cloud->is_dense = false;
        cloud->resize(pcIn.size());

        std::size_t j = 0;
        for (const auto &src: pcIn) {
            if (IS_POS_NAN(src)) { continue; }

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
    LivoxLiDAR::LivoxLiDAR(LidarModelType lidarModel) : LiDARDataLoader(lidarModel) {}

    LivoxLiDAR::Ptr LivoxLiDAR::Create(LidarModelType lidarModel) {
        return std::make_shared<LivoxLiDAR>(lidarModel);
    }

    LiDARFrame::Ptr LivoxLiDAR::UnpackScan(const rosbag::MessageInstance &msgInstance) {
        ikalibr::LivoxCustomMsg::ConstPtr lidarMsg = msgInstance.instantiate<ikalibr::LivoxCustomMsg>();

        CheckMessage<ikalibr::LivoxCustomMsg>(lidarMsg);

        // from nanosecond to second
        double timebase = static_cast<double>(lidarMsg->timebase) * 1E-9;

        IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud);
        cloud->resize(lidarMsg->point_num);
        std::size_t j = 0;

        for (const auto &src: lidarMsg->points) {
            if (IS_POS_NAN(src)) { continue; }

            double depth = std::sqrt(src.x * src.x + src.y * src.y + src.z * src.z);

            if (depth > 100 || depth < 1) {
                continue;
            } else {
                IKalibrPoint dstPoint;
                dstPoint.x = src.x;
                dstPoint.y = src.y;
                dstPoint.z = src.z;
                dstPoint.timestamp = static_cast<double>(src.offset_time) * 1E-9 + timebase;
                // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
                // dstPoint.intensity = src.reflectivity;
                cloud->at(j++) = dstPoint;
            }
        }
        cloud->resize(j);

        return LiDARFrame::Create(timebase, cloud);
    }


}// namespace ns_ikalibr