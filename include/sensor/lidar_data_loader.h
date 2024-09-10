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

#ifndef IKALIBR_LIDAR_DATA_LOADER_H
#define IKALIBR_LIDAR_DATA_LOADER_H

#include "sensor/lidar.h"
#include "rosbag/message_instance.h"
#include "util/enum_cast.hpp"
#include "velodyne_msgs/VelodyneScan.h"
#include "sensor/sensor_model.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

class LiDARDataLoader {
public:
    using Ptr = std::shared_ptr<LiDARDataLoader>;

protected:
    LidarModelType _lidarModel;

public:
    explicit LiDARDataLoader(LidarModelType lidarModel)
        : _lidarModel(lidarModel) {}

    virtual LiDARFrame::Ptr UnpackScan(const rosbag::MessageInstance &msgInstance) = 0;

    static LiDARDataLoader::Ptr GetLoader(const std::string &lidarModelStr);

    [[nodiscard]] LidarModelType GetLiDARModel() const;

    virtual ~LiDARDataLoader() = default;
protected:
    template <class MsgType>
    void CheckMessage(typename MsgType::ConstPtr msg) {
        if (msg == nullptr) {
            throw std::runtime_error(
                "Wrong sensor model: '" + std::string(EnumCast::enumToString(GetLiDARModel())) +
                "' for LiDARs! It's incompatible with the type of ros message to load in!");
        }
    }
};

class Velodyne16 : public LiDARDataLoader {
public:
    using Ptr = std::shared_ptr<Velodyne16>;

public:
    explicit Velodyne16(LidarModelType lidarModel);

    static Velodyne16::Ptr Create(LidarModelType lidarModel);

    // for Velodyne16 ['VLP_16_SIMU' and 'VLP_16_PACKET']
    LiDARFrame::Ptr UnpackScan(const rosbag::MessageInstance &msgInstance) override;

protected:
    [[nodiscard]] LiDARFrame::Ptr UnpackScan(
        const velodyne_msgs::VelodyneScan::ConstPtr &lidarMsg) const;

    [[nodiscard]] double GetExactTime(int dsr, int firing) const;

    void SetParameters();

    [[nodiscard]] bool PointInRange(float range) const;

private:
    static const int RAW_SCAN_SIZE = 3;
    static const int SCANS_PER_BLOCK = 32;
    static const int BLOCK_DATA_SIZE = (SCANS_PER_BLOCK * RAW_SCAN_SIZE);
    constexpr static const float ROTATION_RESOLUTION = 0.01f;
    static const uint16_t ROTATION_MAX_UNITS = 36000u;
    constexpr static const float DISTANCE_RESOLUTION = 0.002f;

    static const uint16_t UPPER_BANK = 0xeeff;
    static const uint16_t LOWER_BANK = 0xddff;

    static const int BLOCKS_PER_PACKET = 12;
    static const int PACKET_STATUS_SIZE = 2;

    int FIRINGS_PER_BLOCK{};
    int SCANS_PER_FIRING{};
    float BLOCK_TDURATION{};
    float DSR_TOFFSET{};
    float FIRING_TOFFSET{};
    float PACKET_TIME{};

    float SIN_ROT_TABLE[ROTATION_MAX_UNITS]{};
    float COS_ROT_TABLE[ROTATION_MAX_UNITS]{};
    float COS_VERT_ANGLE[16]{};
    float SIN_VERT_ANGLE[16]{};
    int SCAN_MAPPING_16[16]{};

    typedef struct RawBlock {
        uint16_t header;    ///< UPPER_BANK or LOWER_BANK
        uint16_t rotation;  ///< 0-35999, divide by 100 to get degrees
        uint8_t data[BLOCK_DATA_SIZE];
    } RAW_BLOCK_T;

    union TwoBytes {
        uint16_t uint;
        uint8_t bytes[2];
    };

    union FourBytes {
        uint32_t uint32;
        float_t float32;
    };

    typedef struct RawPacket {
        RAW_BLOCK_T blocks[BLOCKS_PER_PACKET];
        uint32_t revolution;
        uint8_t status[PACKET_STATUS_SIZE];
    } RAW_PACKET_T;

    /** configuration parameters */
    typedef struct {
        double maxRange;  // maximum range to publish
        double minRange;
        int minAngle;  // minimum angle to publish
        int maxAngle;  // maximum angle to publish
    } Config;
    Config CONFIG{};

    double VLP16_TIME_BLOCK[1824][16]{};
};

class VelodynePoints : public LiDARDataLoader {
public:
    using Ptr = std::shared_ptr<VelodynePoints>;

    explicit VelodynePoints(LidarModelType lidarModel);

    static VelodynePoints::Ptr Create(LidarModelType lidarModel);

    LiDARFrame::Ptr UnpackScan(const rosbag::MessageInstance &msgInstance) override;
};

class OusterLiDAR : public LiDARDataLoader {
public:
    using Ptr = std::shared_ptr<OusterLiDAR>;

public:
    explicit OusterLiDAR(LidarModelType lidarModel);

    static OusterLiDAR::Ptr Create(LidarModelType lidarModel);

    LiDARFrame::Ptr UnpackScan(const rosbag::MessageInstance &msgInstance) override;
};

class PandarXTLiDAR : public LiDARDataLoader {
public:
    using Ptr = std::shared_ptr<PandarXTLiDAR>;

public:
    explicit PandarXTLiDAR(LidarModelType lidarModel);

    static PandarXTLiDAR::Ptr Create(LidarModelType lidarModel);

    LiDARFrame::Ptr UnpackScan(const rosbag::MessageInstance &msgInstance) override;
};

class LivoxLiDAR : public LiDARDataLoader {
public:
    using Ptr = std::shared_ptr<LivoxLiDAR>;

public:
    explicit LivoxLiDAR(LidarModelType lidarModel);

    static LivoxLiDAR::Ptr Create(LidarModelType lidarModel);

    LiDARFrame::Ptr UnpackScan(const rosbag::MessageInstance &msgInstance) override;
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_LIDAR_DATA_LOADER_H
