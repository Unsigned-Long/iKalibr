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

#include "sensor/imu_data_loader.h"
#include "util/enum_cast.hpp"
#include "util/status.hpp"
#include "spdlog/fmt/fmt.h"
#include "config/configor.h"
#include "ikalibr/SbgImuData.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

IMUDataLoader::IMUDataLoader(IMUModelType imuModel)
    : _imuModel(imuModel) {}

IMUDataLoader::Ptr IMUDataLoader::GetLoader(const std::string &imuModelStr) {
    // try extract radar model
    IMUModelType imuModel;
    try {
        imuModel = EnumCast::stringToEnum<IMUModelType>(imuModelStr);
    } catch (...) {
        throw Status(Status::WARNING, IMUModel::UnsupportedIMUModelMsg(imuModelStr));
    }
    IMUDataLoader::Ptr dataLoader;

    switch (imuModel) {
        case IMUModelType::SBG_IMU:
            dataLoader = SbgIMULoader::Create(imuModel);
            break;
        case IMUModelType::SENSOR_IMU:
            dataLoader = SensorIMULoader::Create(imuModel, 1.0, 1.0);
            break;
        case IMUModelType::SENSOR_IMU_G:
            dataLoader = SensorIMULoader::Create(imuModel, 1.0, Configor::Prior::GravityNorm);
            break;
        case IMUModelType::SENSOR_IMU_G_NEG:
            dataLoader = SensorIMULoader::Create(imuModel, 1.0, -Configor::Prior::GravityNorm);
            break;
        case IMUModelType::SENSOR_IMU_DEG:
            dataLoader = SensorIMULoader::Create(imuModel, DEG_TO_RAD, 1.0);
            break;
        case IMUModelType::SENSOR_IMU_DEG_G:
            dataLoader =
                SensorIMULoader::Create(imuModel, DEG_TO_RAD, Configor::Prior::GravityNorm);
            break;
        case IMUModelType::SENSOR_IMU_DEG_G_NEG:
            dataLoader =
                SensorIMULoader::Create(imuModel, DEG_TO_RAD, -Configor::Prior::GravityNorm);
            break;
        default:
            throw Status(Status::WARNING, IMUModel::UnsupportedIMUModelMsg(imuModelStr));
    }
    return dataLoader;
}

IMUModelType IMUDataLoader::GetIMUModel() const { return _imuModel; }

// ---------------
// SensorIMULoader
// ---------------
SensorIMULoader::SensorIMULoader(IMUModelType imuModel, double g2StdUnit, double a2StdUnit)
    : IMUDataLoader(imuModel),
      g2StdUnit(g2StdUnit),
      a2StdUnit(a2StdUnit) {}

SensorIMULoader::Ptr SensorIMULoader::Create(IMUModelType imuModel,
                                             double g2StdUnit,
                                             double a2StdUnit) {
    return std::make_shared<SensorIMULoader>(imuModel, g2StdUnit, a2StdUnit);
}

IMUFrame::Ptr SensorIMULoader::UnpackFrame(const rosbag::MessageInstance &msgInstance) {
    // imu data item
    sensor_msgs::ImuConstPtr msg = msgInstance.instantiate<sensor_msgs::Imu>();

    CheckMessage<sensor_msgs::Imu>(msg);

    Eigen::Vector3d acce =
        a2StdUnit * Eigen::Vector3d(msg->linear_acceleration.x, msg->linear_acceleration.y,
                                    msg->linear_acceleration.z);
    Eigen::Vector3d gyro =
        g2StdUnit *
        Eigen::Vector3d(msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z);

    return IMUFrame::Create(msg->header.stamp.toSec(), gyro, acce);
}

// ------------
// SbgIMULoader
// ------------
SbgIMULoader::SbgIMULoader(IMUModelType imuModel)
    : IMUDataLoader(imuModel) {}

SbgIMULoader::Ptr SbgIMULoader::Create(IMUModelType imuModel) {
    return std::make_shared<SbgIMULoader>(imuModel);
}

IMUFrame::Ptr SbgIMULoader::UnpackFrame(const rosbag::MessageInstance &msgInstance) {
    // imu data item
    ikalibr::SbgImuData::ConstPtr msg = msgInstance.instantiate<ikalibr::SbgImuData>();

    CheckMessage<ikalibr::SbgImuData>(msg);

    Eigen::Vector3d acce = Eigen::Vector3d(msg->accel.x, msg->accel.y, msg->accel.z);
    Eigen::Vector3d gyro = Eigen::Vector3d(msg->gyro.x, msg->gyro.y, msg->gyro.z);

    return IMUFrame::Create(msg->header.stamp.toSec(), gyro, acce);
}
}  // namespace ns_ikalibr