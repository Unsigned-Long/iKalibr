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
#include "spdlog/spdlog.h"
#include "util/status.hpp"
#include "magic_enum_flags.hpp"
#include "ros/package.h"
#include "filesystem"
#include "cereal/types/vector.hpp"
#include "cereal/types/set.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
const static std::map<std::string, OutputOption> OutputOptionMap = {
    {"NONE", OutputOption::NONE},
    {"ParamInEachIter", OutputOption::ParamInEachIter},
    {"BSplines", OutputOption::BSplines},
    {"LiDARMaps", OutputOption::LiDARMaps},
    {"VisualMaps", OutputOption::VisualMaps},
    {"RadarMaps", OutputOption::RadarMaps},
    {"HessianMat", OutputOption::HessianMat},
    {"VisualLiDARCovisibility", OutputOption::VisualLiDARCovisibility},
    {"VisualKinematics", OutputOption::VisualKinematics},
    {"ColorizedLiDARMap", OutputOption::ColorizedLiDARMap},
    {"AlignedInertialMes", OutputOption::AlignedInertialMes},
    {"VisualReprojErrors", OutputOption::VisualReprojErrors},
    {"RadarDopplerErrors", OutputOption::RadarDopplerErrors},
    {"ALL", OutputOption::ALL},
};

// ------------------------
// static initialized filed
// ------------------------
Configor::DataStream Configor::dataStream = {};
Configor::Prior Configor::prior = {};
Configor::Preference Configor::preference = {};

std::map<std::string, Configor::DataStream::IMUConfig> Configor::DataStream::IMUTopics = {};
std::map<std::string, Configor::DataStream::RadarConfig> Configor::DataStream::RadarTopics = {};
std::map<std::string, Configor::DataStream::LiDARConfig> Configor::DataStream::LiDARTopics = {};
std::map<std::string, Configor::DataStream::CameraConfig> Configor::DataStream::CameraTopics = {};
std::map<std::string, Configor::DataStream::RGBDConfig> Configor::DataStream::RGBDTopics = {};
std::string Configor::DataStream::ReferIMU = {};
std::string Configor::DataStream::BagPath = {};
double Configor::DataStream::BeginTime = {};
double Configor::DataStream::Duration = {};
std::string Configor::DataStream::OutputPath = {};
const std::string Configor::DataStream::PkgPath = ros::package::getPath("ikalibr");
const std::string Configor::DataStream::DebugPath = PkgPath + "/debug/";

std::string Configor::Prior::SpatTempPrioriPath = {};
double Configor::Prior::GravityNorm = {};
double Configor::Prior::TimeOffsetPadding = {};
double Configor::Prior::ReadoutTimePadding = {};
double Configor::Prior::KnotTimeDist::SO3Spline = {};
double Configor::Prior::KnotTimeDist::ScaleSpline = {};
double Configor::Prior::NDTLiDAROdometer::Resolution = {};
double Configor::Prior::NDTLiDAROdometer::KeyFrameDownSample = {};
double Configor::Prior::LiDARDataAssociate::MapDownSample = {};
double Configor::Prior::LiDARDataAssociate::PointToSurfelMax = {};
double Configor::Prior::LiDARDataAssociate::PlanarityMin = {};
double Configor::Prior::CauchyLossForRadarFactor = {};
double Configor::Prior::CauchyLossForLiDARFactor = {};
double Configor::Prior::CauchyLossForCameraFactor = {};
bool Configor::Prior::OptTemporalParams = {};

bool Configor::Preference::UseCudaInSolving = {};
OutputOption Configor::Preference::Outputs = OutputOption::NONE;
std::set<std::string> Configor::Preference::OutputsStr = {};
std::string Configor::Preference::OutputDataFormatStr = {};
CerealArchiveType::Enum Configor::Preference::OutputDataFormat = CerealArchiveType::Enum::YAML;
const std::map<CerealArchiveType::Enum, std::string> Configor::Preference::FileExtension = {
    {CerealArchiveType::Enum::YAML, ".yaml"},
    {CerealArchiveType::Enum::JSON, ".json"},
    {CerealArchiveType::Enum::XML, ".xml"},
    {CerealArchiveType::Enum::BINARY, ".bin"}};
int Configor::Preference::ThreadsToUse = {};
const std::string Configor::Preference::SO3_SPLINE = "SO3_SPLINE";
const std::string Configor::Preference::SCALE_SPLINE = "SCALE_SPLINE";
double Configor::Preference::SplineScaleInViewer = {};
double Configor::Preference::CoordSScaleInViewer = {};

std::optional<std::string> Configor::DataStream::CreateImageStoreFolder(
    const std::string &camTopic) {
    std::string path =
        ns_ikalibr::Configor::DataStream::OutputPath + "/images/" + camTopic + "/images";
    if (!std::filesystem::exists(path)) {
        if (!std::filesystem::create_directories(path)) {
            return {};
        }
    }
    return std::filesystem::canonical(path);
}

std::optional<std::string> Configor::DataStream::CreateSfMWorkspace(const std::string &camTopic) {
    std::string path =
        ns_ikalibr::Configor::DataStream::OutputPath + "/images/" + camTopic + "/sfm_ws";
    if (!std::filesystem::exists(path)) {
        if (!std::filesystem::create_directories(path)) {
            return {};
        }
    }
    return std::filesystem::canonical(path);
}

std::string Configor::DataStream::GetImageStoreInfoFile(const std::string &camTopic) {
    return ns_ikalibr::Configor::DataStream::OutputPath + "/images/" + camTopic + "/info" +
           Configor::GetFormatExtension();
}

int Configor::Preference::AvailableThreads() {
    int hardwareConcurrency = static_cast<int>(std::thread::hardware_concurrency());
    if (ThreadsToUse <= 0 || ThreadsToUse > hardwareConcurrency) {
        return hardwareConcurrency;
    } else {
        return ThreadsToUse;
    }
}

Configor::Configor() = default;

void Configor::PrintMainFields() {
    std::stringstream ssIMUTopics, ssRadarTopics, ssLiDARTopics, ssCameraTopics, ssRGBDTopics;

    for (const auto &[topic, _] : DataStream::IMUTopics) {
        ssIMUTopics << topic << " ";
    }
    for (const auto &[topic, _] : DataStream::RadarTopics) {
        ssRadarTopics << topic << " ";
    }
    for (const auto &[topic, _] : DataStream::LiDARTopics) {
        ssLiDARTopics << topic << " ";
    }
    for (const auto &[topic, _] : DataStream::CameraTopics) {
        ssCameraTopics << topic << " ";
    }
    for (const auto &[topic, info] : DataStream::RGBDTopics) {
        ssRGBDTopics << topic << ':' << info.DepthTopic << " ";
    }

    std::string IMUTopics = ssIMUTopics.str();
    std::string RadarTopics = ssRadarTopics.str();
    std::string LiDARTopics = ssLiDARTopics.str();
    std::string CameraTopics = ssCameraTopics.str();
    std::string RGBDTopics = ssRGBDTopics.str();

    auto GetOptString = [](OutputOption opt) -> std::string {
        std::stringstream stringStream;
        stringStream << magic_enum::enum_flags_name(opt);
        return stringStream.str();
    };

#define DESC_FIELD(field) #field, field
#define DESC_FORMAT "\n{:>45}: {}"
    spdlog::info(
        "main fields of configor:" DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT
            DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT
                DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT
                    DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT DESC_FORMAT
                        DESC_FORMAT DESC_FORMAT,
        DESC_FIELD(IMUTopics), DESC_FIELD(RadarTopics), DESC_FIELD(LiDARTopics),
        DESC_FIELD(CameraTopics), DESC_FIELD(RGBDTopics), DESC_FIELD(DataStream::ReferIMU),
        DESC_FIELD(DataStream::BagPath), DESC_FIELD(DataStream::BeginTime),
        DESC_FIELD(DataStream::Duration), DESC_FIELD(DataStream::OutputPath),
        DESC_FIELD(Prior::GravityNorm), DESC_FIELD(Prior::OptTemporalParams),
        DESC_FIELD(Prior::TimeOffsetPadding), DESC_FIELD(Prior::ReadoutTimePadding),
        DESC_FIELD(Prior::KnotTimeDist::SO3Spline), DESC_FIELD(Prior::KnotTimeDist::ScaleSpline),
        DESC_FIELD(Prior::NDTLiDAROdometer::Resolution),
        DESC_FIELD(Prior::NDTLiDAROdometer::KeyFrameDownSample),
        DESC_FIELD(Prior::LiDARDataAssociate::MapDownSample),
        DESC_FIELD(Prior::LiDARDataAssociate::PointToSurfelMax),
        DESC_FIELD(Prior::LiDARDataAssociate::PlanarityMin),
        DESC_FIELD(Prior::CauchyLossForRadarFactor), DESC_FIELD(Prior::CauchyLossForLiDARFactor),
        DESC_FIELD(Prior::CauchyLossForLiDARFactor), DESC_FIELD(Preference::UseCudaInSolving),
        "Preference::OutputDataFormat", Preference::OutputDataFormatStr, "Preference::Outputs",
        GetOptString(Preference::Outputs), DESC_FIELD(Preference::ThreadsToUse));

#undef DESC_FIELD
#undef DESC_FORMAT
}

void Configor::CheckConfigure() {
    // one imu need to be involved in ikalibr
    if (DataStream::IMUTopics.empty()) {
        throw Status(
            Status::ERROR,
            "the imu topic num (i.e., DataStream::IMUTopic) should be larger equal than 1!");
    }
    if (!Configor::IsLiDARIntegrated() && !Configor::IsRadarIntegrated() &&
        !Configor::IsCameraIntegrated() && !Configor::IsRGBDIntegrated() &&
        DataStream::IMUTopics.size() < 2) {
        throw Status(
            Status::ERROR,
            "performing multi-imu calibration requires imus that are more than or equal to 2!");
    }

    // check empty topics
    for (const auto &[topic, config] : DataStream::IMUTopics) {
        if (topic.empty()) {
            throw Status(Status::ERROR, "empty IMU topic exists!");
        }
        if (config.AcceWeight <= 0.0) {
            throw Status(Status::ERROR, "accelerator weight of IMU '{}' should be positive!",
                         topic);
        }
        if (config.GyroWeight <= 0.0) {
            throw Status(Status::ERROR, "gyroscope weight of IMU '{}' should be positive!", topic);
        }
        if (!std::filesystem::exists(config.Intrinsics)) {
            throw Status(Status::ERROR, "IMU intrinsic file for '{}' dose not exist: '{}'", topic,
                         config.Intrinsics);
        }
    }
    for (const auto &[topic, config] : DataStream::RadarTopics) {
        if (topic.empty()) {
            throw Status(Status::ERROR, "empty Radar topic exists!");
        }
        if (config.Weight <= 0.0) {
            throw Status(Status::ERROR, "weight of Radar '{}' should be positive!", topic);
        }
    }
    for (const auto &[topic, config] : DataStream::LiDARTopics) {
        if (topic.empty()) {
            throw Status(Status::ERROR, "empty LiDAR topic exists!");
        }
        if (config.Weight <= 0.0) {
            throw Status(Status::ERROR, "weight of LiDAR '{}' should be positive!", topic);
        }
    }
    for (const auto &[topic, config] : DataStream::CameraTopics) {
        if (topic.empty()) {
            throw Status(Status::ERROR, "empty camera topic exists!");
        }
        if (config.Weight <= 0.0) {
            throw Status(Status::ERROR, "weight of camera '{}' should be positive!", topic);
        }
        if (config.TrackLengthMin < 2) {
            throw Status(Status::ERROR, "track length of camera '{}' should be larger than '1'!",
                         topic);
        }
        if (!std::filesystem::exists(config.Intrinsics)) {
            throw Status(Status::ERROR, "camera intrinsic file for '{}' dose not exist: '{}'",
                         topic, config.Intrinsics);
        }
    }
    for (const auto &[topic, config] : DataStream::RGBDTopics) {
        if (topic.empty()) {
            throw Status(Status::ERROR, "empty color topic (rgbd) topic exists!");
        }
        if (config.DepthTopic.empty()) {
            throw Status(Status::ERROR, "empty depth topic (rgbd) topic exists!");
        }
        if (config.Weight <= 0.0) {
            throw Status(Status::ERROR, "weight of rgbd '{}' should be positive!", topic);
        }
        if (!std::filesystem::exists(config.Intrinsics)) {
            throw Status(Status::ERROR, "rgbd intrinsic file for '{}' dose not exist: '{}'", topic,
                         config.Intrinsics);
        }
    }

    // the reference imu should be one of multiple imus
    if (DataStream::IMUTopics.find(DataStream::ReferIMU) == DataStream::IMUTopics.cend()) {
        throw Status(Status::ERROR, "the reference IMU is not set, it should be one of the IMUs!");
    }

    if (!std::filesystem::exists(DataStream::BagPath)) {
        throw Status(Status::ERROR, "can not find the ros bag (i.e., DataStream::BagPath)!");
    }
    if (DataStream::OutputPath.empty()) {
        throw Status(Status::ERROR, "the output path (i.e., DataStream::OutputPath) is empty!");
    }
    if (!std::filesystem::exists(DataStream::OutputPath) &&
        !std::filesystem::create_directories(DataStream::OutputPath)) {
        // if the output path doesn't exist and create it failed
        throw Status(Status::ERROR,
                     "the output path (i.e., DataStream::OutputPath) can not be created!");
    }

    if (Prior::TimeOffsetPadding <= 0.0) {
        throw Status(
            Status::ERROR,
            "the time offset padding (i.e., Prior::TimeOffsetPadding) should be positive!");
    }
    if (Prior::ReadoutTimePadding <= 0.0) {
        throw Status(
            Status::ERROR,
            "the readout time padding (i.e., Prior::ReadoutTimePadding) should be positive!");
    }
    if (Prior::KnotTimeDist::SO3Spline <= 0.0) {
        throw Status(Status::ERROR,
                     "the knot time distance of so3 spline (i.e., Prior::KnotTimeDist::SO3Spline) "
                     "should be positive!");
    }
    if (Prior::KnotTimeDist::ScaleSpline <= 0.0) {
        throw Status(Status::ERROR,
                     "the knot time distance of scale spline (i.e., "
                     "Prior::KnotTimeDist::ScaleSpline) should be positive!");
    }
    if (Prior::NDTLiDAROdometer::Resolution <= 0.0) {
        throw Status(Status::ERROR,
                     "the resolution for NDT LiDAR odometer (i.e., "
                     "Prior::NDTLiDAROdometer::Resolution) should be positive!");
    }
    if (Prior::NDTLiDAROdometer::KeyFrameDownSample <= 0.0) {
        throw Status(Status::ERROR,
                     "the down sample rate for NDT LiDAR odometer (i.e., "
                     "Prior::NDTLiDAROdometer::KeyFrameDownSample) should be positive!");
    }

    if (Prior::CauchyLossForRadarFactor <= 0.0) {
        throw Status(Status::ERROR, "the Prior::CauchyLossForRadarFactor should be positive!");
    }
    if (Prior::CauchyLossForLiDARFactor <= 0.0) {
        throw Status(Status::ERROR, "the Prior::CauchyLossForLiDARFactor should be positive!");
    }
    if (Prior::CauchyLossForCameraFactor <= 0.0) {
        throw Status(Status::ERROR, "the Prior::CauchyLossForCameraFactor should be positive!");
    }

    if (Preference::SplineScaleInViewer <= 0.0) {
        throw Status(Status::ERROR, "the scale of splines in visualization should be positive!");
    }
    if (Preference::CoordSScaleInViewer <= 0.0) {
        throw Status(Status::ERROR,
                     "the scale of coordinates in visualization should be positive!");
    }
}

Configor::Ptr Configor::Create() { return std::make_shared<Configor>(); }

std::string Configor::GetFormatExtension() {
    return Preference::FileExtension.at(Preference::OutputDataFormat);
}

bool Configor::LoadConfigure(const std::string &filename, CerealArchiveType::Enum archiveType) {
    // load configure info
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    auto archive = GetInputArchiveVariant(file, archiveType);
    auto configor = Configor::Create();
    SerializeByInputArchiveVariant(archive, archiveType, cereal::make_nvp("Configor", *configor));

    // perform internal data transformation
    try {
        Configor::Preference::OutputDataFormat = EnumCast::stringToEnum<CerealArchiveType::Enum>(
            Configor::Preference::OutputDataFormatStr);
    } catch (...) {
        throw Status(Status::CRITICAL, "unsupported data format '{}' for io!!!",
                     Configor::Preference::OutputDataFormatStr);
    }
    for (const auto &output : Preference::OutputsStr) {
        // when the enum is out of range of [MAGIC_ENUM_RANGE_MIN, MAGIC_ENUM_RANGE_MAX],
        // magic_enum would not work
        // try {
        //     Configor::Preference::Outputs |= EnumCast::stringToEnum<OutputOption>(output);
        // } catch (...) {
        //     throw Status(Status::CRITICAL, "unsupported output context: '{}'!!!", output);
        // }
        if (auto iter = OutputOptionMap.find(output); iter == OutputOptionMap.cend()) {
            throw Status(Status::CRITICAL, "unsupported output context: '{}'!!!", output);
        } else {
            Configor::Preference::Outputs |= iter->second;
        }
    }

    // perform checking
    configor->CheckConfigure();
    return true;
}

bool Configor::SaveConfigure(const std::string &filename, CerealArchiveType::Enum archiveType) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    auto archive = GetOutputArchiveVariant(file, archiveType);
    SerializeByOutputArchiveVariant(archive, archiveType, cereal::make_nvp("Configor", *this));
    return true;
}

bool Configor::IsLiDARIntegrated() { return !DataStream::LiDARTopics.empty(); }

bool Configor::IsRadarIntegrated() { return !DataStream::RadarTopics.empty(); }

bool Configor::IsCameraIntegrated() { return !DataStream::CameraTopics.empty(); }

bool Configor::IsRGBDIntegrated() { return !DataStream::RGBDTopics.empty(); }
}  // namespace ns_ikalibr
