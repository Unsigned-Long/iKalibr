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

#include "calib/calib_param_manager.h"
#include "util/status.hpp"
#include "util/utils_tpl.hpp"
#include "spdlog/fmt/bundled/color.h"
#include "opencv2/calib3d.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "spdlog/spdlog.h"
#include "sensor/lidar_data_loader.h"
#include "tiny-viewer/object/radar.h"
#include "tiny-viewer/object/lidar.h"
#include "tiny-viewer/object/imu.h"
#include "tiny-viewer/object/camera.h"
#include "tiny-viewer/core/pose.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
// ------------------------
// static initialized filed
// ------------------------
const std::string CalibParamManager::Header::Software = "iKalibr";
const std::string CalibParamManager::Header::Version = "2.0.0";
const std::string CalibParamManager::Header::Address =
    "https://github.com/Unsigned-Long/iKalibr.git";

CalibParamManager::CalibParamManager(const std::vector<std::string> &imuTopics,
                                     const std::vector<std::string> &radarTopics,
                                     const std::vector<std::string> &lidarTopics,
                                     const std::vector<std::string> &cameraTopics,
                                     const std::vector<std::string> &rgbdTopics)
    : EXTRI(),
      TEMPORAL(),
      INTRI(),
      GRAVITY() {
    for (const auto &topic : imuTopics) {
        EXTRI.SO3_BiToBr[topic] = Sophus::SO3d();
        EXTRI.POS_BiInBr[topic] = Eigen::Vector3d::Zero();
        TEMPORAL.TO_BiToBr[topic] = 0.0;
        INTRI.IMU[topic] = IMUIntrinsics::Create();
    }
    for (const auto &topic : radarTopics) {
        EXTRI.SO3_RjToBr[topic] = Sophus::SO3d();
        EXTRI.POS_RjInBr[topic] = Eigen::Vector3d::Zero();
        TEMPORAL.TO_RjToBr[topic] = 0.0;
    }
    for (const auto &topic : lidarTopics) {
        EXTRI.SO3_LkToBr[topic] = Sophus::SO3d();
        EXTRI.POS_LkInBr[topic] = Eigen::Vector3d::Zero();
        TEMPORAL.TO_LkToBr[topic] = 0.0;
    }
    for (const auto &topic : cameraTopics) {
        EXTRI.SO3_CmToBr[topic] = Sophus::SO3d();
        EXTRI.POS_CmInBr[topic] = Eigen::Vector3d::Zero();
        TEMPORAL.TO_CmToBr[topic] = 0.0;
        TEMPORAL.RS_READOUT[topic] = 0.0;
        INTRI.Camera[topic] = nullptr;
    }
    for (const auto &topic : rgbdTopics) {
        EXTRI.SO3_DnToBr[topic] = Sophus::SO3d();
        EXTRI.POS_DnInBr[topic] = Eigen::Vector3d::Zero();
        TEMPORAL.TO_DnToBr[topic] = 0.0;
        TEMPORAL.RS_READOUT[topic] = 0.0;
        INTRI.RGBD[topic] = nullptr;
    }
    GRAVITY = Eigen::Vector3d(0.0, 0.0, -9.8);
}

CalibParamManager::Ptr CalibParamManager::Create(const std::vector<std::string> &imuTopics,
                                                 const std::vector<std::string> &radarTopics,
                                                 const std::vector<std::string> &lidarTopics,
                                                 const std::vector<std::string> &cameraTopics,
                                                 const std::vector<std::string> &rgbdTopics) {
    return std::make_shared<CalibParamManager>(imuTopics, radarTopics, lidarTopics, cameraTopics,
                                               rgbdTopics);
}

CalibParamManager::Ptr CalibParamManager::InitParamsFromConfigor() {
    spdlog::info("initialize calibration parameter manager using configor...");

    auto parMarg = CalibParamManager::Create(ExtractKeysAsVec(Configor::DataStream::IMUTopics),
                                             ExtractKeysAsVec(Configor::DataStream::RadarTopics),
                                             ExtractKeysAsVec(Configor::DataStream::LiDARTopics),
                                             ExtractKeysAsVec(Configor::DataStream::CameraTopics),
                                             ExtractKeysAsVec(Configor::DataStream::RGBDTopics));

    // intrinsics
    for (auto &[topic, intri] : parMarg->INTRI.IMU) {
        intri = ParIntri::LoadIMUIntri(Configor::DataStream::IMUTopics.at(topic).Intrinsics,
                                       Configor::Preference::OutputDataFormat);
    }
    for (auto &[topic, intri] : parMarg->INTRI.Camera) {
        intri = ParIntri::LoadCameraIntri(Configor::DataStream::CameraTopics.at(topic).Intrinsics,
                                          Configor::Preference::OutputDataFormat);
    }
    for (auto &[topic, intri] : parMarg->INTRI.RGBD) {
        intri = RGBDIntrinsics::Create(
            ParIntri::LoadCameraIntri(Configor::DataStream::RGBDTopics.at(topic).Intrinsics,
                                      Configor::Preference::OutputDataFormat),
            std::abs(Configor::DataStream::RGBDTopics.at(topic).DepthFactor), 0.0);
    }

    // align to the negative 'z' axis
    parMarg->GRAVITY = Eigen::Vector3d(0.0, 0.0, -Configor::Prior::GravityNorm);

    spdlog::info("initialize calibration parameter manager using configor finished.");
    return parMarg;
}

void CalibParamManager::ShowParamStatus() {
    std::stringstream stream;
#define ITEM(name) fmt::format(fmt::emphasis::bold | fmt::fg(fmt::color::green), name)
#define PARAM(name) fmt::format(fmt::emphasis::bold | fmt::fg(fmt::color::steel_blue), name)
#define STREAM_PACK(obj) stream << "-- " << obj << std::endl;

    constexpr std::size_t n = 74;

    STREAM_PACK(std::string(25, '-'))
    STREAM_PACK(ITEM("calibration parameters") << " --")
    STREAM_PACK(std::string(n, '-'))

    // -------------------------
    STREAM_PACK(ITEM("EXTRI"))
    // -------------------------
    STREAM_PACK("")

#define OUTPUT_EXTRINSICS(SENSOR1, IDX1, SENSOR2, IDX2)                                           \
    const auto EULER = EXTRI.EULER_##SENSOR1##IDX1##To##SENSOR2##IDX2##_DEG(topic);               \
    STREAM_PACK(PARAM("EULER_" #SENSOR1 #IDX1 "To" #SENSOR2 #IDX2 ": ")                           \
                << FormatValueVector<double>({"Xr", "Yp", "Zy"}, {EULER(0), EULER(1), EULER(2)})) \
                                                                                                  \
    const auto POS = EXTRI.POS_##SENSOR1##IDX1##In##SENSOR2##IDX2.at(topic);                      \
    STREAM_PACK(PARAM("  POS_" #SENSOR1 #IDX1 "In" #SENSOR2 #IDX2 ": ")                           \
                << FormatValueVector<double>({"Px", "Py", "Pz"}, {POS(0), POS(1), POS(2)}))

    // imus
    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        STREAM_PACK("IMU: '" << topic << "'")
        OUTPUT_EXTRINSICS(B, i, B, r)
        STREAM_PACK("")
    }

    // radars
    for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
        STREAM_PACK("Radar: '" << topic << "'")
        OUTPUT_EXTRINSICS(R, j, B, r)
        STREAM_PACK("")
    }

    // lidars
    for (const auto &[topic, _] : Configor::DataStream::LiDARTopics) {
        STREAM_PACK("LiDAR: '" << topic << "'")
        OUTPUT_EXTRINSICS(L, k, B, r)
        STREAM_PACK("")
    }

    // cameras
    for (const auto &[topic, _] : Configor::DataStream::CameraTopics) {
        STREAM_PACK("Camera: '" << topic << "'")
        OUTPUT_EXTRINSICS(C, m, B, r)
        STREAM_PACK("")
    }

    // rgbds
    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        STREAM_PACK("RGBD: '" << topic << "'")
        OUTPUT_EXTRINSICS(D, n, B, r)
        STREAM_PACK("")
    }

    STREAM_PACK(std::string(n, '-'))

    // ----------------------------
    STREAM_PACK(ITEM("TEMPORAL"))
    // ----------------------------
    STREAM_PACK("")

#define OUTPUT_TEMPORAL(SENSOR1, IDX1, SENSOR2, IDX2)                                           \
    const auto TO = TEMPORAL.TO_##SENSOR1##IDX1##To##SENSOR2##IDX2.at(topic);                   \
    STREAM_PACK(                                                                                \
        fmt::format("{}: {:+011.6f} (s)", PARAM("TO_" #SENSOR1 #IDX1 "To" #SENSOR2 #IDX2), TO)) \
                                                                                                \
    // imus
    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        STREAM_PACK("IMU: '" << topic << "'")
        OUTPUT_TEMPORAL(B, i, B, r)
        STREAM_PACK("")
    }

    // radars
    for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
        STREAM_PACK("Radar: '" << topic << "'")
        OUTPUT_TEMPORAL(R, j, B, r)
        STREAM_PACK("")
    }

    // lidars
    for (const auto &[topic, _] : Configor::DataStream::LiDARTopics) {
        STREAM_PACK("LiDAR: '" << topic << "'")
        OUTPUT_TEMPORAL(L, k, B, r)
        STREAM_PACK("")
    }

    // cameras
    for (const auto &[topic, _] : Configor::DataStream::CameraTopics) {
        STREAM_PACK("Camera: '" << topic << "'")
        OUTPUT_TEMPORAL(C, m, B, r)
        const auto RS_READOUT = TEMPORAL.RS_READOUT.at(topic);
        STREAM_PACK(fmt::format("{}: {:+010.6f} (s)", PARAM("RS_READOUT"), RS_READOUT))
        STREAM_PACK("")
    }

    // rgbds
    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        STREAM_PACK("RGBD: '" << topic << "'")
        OUTPUT_TEMPORAL(D, n, B, r)
        const auto RS_READOUT = TEMPORAL.RS_READOUT.at(topic);
        STREAM_PACK(fmt::format("{}: {:+010.6f} (s)", PARAM("RS_READOUT"), RS_READOUT))
        STREAM_PACK("")
    }

    STREAM_PACK(std::string(n, '-'))

    // -------------------------
    STREAM_PACK(ITEM("INTRI"))
    // -------------------------
    STREAM_PACK("")

    // imus
    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        STREAM_PACK("IMU: '" << topic << "'")
        const auto &ACCE = INTRI.IMU.at(topic)->ACCE;
        const auto &GYRO = INTRI.IMU.at(topic)->GYRO;
        // imu
        STREAM_PACK(PARAM("ACCE      BIAS: ") << FormatValueVector<double>(
                        {"Bx", "By", "Bz"}, {ACCE.BIAS(0), ACCE.BIAS(1), ACCE.BIAS(2)}))
        STREAM_PACK(
            PARAM("ACCE MAP COEFF: ") << FormatValueVector<double>(
                {"00", "11", "22"}, {ACCE.MAP_COEFF(0), ACCE.MAP_COEFF(1), ACCE.MAP_COEFF(2)}))
        STREAM_PACK(
            PARAM("                ") << FormatValueVector<double>(
                {"01", "02", "12"}, {ACCE.MAP_COEFF(3), ACCE.MAP_COEFF(4), ACCE.MAP_COEFF(5)}))

        STREAM_PACK("")

        STREAM_PACK(PARAM("GYRO      BIAS: ") << FormatValueVector<double>(
                        {"Bx", "By", "Bz"}, {GYRO.BIAS(0), GYRO.BIAS(1), GYRO.BIAS(2)}))
        STREAM_PACK(
            PARAM("GYRO MAP COEFF: ") << FormatValueVector<double>(
                {"00", "11", "22"}, {GYRO.MAP_COEFF(0), GYRO.MAP_COEFF(1), GYRO.MAP_COEFF(2)}))
        STREAM_PACK(
            PARAM("                ") << FormatValueVector<double>(
                {"01", "02", "12"}, {GYRO.MAP_COEFF(3), GYRO.MAP_COEFF(4), GYRO.MAP_COEFF(5)}))

        STREAM_PACK("")

        const auto EULER_AtoG = INTRI.IMU.at(topic)->EULER_AtoG_DEG();
        STREAM_PACK(PARAM("EULER AtoG DEG: ") << FormatValueVector<double>(
                        {"Xr", "Yp", "Zy"}, {EULER_AtoG(0), EULER_AtoG(1), EULER_AtoG(2)}))
        STREAM_PACK("")
    }

    // cameras
    for (const auto &[topic, _] : Configor::DataStream::CameraTopics) {
        STREAM_PACK("Camera: '" << topic << "'")
        const auto &intri = INTRI.Camera.at(topic);
        const auto &pars = intri->GetParams();

        STREAM_PACK(
            PARAM("IMAGE     SIZE: ") << FormatValueVector<double>(
                {" w", " h"}, {(double)intri->imgWidth, (double)intri->imgHeight}, "{:+011.5f}"))

        STREAM_PACK(PARAM("FOCAL   LENGTH: ")
                    << FormatValueVector<double>({"fx", "fy"}, {pars.at(0), pars.at(1)}))

        STREAM_PACK(PARAM("PRINCIP  POINT: ")
                    << FormatValueVector<double>({"cx", "cy"}, {pars.at(2), pars.at(3)}))

        STREAM_PACK("")
        if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri)) {
            STREAM_PACK(PARAM("DISTO   PARAMS: ") << FormatValueVector<double>(
                            {"k1", "k2", "k3"}, {pars.at(4), pars.at(5), pars.at(6)}))
            STREAM_PACK(PARAM("                ")
                        << FormatValueVector<double>({"p1", "p2"}, {pars.at(7), pars.at(8)}))
        } else if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicFisheye>(intri)) {
            STREAM_PACK(PARAM("DISTO   PARAMS: ")
                        << FormatValueVector<double>({"k1", "k2"}, {pars.at(4), pars.at(5)}))
            STREAM_PACK(PARAM("                ")
                        << FormatValueVector<double>({"k3", "k4"}, {pars.at(6), pars.at(7)}))
        } else {
            throw Status(Status::CRITICAL,
                         "unknown camera intrinsic model! supported models:\n"
                         "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                         "(b)  pinhole_fisheye (k1, k2, k3, k4)");
        }
        STREAM_PACK("")
    }

    // rgbds
    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        STREAM_PACK("RGBD: '" << topic << "'")
        const auto &rgbdIntri = INTRI.RGBD.at(topic);
        const auto &intri = rgbdIntri->intri;
        const auto &pars = intri->GetParams();

        STREAM_PACK(
            PARAM("IMAGE     SIZE: ") << FormatValueVector<double>(
                {" w", " h"}, {(double)intri->imgWidth, (double)intri->imgHeight}, "{:+011.5f}"))

        STREAM_PACK(PARAM("FOCAL   LENGTH: ")
                    << FormatValueVector<double>({"fx", "fy"}, {pars.at(0), pars.at(1)}))

        STREAM_PACK(PARAM("PRINCIP  POINT: ")
                    << FormatValueVector<double>({"cx", "cy"}, {pars.at(2), pars.at(3)}))

        STREAM_PACK(PARAM("DEPTH   FACTOR: ")
                    << FormatValueVector<double>({" a", " b"}, {rgbdIntri->alpha, rgbdIntri->beta}))

        STREAM_PACK("")
        if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri)) {
            STREAM_PACK(PARAM("DISTO   PARAMS: ") << FormatValueVector<double>(
                            {"k1", "k2", "k3"}, {pars.at(4), pars.at(5), pars.at(6)}))
            STREAM_PACK(PARAM("                ")
                        << FormatValueVector<double>({"p1", "p2"}, {pars.at(7), pars.at(8)}))
        } else if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicFisheye>(intri)) {
            STREAM_PACK(PARAM("DISTO   PARAMS: ")
                        << FormatValueVector<double>({"k1", "k2"}, {pars.at(4), pars.at(5)}))
            STREAM_PACK(PARAM("                ")
                        << FormatValueVector<double>({"k3", "k4"}, {pars.at(6), pars.at(7)}))
        } else {
            throw Status(Status::CRITICAL,
                         "unknown camera intrinsic model! supported models:\n"
                         "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                         "(b)  pinhole_fisheye (k1, k2, k3, k4)");
        }
        STREAM_PACK("")
    }

    STREAM_PACK(std::string(n, '-'))

    // ------------------------------
    STREAM_PACK(ITEM("OTHER FIELDS"))
    // ------------------------------
    STREAM_PACK("")
    STREAM_PACK(PARAM("GRAVITY IN MAP: ") << FormatValueVector<double>(
                    {"Gx", "Gy", "Gz"}, {GRAVITY(0), GRAVITY(1), GRAVITY(2)}))

    STREAM_PACK(std::string(n, '-'))

    spdlog::info("the detail calibration parameters are below: \n{}", stream.str());

#undef ITEM
#undef PARAM
}

std::vector<std::size_t> CalibParamManager::VisualizationSensors(ns_viewer::Viewer &viewer) const {
    return viewer.AddEntity(EntitiesForVisualization());
}

std::vector<std::size_t> CalibParamManager::VisualizationSensors(ns_viewer::MultiViewer &viewer,
                                                                 const std::string &win) const {
    return viewer.AddEntity(EntitiesForVisualization(), win);
}

void CalibParamManager::Save(const std::string &filename,
                             CerealArchiveType::Enum archiveType) const {
    std::ofstream file(filename, std::ios::out);
    auto ar = GetOutputArchiveVariant(file, archiveType);
    SerializeByOutputArchiveVariant(ar, archiveType, cereal::make_nvp("CalibParam", *this));
}

CalibParamManager::Ptr CalibParamManager::Load(const std::string &filename,
                                               CerealArchiveType::Enum archiveType) {
    auto calibParamManager = CalibParamManager::Create();
    std::ifstream file(filename, std::ios::in);
    auto ar = GetInputArchiveVariant(file, archiveType);
    SerializeByInputArchiveVariant(ar, archiveType,
                                   cereal::make_nvp("CalibParam", *calibParamManager));
    return calibParamManager;
}

std::vector<ns_viewer::Entity::Ptr> CalibParamManager::EntitiesForVisualization() const {
#define IMU_SIZE 0.02
#define RADAR_SIZE 0.05
#define LiDAR_SIZE 0.04
#define CAMERA_SIZE 0.04
#define RGBD_SIZE CAMERA_SIZE

    std::vector<ns_viewer::Entity::Ptr> entities;

    // reference imu
    auto SE3_BrToBr = Sophus::SE3f();
    auto refIMU = ns_viewer::IMU::Create(
        ns_viewer::Posef(SE3_BrToBr.so3().matrix(), SE3_BrToBr.translation()), IMU_SIZE,
        ns_viewer::Colour(0.3f, 0.3f, 0.3f, 1.0f));
    entities.push_back(refIMU);

    // imus
    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        auto SE3_BiToBr = EXTRI.SE3_BiToBr(topic).cast<float>();
        auto imu = ns_viewer::IMU::Create(
            ns_viewer::Posef(SE3_BiToBr.so3().matrix(), SE3_BiToBr.translation()), IMU_SIZE,
            ns_viewer::Colour::Red());
        auto line =
            ns_viewer::Line::Create(Eigen::Vector3f::Zero(), SE3_BiToBr.translation().cast<float>(),
                                    ns_viewer::Colour::Black());
        entities.push_back(imu);
        entities.push_back(line);
    }

    // radars
    for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
        auto SE3_RjToBr = EXTRI.SE3_RjToBr(topic).cast<float>();
        auto radar = ns_viewer::Radar::Create(
            ns_viewer::Posef(SE3_RjToBr.so3().matrix(), SE3_RjToBr.translation()), RADAR_SIZE,
            ns_viewer::Colour::Green());
        auto line =
            ns_viewer::Line::Create(Eigen::Vector3f::Zero(), SE3_RjToBr.translation().cast<float>(),
                                    ns_viewer::Colour::Black());
        entities.push_back(radar);
        entities.push_back(line);
    }

    // lidars
    for (const auto &[topic, _] : Configor::DataStream::LiDARTopics) {
        auto SE3_LkToBr = EXTRI.SE3_LkToBr(topic).cast<float>();
        auto line =
            ns_viewer::Line::Create(Eigen::Vector3f::Zero(), SE3_LkToBr.translation().cast<float>(),
                                    ns_viewer::Colour::Black());
        entities.push_back(line);
        ns_viewer::Entity::Ptr lidar;
        if (EnumCast::stringToEnum<LidarModelType>(
                Configor::DataStream::LiDARTopics.at(topic).Type) == LidarModelType::LIVOX_CUSTOM) {
            lidar = ns_viewer::LivoxLiDAR::Create(
                ns_viewer::Posef(SE3_LkToBr.so3().matrix(), SE3_LkToBr.translation()), LiDAR_SIZE,
                ns_viewer::Colour(0.33f, 0.33f, 0.5f, 1.0f));
        } else {
            lidar = ns_viewer::LiDAR::Create(
                ns_viewer::Posef(SE3_LkToBr.so3().matrix(), SE3_LkToBr.translation()), LiDAR_SIZE,
                ns_viewer::Colour(0.33f, 0.33f, 0.5f, 1.0f));
        }
        entities.push_back(lidar);
    }

    // pos cameras
    for (const auto &[topic, _] : Configor::DataStream::PosCameraTopics()) {
        auto SE3_CmToBr = EXTRI.SE3_CmToBr(topic).cast<float>();
        auto camera = ns_viewer::CubeCamera::Create(
            ns_viewer::Posef(SE3_CmToBr.so3().matrix(), SE3_CmToBr.translation()), CAMERA_SIZE,
            ns_viewer::Colour::Blue());
        auto line =
            ns_viewer::Line::Create(Eigen::Vector3f::Zero(), SE3_CmToBr.translation().cast<float>(),
                                    ns_viewer::Colour::Black());
        entities.push_back(camera);
        entities.push_back(line);
    }

    // vel cameras
    for (const auto &[topic, _] : Configor::DataStream::VelCameraTopics()) {
        auto SE3_CmToBr = EXTRI.SE3_CmToBr(topic).cast<float>();
        auto camera = ns_viewer::CubeCamera::Create(
            ns_viewer::Posef(SE3_CmToBr.so3().matrix(), SE3_CmToBr.translation()), CAMERA_SIZE,
            ns_viewer::Colour::Green());
        auto line =
            ns_viewer::Line::Create(Eigen::Vector3f::Zero(), SE3_CmToBr.translation().cast<float>(),
                                    ns_viewer::Colour::Black());
        entities.push_back(camera);
        entities.push_back(line);
    }

    // rgbds
    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        auto SE3_DnToBr = EXTRI.SE3_DnToBr(topic).cast<float>();
        auto rgbd = ns_viewer::CubeCamera::Create(
            ns_viewer::Posef(SE3_DnToBr.so3().matrix(), SE3_DnToBr.translation()), RGBD_SIZE,
            ns_viewer::Colour(1.0f, 0.5f, 0.0f, 1.0f));
        auto line =
            ns_viewer::Line::Create(Eigen::Vector3f::Zero(), SE3_DnToBr.translation().cast<float>(),
                                    ns_viewer::Colour::Black());
        entities.push_back(rgbd);
        entities.push_back(line);
    }

#undef IMU_SIZE
#undef RADAR_SIZE
#undef LiDAR_SIZE
#undef CAMERA_SIZE
#undef RGBD_SIZE

    return entities;
}

ns_veta::PinholeIntrinsic::Ptr CalibParamManager::ParIntri::LoadCameraIntri(
    const std::string &filename, CerealArchiveType::Enum archiveType) {
    auto intri = ns_veta::PinholeIntrinsic::Create(0, 0, 0, 0, 0, 0);
    std::ifstream file(filename, std::ios::in);
    auto ar = GetInputArchiveVariant(file, archiveType);
    try {
        SerializeByInputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", intri));
    } catch (const cereal::Exception &exception) {
        throw Status(Status::CRITICAL,
                     "The configuration file '{}' for 'PinholeIntrinsic' is "
                     "outdated or broken, and can not be loaded in iKalibr using cereal!!! "
                     "To make it right, please refer to our latest configuration file "
                     "template released at "
                     "https://github.com/Unsigned-Long/iKalibr/blob/master/config/"
                     "cam-intri-pinhole-fisheye.yaml and "
                     "https://github.com/Unsigned-Long/iKalibr/blob/master/config/"
                     "cam-intri-pinhole-brown.yaml, and then fix your custom configuration file. "
                     "Detailed cereal exception information: \n'{}'",
                     filename, exception.what());
    }
    return intri;
}

IMUIntrinsics::Ptr CalibParamManager::ParIntri::LoadIMUIntri(const std::string &filename,
                                                             CerealArchiveType::Enum archiveType) {
    auto intri = IMUIntrinsics::Create();
    std::ifstream file(filename, std::ios::in);
    auto ar = GetInputArchiveVariant(file, archiveType);
    try {
        SerializeByInputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", intri));
    } catch (const cereal::Exception &exception) {
        throw Status(Status::CRITICAL,
                     "The configuration file '{}' for 'IMUIntrinsics' is "
                     "outdated or broken, and can not be loaded in iKalibr using cereal!!! "
                     "To make it right, please refer to our latest configuration file "
                     "template released at "
                     "https://github.com/Unsigned-Long/iKalibr/blob/master/config/imu-intri.yaml, "
                     "and then fix your custom configuration file. Detailed cereal exception "
                     "information: \n'{}'",
                     filename, exception.what());
    }
    return intri;
}

void CalibParamManager::ParIntri::SaveCameraIntri(const ns_veta::PinholeIntrinsic::Ptr &intri,
                                                  const std::string &filename,
                                                  CerealArchiveType::Enum archiveType) {
    std::ofstream file(filename, std::ios::out);
    auto ar = GetOutputArchiveVariant(file, archiveType);
    SerializeByOutputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", intri));
}

void CalibParamManager::ParIntri::SaveIMUIntri(const IMUIntrinsics::Ptr &intri,
                                               const std::string &filename,
                                               CerealArchiveType::Enum archiveType) {
    std::ofstream file(filename, std::ios::out);
    auto ar = GetOutputArchiveVariant(file, archiveType);
    SerializeByOutputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", intri));
}

cv::Mat CalibParamManager::ParIntri::UndistortImage(const ns_veta::PinholeIntrinsic::Ptr &intri,
                                                    const cv::Mat &src) {
    auto [K, D] = ObtainKDMatForUndisto(intri);
    if (cv::countNonZero(D) == 0) {
        // don't need to undedistort this image
        return src.clone();
    }
    cv::Mat undistImg;
    if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri)) {
        cv::undistort(src, undistImg, K, D);
    } else if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicFisheye>(intri)) {
        cv::Mat E = cv::Mat::eye(3, 3, cv::DataType<double>::type);
        cv::Size size(src.cols, src.rows);
        cv::Mat map1, map2;
        cv::fisheye::initUndistortRectifyMap(K, D, E, K, size, CV_16SC2, map1, map2);
        cv::remap(src, undistImg, map1, map2, cv::INTER_LINEAR, CV_HAL_BORDER_CONSTANT);
    } else {
        throw Status(Status::CRITICAL,
                     "unknown camera intrinsic model! supported models:\n"
                     "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                     "(b)  pinhole_fisheye (k1, k2, k3, k4)");
    }
    return undistImg;
}

std::pair<cv::Mat, cv::Mat> CalibParamManager::ParIntri::ObtainKDMatForUndisto(
    const ns_veta::PinholeIntrinsic::Ptr &intri) {
    const auto &par = intri->GetParams();
    cv::Mat K, D;
    K = (cv::Mat_<double>(3, 3) << par[0], 0.0, par[2], 0.0, par[1], par[3], 0.0, 0.0, 1.0);

    if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri)) {
        // k_1, k_2, p_1, p_2[, k_3[, k_4, k_5, k_6[, s_1, s_2, s_3, s_4[, t_x, t_y]]]]
        D = (cv::Mat_<double>(5, 1) << par[4], par[5], par[7], par[8], par[6]);
    } else if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicFisheye>(intri)) {
        // k_1, k_2, k_3, k_4
        D = (cv::Mat_<double>(4, 1) << par[4], par[5], par[6], par[7]);
    } else {
        throw Status(Status::CRITICAL,
                     "unknown camera intrinsic model! supported models:\n"
                     "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                     "(b)  pinhole_fisheye (k1, k2, k3, k4)");
    }

    return {K, D};
}
}  // namespace ns_ikalibr
