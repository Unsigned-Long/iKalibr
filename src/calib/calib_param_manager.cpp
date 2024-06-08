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

#include "filesystem"
#include "calib/calib_param_manager.h"
#include "spdlog/fmt/bundled/color.h"
#include "util/utils.hpp"
#include <algorithm>
#include "opencv2/calib3d.hpp"
#include <opencv2/imgproc/imgproc.hpp>

_3_

namespace ns_ikalibr {

    CalibParamManager::CalibParamManager(
            const std::vector<std::string> &imuTopics,
            const std::vector<std::string> &radarTopics,
            const std::vector<std::string> &lidarTopics,
            const std::vector<std::string> &cameraTopics
    ) : EXTRI(), TEMPORAL(), INTRI(), GRAVITY() {
        for (const auto &topic: imuTopics) {
            EXTRI.SO3_BiToBr[topic] = Sophus::SO3d();
            EXTRI.POS_BiInBr[topic] = Eigen::Vector3d::Zero();
            TEMPORAL.TO_BiToBr[topic] = 0.0;
            INTRI.IMU[topic] = IMUIntrinsics::Create();
        }
        for (const auto &topic: radarTopics) {
            EXTRI.SO3_RjToBr[topic] = Sophus::SO3d();
            EXTRI.POS_RjInBr[topic] = Eigen::Vector3d::Zero();
            TEMPORAL.TO_RjToBr[topic] = 0.0;
        }
        for (const auto &topic: lidarTopics) {
            EXTRI.SO3_LkToBr[topic] = Sophus::SO3d();
            EXTRI.POS_LkInBr[topic] = Eigen::Vector3d::Zero();
            TEMPORAL.TO_LkToBr[topic] = 0.0;
        }
        for (const auto &topic: cameraTopics) {
            EXTRI.SO3_CmToBr[topic] = Sophus::SO3d();
            EXTRI.POS_CmInBr[topic] = Eigen::Vector3d::Zero();
            TEMPORAL.TO_CmToBr[topic] = 0.0;
            TEMPORAL.RS_READOUT[topic] = 0.0;
            INTRI.Camera[topic] = nullptr;
        }
        GRAVITY = Eigen::Vector3d(0.0, 0.0, -9.8);
    }

    CalibParamManager::Ptr CalibParamManager::Create(const std::vector<std::string> &imuTopics,
                                                     const std::vector<std::string> &radarTopics,
                                                     const std::vector<std::string> &lidarTopics,
                                                     const std::vector<std::string> &cameraTopics) {
        return std::make_shared<CalibParamManager>(imuTopics, radarTopics, lidarTopics, cameraTopics);
    }

    CalibParamManager::Ptr CalibParamManager::InitParamsFromConfigor() {
        spdlog::info("initialize calibration parameter manager using configor.");

        auto parMarg = CalibParamManager::Create(
                ExtractKeysAsVec(Configor::DataStream::IMUTopics),
                ExtractKeysAsVec(Configor::DataStream::RadarTopics),
                ExtractKeysAsVec(Configor::DataStream::LiDARTopics),
                ExtractKeysAsVec(Configor::DataStream::CameraTopics)
        );

        // intrinsics
        for (auto &[topic, intri]: parMarg->INTRI.IMU) {
            intri = ParIntri::LoadIMUIntri(
                    Configor::DataStream::IMUTopics.at(topic).Intrinsics, Configor::Preference::DataIOFormat()
            );
        }
        for (auto &[topic, intri]: parMarg->INTRI.Camera) {
            intri = ParIntri::LoadCameraIntri(
                    Configor::DataStream::CameraTopics.at(topic).Intrinsics, Configor::Preference::DataIOFormat()
            );
        }

        // align to the negative 'z' axis
        parMarg->GRAVITY = Eigen::Vector3d(0.0, 0.0, -Configor::Prior::GravityNorm);

        return parMarg;
    }

    void CalibParamManager::ShowParamStatus() {
        std::stringstream stream;
#define ITEM(name) fmt::format(fmt::emphasis::bold | fmt::fg(fmt::color::green), name)
#define PARAM(name) fmt::format(fmt::emphasis::bold | fmt::fg(fmt::color::black), name)
#define STREAM_PACK(obj) stream << "-- " << obj << std::endl;

        constexpr std::size_t n = 74;

        STREAM_PACK(std::string(25, '-'))
        STREAM_PACK(ITEM("calibration parameters") << " --")
        STREAM_PACK(std::string(n, '-'))

        // -------------------------
        STREAM_PACK(ITEM("EXTRI"))
        // -------------------------
        STREAM_PACK("")

        for (const auto &[topic, _]: EXTRI.SO3_BiToBr) {
            STREAM_PACK("IMU: '" << topic << "'")
            const auto EULER_BiToBr = EXTRI.EULER_BiToBr_DEG(topic);
            STREAM_PACK(PARAM("EULER_BiToBr: ") << FormatValueVector<double>(
                    {"Xr", "Yp", "Zy"}, {EULER_BiToBr(0), EULER_BiToBr(1), EULER_BiToBr(2)}))

            const auto POS_BiInBc = EXTRI.POS_BiInBr.at(topic);
            STREAM_PACK(PARAM("  POS_BiInBr: ") << FormatValueVector<double>(
                    {"Px", "Py", "Pz"}, {POS_BiInBc(0), POS_BiInBc(1), POS_BiInBc(2)}))
            STREAM_PACK("")
        }

        for (const auto &[topic, _]: EXTRI.SO3_RjToBr) {
            STREAM_PACK("Radar: '" << topic << "'")
            const auto EULER_RjToBr = EXTRI.EULER_RjToBr_DEG(topic);
            STREAM_PACK(PARAM("EULER_RjToBr: ") << FormatValueVector<double>(
                    {"Xr", "Yp", "Zy"}, {EULER_RjToBr(0), EULER_RjToBr(1), EULER_RjToBr(2)}))

            const auto POS_RjInBr = EXTRI.POS_RjInBr.at(topic);
            STREAM_PACK(PARAM("  POS_RjInBr: ") << FormatValueVector<double>(
                    {"Px", "Py", "Pz"}, {POS_RjInBr(0), POS_RjInBr(1), POS_RjInBr(2)}))
            STREAM_PACK("")
        }

        for (const auto &[topic, _]: EXTRI.SO3_LkToBr) {
            STREAM_PACK("LiDAR: '" << topic << "'")
            const auto EULER_LkToBr = EXTRI.EULER_LkToBr_DEG(topic);
            STREAM_PACK(PARAM("EULER_LkToBr: ") << FormatValueVector<double>(
                    {"Xr", "Yp", "Zy"}, {EULER_LkToBr(0), EULER_LkToBr(1), EULER_LkToBr(2)}))

            const auto POS_LkInBr = EXTRI.POS_LkInBr.at(topic);
            STREAM_PACK(PARAM("  POS_LkInBr: ") << FormatValueVector<double>(
                    {"Px", "Py", "Pz"}, {POS_LkInBr(0), POS_LkInBr(1), POS_LkInBr(2)}))
            STREAM_PACK("")
        }

        for (const auto &[topic, _]: EXTRI.SO3_CmToBr) {
            STREAM_PACK("Camera: '" << topic << "'")
            const auto EULER_CmToBr = EXTRI.EULER_CmToBr_DEG(topic);
            STREAM_PACK(PARAM("EULER_CmToBr: ") << FormatValueVector<double>(
                    {"Xr", "Yp", "Zy"}, {EULER_CmToBr(0), EULER_CmToBr(1), EULER_CmToBr(2)}))

            const auto POS_CmInBr = EXTRI.POS_CmInBr.at(topic);
            STREAM_PACK(PARAM("  POS_CmInBr: ") << FormatValueVector<double>(
                    {"Px", "Py", "Pz"}, {POS_CmInBr(0), POS_CmInBr(1), POS_CmInBr(2)}))
            STREAM_PACK("")
        }

        STREAM_PACK(std::string(n, '-'))

        // ----------------------------
        STREAM_PACK(ITEM("TEMPORAL"))
        // ----------------------------
        STREAM_PACK("")

        for (const auto &[topic, _]: EXTRI.SO3_BiToBr) {
            STREAM_PACK("IMU: '" << topic << "'")
            const auto TO_BiToBr = TEMPORAL.TO_BiToBr.at(topic);
            STREAM_PACK(fmt::format("{}: {:+011.6f} (s)", PARAM("TO_BiToBr"), TO_BiToBr))
            STREAM_PACK("")
        }

        for (const auto &[topic, _]: EXTRI.SO3_RjToBr) {
            STREAM_PACK("Radar: '" << topic << "'")
            const auto TO_RjToBr = TEMPORAL.TO_RjToBr.at(topic);
            STREAM_PACK(fmt::format("{}: {:+011.6f} (s)", PARAM("TO_RjToBr"), TO_RjToBr))
            STREAM_PACK("")
        }

        for (const auto &[topic, _]: EXTRI.SO3_LkToBr) {
            STREAM_PACK("LiDAR: '" << topic << "'")
            const auto TO_LkToBr = TEMPORAL.TO_LkToBr.at(topic);
            STREAM_PACK(fmt::format("{}: {:+011.6f} (s)", PARAM("TO_LkToBr"), TO_LkToBr))
            STREAM_PACK("")
        }

        for (const auto &[topic, _]: EXTRI.SO3_CmToBr) {
            STREAM_PACK("Camera: '" << topic << "'")
            const auto TO_CmToBr = TEMPORAL.TO_CmToBr.at(topic);
            const auto RS_READOUT = TEMPORAL.RS_READOUT.at(topic);
            STREAM_PACK(fmt::format("{}: {:+011.6f} (s)", PARAM("TO_CmToBr"), TO_CmToBr))
            STREAM_PACK(fmt::format("{}: {:+010.6f} (s)", PARAM("RS_READOUT"), RS_READOUT))
            STREAM_PACK("")
        }


        STREAM_PACK(std::string(n, '-'))

        // -------------------------
        STREAM_PACK(ITEM("INTRI"))
        // -------------------------
        STREAM_PACK("")

        for (const auto &[topic, _]: EXTRI.SO3_BiToBr) {
            STREAM_PACK("IMU: '" << topic << "'")
            const auto &ACCE = INTRI.IMU.at(topic)->ACCE;
            const auto &GYRO = INTRI.IMU.at(topic)->GYRO;
            // imu
            STREAM_PACK(PARAM("ACCE      BIAS: ") << FormatValueVector<double>(
                    {"Bx", "By", "Bz"}, {ACCE.BIAS(0), ACCE.BIAS(1), ACCE.BIAS(2)}))
            STREAM_PACK(PARAM("ACCE MAP COEFF: ") << FormatValueVector<double>(
                    {"00", "11", "22"},
                    {ACCE.MAP_COEFF(0), ACCE.MAP_COEFF(1), ACCE.MAP_COEFF(2)}))
            STREAM_PACK(PARAM("                ") << FormatValueVector<double>(
                    {"01", "02", "12"}, {ACCE.MAP_COEFF(3), ACCE.MAP_COEFF(4), ACCE.MAP_COEFF(5)}))

            STREAM_PACK("")

            STREAM_PACK(PARAM("GYRO      BIAS: ") << FormatValueVector<double>(
                    {"Bx", "By", "Bz"}, {GYRO.BIAS(0), GYRO.BIAS(1), GYRO.BIAS(2)}))
            STREAM_PACK(PARAM("GYRO MAP COEFF: ") << FormatValueVector<double>(
                    {"00", "11", "22"},
                    {GYRO.MAP_COEFF(0), GYRO.MAP_COEFF(1), GYRO.MAP_COEFF(2)}))
            STREAM_PACK(PARAM("                ") << FormatValueVector<double>(
                    {"01", "02", "12"}, {GYRO.MAP_COEFF(3), GYRO.MAP_COEFF(4), GYRO.MAP_COEFF(5)}))

            STREAM_PACK("")

            const auto EULER_AtoG = INTRI.IMU.at(topic)->EULER_AtoG_DEG();
            STREAM_PACK(PARAM("EULER AtoG DEG: ") << FormatValueVector<double>(
                    {"Xr", "Yp", "Zy"}, {EULER_AtoG(0), EULER_AtoG(1), EULER_AtoG(2)}))
            STREAM_PACK("")
        }

        for (const auto &[topic, _]: EXTRI.SO3_CmToBr) {
            STREAM_PACK("Camera: '" << topic << "'")
            const auto &intri = INTRI.Camera.at(topic);
            const auto &pars = intri->GetParams();

            STREAM_PACK(PARAM("IMAGE     SIZE: ") << FormatValueVector<double>(
                    {" w", " h"}, {(double) intri->imgWidth, (double) intri->imgHeight}, "{:+011.5f}"))

            STREAM_PACK(PARAM("FOCAL   LENGTH: ") << FormatValueVector<double>(
                    {"fx", "fy"}, {pars.at(0), pars.at(1)}))

            STREAM_PACK(PARAM("PRINCIP  POINT: ") << FormatValueVector<double>(
                    {"cx", "cy"}, {pars.at(2), pars.at(3)}))

            STREAM_PACK("")
            if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri)) {
                STREAM_PACK(PARAM("DISTO   PARAMS: ") << FormatValueVector<double>(
                        {"k1", "k2", "k3"}, {pars.at(4), pars.at(5), pars.at(6)}))
                STREAM_PACK(PARAM("                ") << FormatValueVector<double>(
                        {"p1", "p2"}, {pars.at(7), pars.at(8)}))
            } else if (std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicFisheye>(intri)) {
                STREAM_PACK(PARAM("DISTO   PARAMS: ") << FormatValueVector<double>(
                        {"k1", "k2"}, {pars.at(4), pars.at(5)}))
                STREAM_PACK(PARAM("                ") << FormatValueVector<double>(
                        {"k3", "k4"}, {pars.at(6), pars.at(7)}))
            } else {
                throw Status(
                        Status::CRITICAL, "unknown camera intrinsic model! supported models:\n"
                                          "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                                          "(b)  pinhole_fisheye (k1, k2, k3, k4)"
                );
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

    std::vector<std::size_t>
    CalibParamManager::VisualizationSensors(ns_viewer::MultiViewer &viewer, const std::string &win) const {
        return viewer.AddEntity(EntitiesForVisualization(), win);
    }

    void CalibParamManager::Save(const std::string &filename, CerealArchiveType::Enum archiveType) const {
        std::ofstream file(filename, std::ios::out);
        auto ar = GetOutputArchiveVariant(file, archiveType);
        SerializeByOutputArchiveVariant(ar, archiveType, cereal::make_nvp("CalibParam", *this));
    }

    CalibParamManager::Ptr CalibParamManager::Load(const std::string &filename, CerealArchiveType::Enum archiveType) {
        auto calibParamManager = CalibParamManager::Create();
        std::ifstream file(filename, std::ios::in);
        auto ar = GetInputArchiveVariant(file, archiveType);
        SerializeByInputArchiveVariant(ar, archiveType, cereal::make_nvp("CalibParam", *calibParamManager));
        return calibParamManager;
    }

    std::vector<ns_viewer::Entity::Ptr> CalibParamManager::EntitiesForVisualization() const {
#define IMU_SIZE 0.02
#define RADAR_SIZE 0.05
#define LiDAR_SIZE 0.04
#define CAMERA_SIZE 0.04

        std::vector<ns_viewer::Entity::Ptr> entities;

        auto SE3_BrToBr = Sophus::SE3f();
        auto refIMU = ns_viewer::IMU::Create(
                ns_viewer::Posef(SE3_BrToBr.so3().matrix(), SE3_BrToBr.translation()), IMU_SIZE,
                ns_viewer::Colour(0.3f, 0.3f, 0.3f, 1.0f)
        );
        entities.push_back(refIMU);

        for (const auto &[topic, _]: CalibParamManager::EXTRI.SO3_BiToBr) {
            auto SE3_BiToBr = EXTRI.SE3_BiToBr(topic).cast<float>();
            auto imu = ns_viewer::IMU::Create(
                    ns_viewer::Posef(SE3_BiToBr.so3().matrix(), SE3_BiToBr.translation()),
                    IMU_SIZE, ns_viewer::Colour::Red()
            );
            auto line = ns_viewer::Line::Create(
                    Eigen::Vector3f::Zero(), SE3_BiToBr.translation().cast<float>(), ns_viewer::Colour::Black()
            );
            entities.push_back(imu);
            entities.push_back(line);
        }

        for (const auto &[topic, _]: CalibParamManager::EXTRI.SO3_RjToBr) {
            auto SE3_RjToBr = EXTRI.SE3_RjToBr(topic).cast<float>();
            auto radar = ns_viewer::Radar::Create(
                    ns_viewer::Posef(SE3_RjToBr.so3().matrix(), SE3_RjToBr.translation()),
                    RADAR_SIZE, ns_viewer::Colour::Green()
            );
            auto line = ns_viewer::Line::Create(
                    Eigen::Vector3f::Zero(), SE3_RjToBr.translation().cast<float>(), ns_viewer::Colour::Black()
            );
            entities.push_back(radar);
            entities.push_back(line);
        }

        for (const auto &[topic, _]: CalibParamManager::EXTRI.SO3_LkToBr) {
            auto SE3_LkToBr = EXTRI.SE3_LkToBr(topic).cast<float>();
            auto lidar = ns_viewer::LiDAR::Create(
                    ns_viewer::Posef(SE3_LkToBr.so3().matrix(), SE3_LkToBr.translation()),
                    LiDAR_SIZE, ns_viewer::Colour(1.0f, 0.5f, 0.0f, 1.0f)
            );
            auto line = ns_viewer::Line::Create(
                    Eigen::Vector3f::Zero(), SE3_LkToBr.translation().cast<float>(), ns_viewer::Colour::Black()
            );
            entities.push_back(lidar);
            entities.push_back(line);
        }

        for (const auto &[topic, _]: CalibParamManager::EXTRI.SO3_CmToBr) {
            auto SE3_CmToBr = EXTRI.SE3_CmToBr(topic).cast<float>();
            auto camera = ns_viewer::CubeCamera::Create(
                    ns_viewer::Posef(SE3_CmToBr.so3().matrix(), SE3_CmToBr.translation()),
                    CAMERA_SIZE, ns_viewer::Colour::Blue()
            );
            auto line = ns_viewer::Line::Create(
                    Eigen::Vector3f::Zero(), SE3_CmToBr.translation().cast<float>(), ns_viewer::Colour::Black()
            );
            entities.push_back(camera);
            entities.push_back(line);
        }

#undef IMU_SIZE
#undef RADAR_SIZE
#undef LiDAR_SIZE
#undef CAMERA_SIZE

        return entities;
    }

    ns_veta::PinholeIntrinsic::Ptr
    CalibParamManager::ParIntri::LoadCameraIntri(const std::string &filename, CerealArchiveType::Enum archiveType) {
        auto intri = ns_veta::PinholeIntrinsic::Create(0, 0, 0, 0, 0, 0);
        std::ifstream file(filename, std::ios::in);
        auto ar = GetInputArchiveVariant(file, archiveType);
        SerializeByInputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", intri));
        return intri;
    }

    IMUIntrinsics::Ptr
    CalibParamManager::ParIntri::LoadIMUIntri(const std::string &filename, CerealArchiveType::Enum archiveType) {
        auto intri = IMUIntrinsics::Create();
        std::ifstream file(filename, std::ios::in);
        auto ar = GetInputArchiveVariant(file, archiveType);
        SerializeByInputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", intri));
        return intri;
    }

    cv::Mat CalibParamManager::ParIntri::UndistortImage(const ns_veta::PinholeIntrinsic::Ptr &intri,
                                                        const cv::Mat &src) {
        auto [K, D] = ObtainKDMatForUndisto(intri);
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
            throw Status(
                    Status::CRITICAL, "unknown camera intrinsic model! supported models:\n"
                                      "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                                      "(b)  pinhole_fisheye (k1, k2, k3, k4)"
            );
        }
        return undistImg;
    }

    std::pair<cv::Mat, cv::Mat>
    CalibParamManager::ParIntri::ObtainKDMatForUndisto(const ns_veta::PinholeIntrinsic::Ptr &intri) {
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
            throw Status(
                    Status::CRITICAL, "unknown camera intrinsic model! supported models:\n"
                                      "(a) pinhole_brown_t2 (k1, k2, k3, p1, p2)\n"
                                      "(b)  pinhole_fisheye (k1, k2, k3, k4)"
            );
        }

        return {K, D};
    }
}
