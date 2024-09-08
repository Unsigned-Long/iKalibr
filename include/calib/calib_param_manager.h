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

#ifndef IKALIBR_CALIB_PARAM_MANAGER_H
#define IKALIBR_CALIB_PARAM_MANAGER_H

#include "config/configor.h"
#include "sensor/imu_intrinsic.hpp"
#include "sensor/rgbd_intrinsic.hpp"
#include "tiny-viewer/core/viewer.h"
#include "tiny-viewer/core/multi_viewer.h"
#include "veta/camera/pinhole_brown.h"
#include "veta/camera/pinhole_fisheye.h"
#include "opencv2/core.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
#define SE3_SEN_TO_REF(SENSOR1, IDX1, SENSOR2, IDX2)                                            \
    [[nodiscard]] Sophus::SE3d SE3_##SENSOR1##IDX1##To##SENSOR2##IDX2(const std::string &topic) \
        const {                                                                                 \
        return {SO3_##SENSOR1##IDX1##To##SENSOR2##IDX2.at(topic),                               \
                POS_##SENSOR1##IDX1##In##SENSOR2##IDX2.at(topic)};                              \
    }

#define Q_SEN_TO_REF(SENSOR1, IDX1, SENSOR2, IDX2)                                 \
    [[nodiscard]] Eigen::Quaterniond Q_##SENSOR1##IDX1##To##SENSOR2##IDX2(         \
        const std::string &topic) const {                                          \
        return SO3_##SENSOR1##IDX1##To##SENSOR2##IDX2.at(topic).unit_quaternion(); \
    }

#define EULER_SEN_TO_REF_RAD(SENSOR1, IDX1, SENSOR2, IDX2)                                      \
    [[nodiscard]] Eigen::Vector3d EULER_##SENSOR1##IDX1##To##SENSOR2##IDX2##_RAD(               \
        const std::string &topic) const {                                                       \
        return Q_##SENSOR1##IDX1##To##SENSOR2##IDX2(topic).toRotationMatrix().eulerAngles(2, 1, \
                                                                                          0);   \
    }

#define EULER_SEN_TO_REF_DEG(SENSOR1, IDX1, SENSOR2, IDX2)                        \
    [[nodiscard]] Eigen::Vector3d EULER_##SENSOR1##IDX1##To##SENSOR2##IDX2##_DEG( \
        const std::string &topic) const {                                         \
        auto euler = EULER_##SENSOR1##IDX1##To##SENSOR2##IDX2##_RAD(topic);       \
        for (int i = 0; i != 3; ++i) {                                            \
            euler(i) *= CalibParamManager::RAD_TO_DEG;                            \
        }                                                                         \
        return euler;                                                             \
    }

struct CalibParamManager {
public:
    using Ptr = std::shared_ptr<CalibParamManager>;

public:
    // trans radian angle to degree angle
    constexpr static double RAD_TO_DEG = 180.0 / M_PI;

    static struct Header {
        const static std::string Software;
        const static std::string Version;
        const static std::string Address;

    public:
        template <class Archive>
        void serialize(Archive &ar) {
            ar(CEREAL_NVP(Software), CEREAL_NVP(Version), CEREAL_NVP(Address));
        }
    } header;

    // ---------
    // extrinsic
    // ---------
    struct ParExtri {
        // topic, SO3
        std::map<std::string, Sophus::SO3d> SO3_BiToBr;
        std::map<std::string, Sophus::SO3d> SO3_RjToBr;
        std::map<std::string, Sophus::SO3d> SO3_LkToBr;
        std::map<std::string, Sophus::SO3d> SO3_CmToBr;
        std::map<std::string, Sophus::SO3d> SO3_DnToBr;

        // topic, POS
        std::map<std::string, Eigen::Vector3d> POS_BiInBr;
        std::map<std::string, Eigen::Vector3d> POS_RjInBr;
        std::map<std::string, Eigen::Vector3d> POS_LkInBr;
        std::map<std::string, Eigen::Vector3d> POS_CmInBr;
        std::map<std::string, Eigen::Vector3d> POS_DnInBr;

        SE3_SEN_TO_REF(B, i, B, r)

        SE3_SEN_TO_REF(R, j, B, r)

        SE3_SEN_TO_REF(L, k, B, r)

        SE3_SEN_TO_REF(C, m, B, r)

        SE3_SEN_TO_REF(D, n, B, r)

        Q_SEN_TO_REF(B, i, B, r)

        Q_SEN_TO_REF(R, j, B, r)

        Q_SEN_TO_REF(L, k, B, r)

        Q_SEN_TO_REF(C, m, B, r)

        Q_SEN_TO_REF(D, n, B, r)

        EULER_SEN_TO_REF_RAD(B, i, B, r)

        EULER_SEN_TO_REF_RAD(R, j, B, r)

        EULER_SEN_TO_REF_RAD(L, k, B, r)

        EULER_SEN_TO_REF_RAD(C, m, B, r)

        EULER_SEN_TO_REF_RAD(D, n, B, r)

        EULER_SEN_TO_REF_DEG(B, i, B, r)

        EULER_SEN_TO_REF_DEG(R, j, B, r)

        EULER_SEN_TO_REF_DEG(L, k, B, r)

        EULER_SEN_TO_REF_DEG(C, m, B, r)

        EULER_SEN_TO_REF_DEG(D, n, B, r)

    public:
        // Serialization
        template <class Archive>
        void serialize(Archive &archive) {
            archive(CEREAL_NVP(SO3_BiToBr), CEREAL_NVP(POS_BiInBr), CEREAL_NVP(SO3_RjToBr),
                    CEREAL_NVP(POS_RjInBr), CEREAL_NVP(SO3_LkToBr), CEREAL_NVP(POS_LkInBr),
                    CEREAL_NVP(SO3_CmToBr), CEREAL_NVP(POS_CmInBr), CEREAL_NVP(SO3_DnToBr),
                    CEREAL_NVP(POS_DnInBr));
        }

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    } EXTRI;

    // --------
    // temporal
    // --------
    struct ParTemporal {
        // topic, time offset
        std::map<std::string, double> TO_BiToBr;
        std::map<std::string, double> TO_RjToBr;
        std::map<std::string, double> TO_LkToBr;
        std::map<std::string, double> TO_CmToBr;
        std::map<std::string, double> TO_DnToBr;
        std::map<std::string, double> RS_READOUT;

    public:
        // Serialization
        template <class Archive>
        void serialize(Archive &ar) {
            ar(CEREAL_NVP(TO_BiToBr), CEREAL_NVP(TO_RjToBr), CEREAL_NVP(TO_LkToBr),
               CEREAL_NVP(TO_CmToBr), CEREAL_NVP(TO_DnToBr), CEREAL_NVP(RS_READOUT));
        }
    } TEMPORAL;

    // ---------
    // intrinsic
    // ---------
    struct ParIntri {
        // topic, param pack
        std::map<std::string, IMUIntrinsics::Ptr> IMU;
        std::map<std::string, ns_veta::PinholeIntrinsic::Ptr> Camera;
        std::map<std::string, RGBDIntrinsics::Ptr> RGBD;

        static ns_veta::PinholeIntrinsic::Ptr LoadCameraIntri(
            const std::string &filename,
            CerealArchiveType::Enum archiveType = CerealArchiveType::Enum::YAML);

        static IMUIntrinsics::Ptr LoadIMUIntri(
            const std::string &filename,
            CerealArchiveType::Enum archiveType = CerealArchiveType::Enum::YAML);

        static void SaveCameraIntri(
            const ns_veta::PinholeIntrinsic::Ptr &intri,
            const std::string &filename,
            CerealArchiveType::Enum archiveType = CerealArchiveType::Enum::YAML);

        static void SaveIMUIntri(
            const IMUIntrinsics::Ptr &intri,
            const std::string &filename,
            CerealArchiveType::Enum archiveType = CerealArchiveType::Enum::YAML);

        static cv::Mat UndistortImage(const ns_veta::PinholeIntrinsic::Ptr &intri,
                                      const cv::Mat &src);

    private:
        [[nodiscard]] static std::pair<cv::Mat, cv::Mat> ObtainKDMatForUndisto(
            const ns_veta::PinholeIntrinsic::Ptr &intri);

    public:
        // Serialization
        template <class Archive>
        void serialize(Archive &archive) {
            archive(CEREAL_NVP(IMU), CEREAL_NVP(Camera), CEREAL_NVP(RGBD));
        }

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    } INTRI;

    // S2Manifold
    Eigen::Vector3d GRAVITY;

public:
    // the constructor
    explicit CalibParamManager(const std::vector<std::string> &imuTopics = {},
                               const std::vector<std::string> &radarTopics = {},
                               const std::vector<std::string> &lidarTopics = {},
                               const std::vector<std::string> &cameraTopics = {},
                               const std::vector<std::string> &rgbdTopics = {});

    // the creator
    static CalibParamManager::Ptr Create(const std::vector<std::string> &imuTopics = {},
                                         const std::vector<std::string> &radarTopics = {},
                                         const std::vector<std::string> &lidarTopics = {},
                                         const std::vector<std::string> &cameraTopics = {},
                                         const std::vector<std::string> &rgbdTopics = {});

    // save the parameters to file using cereal library
    template <class CerealArchiveType = CerealArchiveType::YAML>
    void Save(const std::string &filename) const {
        std::ofstream file(filename, std::ios::out);
        auto ar = GetOutputArchive<CerealArchiveType>(file);

        (*ar)(cereal::make_nvp("CalibParam", *this));
    }

    // load the parameters from file using cereal library
    template <class CerealArchiveType = CerealArchiveType::YAML>
    static CalibParamManager::Ptr Load(const std::string &filename) {
        auto calibParamManager = CalibParamManager::Create();
        std::ifstream file(filename, std::ios::in);
        auto ar = GetInputArchive<CerealArchiveType>(file);

        (*ar)(cereal::make_nvp("CalibParam", *calibParamManager));
        return calibParamManager;
    }

    // save the parameters to file using cereal library
    void Save(const std::string &filename, CerealArchiveType::Enum archiveType) const;

    // load the parameters from file using cereal library
    static CalibParamManager::Ptr Load(const std::string &filename,
                                       CerealArchiveType::Enum archiveType);

    // print the parameters in the console
    void ShowParamStatus();

    std::vector<std::size_t> VisualizationSensors(ns_viewer::Viewer &viewer) const;

    std::vector<std::size_t> VisualizationSensors(ns_viewer::MultiViewer &viewer,
                                                  const std::string &win) const;

    // set the params to the init values, the intrinsic coeff of camera will load from the config
    // file make sure load and check config before initialize the parameters
    static CalibParamManager::Ptr InitParamsFromConfigor();

public:
    // Serialization
    template <class Archive>
    void save(Archive &archive) const {
        archive(cereal::make_nvp("Header", header), CEREAL_NVP(EXTRI), CEREAL_NVP(TEMPORAL),
                CEREAL_NVP(INTRI), CEREAL_NVP(GRAVITY));
    }

    template <class Archive>
    void load(Archive &archive) {
        archive(CEREAL_NVP(EXTRI), CEREAL_NVP(TEMPORAL), CEREAL_NVP(INTRI), CEREAL_NVP(GRAVITY));
    }

private:
    [[nodiscard]] std::vector<ns_viewer::EntityPtr> EntitiesForVisualization() const;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#undef SE3_SEN_TO_REF
#undef Q_SEN_TO_REF
#undef EULER_SEN_TO_REF_RAD
#undef EULER_SEN_TO_REF_DEG
}  // namespace ns_ikalibr

#endif  // IKALIBR_CALIB_PARAM_MANAGER_H
