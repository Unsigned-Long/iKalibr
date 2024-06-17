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

#ifndef IKALIBR_IMU_INTRINSIC_HPP
#define IKALIBR_IMU_INTRINSIC_HPP

#include "util/utils.h"
#include "sensor/imu.h"
#include "util/cereal_archive_helper.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct BiasMapCoeff {
    Eigen::Vector3d BIAS;
    /**
     * MAP_COEFF: [v1, v2, v3, v4, v5, v6]^T
     * mapMatrix:
     *   v1 & v4 & v5
     *    0 & v2 & v6
     *    0 &  0 & v3
     * f(measure) = mapMat * f(real) + bias
     */
    Eigen::Vector6d MAP_COEFF;

    // organize the vector to a matrix
    [[nodiscard]] Eigen::Matrix3d MapMatrix() const {
        Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
        mat(0, 0) = MAP_COEFF(0), mat(1, 1) = MAP_COEFF(1), mat(2, 2) = MAP_COEFF(2);
        mat(0, 1) = MAP_COEFF(3);
        mat(0, 2) = MAP_COEFF(4);
        mat(1, 2) = MAP_COEFF(5);
        return mat;
    }

    BiasMapCoeff() { Clear(); }

    void Clear() {
        BIAS = Eigen::Vector3d::Zero();
        MAP_COEFF = Eigen::Vector6d::Zero();
        MAP_COEFF(0) = 1.0;
        MAP_COEFF(1) = 1.0;
        MAP_COEFF(2) = 1.0;
    }

    // Serialization
    template <class Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(BIAS), CEREAL_NVP(MAP_COEFF));
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct IMUIntrinsics {
    using Ptr = std::shared_ptr<IMUIntrinsics>;

    // trans radian angle to degree angle
    constexpr static double RAD_TO_DEG = 180.0 / M_PI;
    // trans degree angle to radian angle
    constexpr static double DEG_TO_RAD = M_PI / 180.0;

    BiasMapCoeff ACCE, GYRO;
    Sophus::SO3d SO3_AtoG;

    IMUIntrinsics() { Clear(); }

    static Ptr Create() { return std::make_shared<IMUIntrinsics>(); }

    void Clear() {
        ACCE.Clear();
        GYRO.Clear();
        SO3_AtoG = Sophus::SO3d();
    }

    [[nodiscard]] static IMUFrame::Ptr KinematicsToInertialMes(double time,
                                                               const Eigen::Vector3d &linAcceInW,
                                                               const Eigen::Vector3d &angVelInW,
                                                               const Sophus::SO3d &so3CurToW,
                                                               const Eigen::Vector3d &gravityInW) {
        Eigen::Vector3d force = so3CurToW.inverse() * (linAcceInW - gravityInW);
        Eigen::Vector3d angVel = so3CurToW.inverse() * angVelInW;
        return IMUFrame::Create(time, angVel, force);
    }

    [[nodiscard]] static std::tuple<double, Eigen::Vector3d, Eigen::Vector3d>
    InertialMesToKinematics(double time,
                            const IMUFrame::Ptr &frame,
                            const Sophus::SO3d &so3CurToW,
                            const Eigen::Vector3d &gravityInW) {
        Eigen::Vector3d acceInW = so3CurToW * frame->GetAcce() + gravityInW;
        Eigen::Vector3d angVelInW = so3CurToW * frame->GetGyro();
        return {time, acceInW, angVelInW};
    }

    [[nodiscard]] Eigen::Vector3d InvolveForceIntri(const Eigen::Vector3d &force) const {
        return ACCE.MapMatrix() * force + ACCE.BIAS;
    }

    [[nodiscard]] Eigen::Vector3d InvolveGyroIntri(const Eigen::Vector3d &gyro) const {
        return GYRO.MapMatrix() * (SO3_AtoG * gyro) + GYRO.BIAS;
    }

    [[nodiscard]] IMUFrame::Ptr InvolveIntri(const IMUFrame::Ptr &frame) const {
        return IMUFrame::Create(frame->GetTimestamp(), InvolveGyroIntri(frame->GetGyro()),
                                InvolveForceIntri(frame->GetAcce()));
    }

    [[nodiscard]] Eigen::Vector3d RemoveForceIntri(const Eigen::Vector3d &force) const {
        return ACCE.MapMatrix().inverse() * (force - ACCE.BIAS);
    }

    [[nodiscard]] Eigen::Vector3d RemoveGyroIntri(const Eigen::Vector3d &gyro) const {
        return SO3_AtoG.inverse().matrix() * GYRO.MapMatrix().inverse() * (gyro - GYRO.BIAS);
    }

    [[nodiscard]] IMUFrame::Ptr RemoveIntri(const IMUFrame::Ptr &frame) const {
        return IMUFrame::Create(frame->GetTimestamp(), RemoveGyroIntri(frame->GetGyro()),
                                RemoveForceIntri(frame->GetAcce()));
    }

    // quaternion
    [[nodiscard]] Eigen::Quaterniond Q_AtoG() const { return SO3_AtoG.unit_quaternion(); }

    // euler angles
    [[nodiscard]] Eigen::Vector3d EULER_AtoG_RAD() const {
        return Q_AtoG().toRotationMatrix().eulerAngles(0, 1, 2);
    }

    [[nodiscard]] Eigen::Vector3d EULER_AtoG_DEG() const {
        auto euler = EULER_AtoG_RAD();
        for (int i = 0; i != 3; ++i) {
            euler(i) *= IMUIntrinsics::RAD_TO_DEG;
        }
        return euler;
    }

    // save the parameters to file using cereal library
    void Save(const std::string &filename,
              CerealArchiveType::Enum archiveType = CerealArchiveType::Enum::YAML) const {
        std::ofstream file(filename, std::ios::out);
        auto ar = GetOutputArchiveVariant(file, archiveType);
        SerializeByOutputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", *this));
    }

    // load the parameters from file using cereal library
    static IMUIntrinsics::Ptr Load(
        const std::string &filename,
        CerealArchiveType::Enum archiveType = CerealArchiveType::Enum::YAML) {
        auto intri = IMUIntrinsics::Create();
        std::ifstream file(filename, std::ios::in);
        auto ar = GetInputArchiveVariant(file, archiveType);
        SerializeByInputArchiveVariant(ar, archiveType, cereal::make_nvp("Intrinsics", *intri));
        return intri;
    }

    // Serialization
    template <class Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(ACCE), CEREAL_NVP(GYRO), CEREAL_NVP(SO3_AtoG));
    }
};

}  // namespace ns_ikalibr
#endif  // IKALIBR_IMU_INTRINSIC_HPP
