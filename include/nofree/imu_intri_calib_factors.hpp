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

#ifndef IKALIBR_IMU_INTRI_CALIB_FACTORS_HPP
#define IKALIBR_IMU_INTRI_CALIB_FACTORS_HPP

#include "ctraj/utils/eigen_utils.hpp"
#include "ctraj/utils/sophus_utils.hpp"
#include "ctraj/spline/spline_segment.h"
#include "ctraj/spline/ceres_spline_helper.h"
#include "ctraj/spline/ceres_spline_helper_jet.h"
#include "ceres/ceres.h"
#include "sensor/imu.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct IMUIntriAcceFactor {
private:
    IMUFrame::Ptr _imuFrame{};
    double _weight;

public:
    explicit IMUIntriAcceFactor(IMUFrame::Ptr imuFrame, double weight)
        : _imuFrame(std::move(imuFrame)),
          _weight(weight) {}

    static auto Create(const IMUFrame::Ptr &imuFrame, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<IMUIntriAcceFactor>(
            new IMUIntriAcceFactor(imuFrame, weight));
    }

    static std::size_t TypeHashCode() { return typeid(IMUIntriAcceFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ ACCE_BIAS | ACCE_MAP_COEFF | GRAVITY ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        std::size_t ACCE_BIAS_OFFSET = 0;
        std::size_t ACCE_MAP_COEFF_OFFSET = ACCE_BIAS_OFFSET + 1;
        std::size_t GRAVITY_OFFSET = ACCE_MAP_COEFF_OFFSET + 1;

        Eigen::Map<const Eigen::Vector3<T>> acceBias(sKnots[ACCE_BIAS_OFFSET]);
        Eigen::Map<const Eigen::Vector3<T>> gravity(sKnots[GRAVITY_OFFSET]);

        auto acceCoeff = sKnots[ACCE_MAP_COEFF_OFFSET];

        Eigen::Matrix33<T> acceMapMat = Eigen::Matrix33<T>::Zero();

        acceMapMat.diagonal() = Eigen::Map<const Eigen::Vector3<T>>(acceCoeff, 3);
        acceMapMat(0, 1) = *(acceCoeff + 3);
        acceMapMat(0, 2) = *(acceCoeff + 4);
        acceMapMat(1, 2) = *(acceCoeff + 5);

        Eigen::Vector3<T> accePred = -acceMapMat * gravity + acceBias;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = accePred - _imuFrame->GetAcce().template cast<T>();
        residuals = T(_weight) * residuals;

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct IMUIntriGyroFactor {
private:
    IMUFrame::Ptr _frame{};
    double _weight;

public:
    explicit IMUIntriGyroFactor(IMUFrame::Ptr frame, double weight)
        : _frame(std::move(frame)),
          _weight(weight) {}

    static auto Create(const IMUFrame::Ptr &frame, double weight) {
        return new ceres::DynamicAutoDiffCostFunction<IMUIntriGyroFactor>(
            new IMUIntriGyroFactor(frame, weight));
    }

    static std::size_t TypeHashCode() { return typeid(IMUIntriGyroFactor).hash_code(); }

public:
    /**
     * param blocks:
     * [ GYRO_BIAS | GYRO_MAP_COEFF | SO3_AtoG ]
     */
    template <class T>
    bool operator()(T const *const *sKnots, T *sResiduals) const {
        // array offset
        std::size_t GYRO_BIAS_OFFSET = 0;
        std::size_t GYRO_MAP_COEFF_OFFSET = GYRO_BIAS_OFFSET + 1;
        std::size_t SO3_AtoG_OFFSET = GYRO_MAP_COEFF_OFFSET + 1;

        Eigen::Map<const Eigen::Vector3<T>> gyroBias(sKnots[GYRO_BIAS_OFFSET]);
        // thi value would be fixed as identities
        auto gyroCoeff = sKnots[GYRO_MAP_COEFF_OFFSET];
        Eigen::Matrix33<T> gyroMapMat = Eigen::Matrix33<T>::Zero();
        gyroMapMat.diagonal() = Eigen::Map<const Eigen::Vector3<T>>(gyroCoeff, 3);
        gyroMapMat(0, 1) = *(gyroCoeff + 3);
        gyroMapMat(0, 2) = *(gyroCoeff + 4);
        gyroMapMat(1, 2) = *(gyroCoeff + 5);

        // thi value would be fixed as identities
        Eigen::Map<Sophus::SO3<T> const> const SO3_AtoG(sKnots[SO3_AtoG_OFFSET]);

        // we do not consider the earth rotation
        Eigen::Vector3<T> gyroPred = gyroBias;

        Eigen::Map<Eigen::Vector3<T>> residuals(sResiduals);
        residuals = gyroPred - _frame->GetGyro().template cast<T>();
        residuals = T(_weight) * residuals;

        return true;
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace ns_ikalibr
#endif  // IKALIBR_IMU_INTRI_CALIB_FACTORS_HPP
