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

#include "solver/calib_solver.h"
#include "util/utils_tpl.hpp"
#include "spdlog/spdlog.h"
#include "calib/estimator.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void CalibSolver::InitSO3Spline() const {
    /**
     * this function would initialize the rotation spline, as well as the extrinsic rotations and
     * time offsets between multiple imus, if they are integrated
     */
    spdlog::info("fitting rotation b-spline...");

    /**
     * here we recover the so3 spline first using only the angular velocities from the reference
     * IMU, then estimates other quantities by involving angular velocity measurements from other
     * IMUs. For better readability, we could also optimize them together as follows.
     * ----------------------------------------------------------------
     * auto estimator = Estimator::Create(_splines, _parMagr);
     * auto optOption = OptOption::OPT_SO3_BiToBr;
     * if (Configor::Prior::OptTemporalParams) {
     *     optOption |= OptOption::OPT_TO_BiToBr;
     * }
     * for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
     *     this->AddGyroFactor(estimator, topic, optOption);
     * }
     * estimator->SetRefIMUParamsConstant();
     * auto sum = estimator->Solve(_ceresOption, this->_priori);
     * spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
     * ----------------------------------------------------------------
     */

    auto estimator = Estimator::Create(_splines, _parMagr);
    // we initialize the rotation spline first use only the measurements from the reference imu
    this->AddGyroFactor(estimator, Configor::DataStream::ReferIMU, OptOption::OPT_SO3_SPLINE);
    auto sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

    if (Configor::DataStream::IMUTopics.size() > 1) {
        // if multiple imus involved, we continue to recover extrinsic rotations and time offsets
        estimator = Estimator::Create(_splines, _parMagr);
        auto optOption = OptOption::OPT_SO3_BiToBr;
        if (Configor::Prior::OptTemporalParams) {
            optOption |= OptOption::OPT_TO_BiToBr;
        }
        for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
            this->AddGyroFactor(estimator, topic, optOption);
        }
        // make this problem full rank
        estimator->SetRefIMUParamsConstant();

        sum = estimator->Solve(_ceresOption, this->_priori);
        spdlog::info("here is the summary:\n{}\n", sum.BriefReport());
    }
}

}  // namespace ns_ikalibr