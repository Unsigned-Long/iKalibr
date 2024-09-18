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

#ifndef IKALIBR_BATCH_OPT_OPTION_HPP
#define IKALIBR_BATCH_OPT_OPTION_HPP

#include "calib/estimator.h"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct BatchOptOption {
public:
    using Opt = OptOption;

protected:
    constexpr static std::array<Opt, 1> MultiIMU = {Opt::OPT_SO3_SPLINE | Opt::OPT_SCALE_SPLINE |
                                                    Opt::OPT_GRAVITY | Opt::OPT_SO3_BiToBr |
                                                    Opt::OPT_POS_BiInBr | Opt::OPT_TO_BiToBr};

    constexpr static std::array<Opt, 2> MultiRadarIMU = {
        // first batch optimization
        Opt::OPT_SO3_SPLINE | Opt::OPT_SCALE_SPLINE | Opt::OPT_GRAVITY | Opt::OPT_SO3_RjToBr |
            Opt::OPT_POS_RjInBr | Opt::OPT_TO_RjToBr,
        // second batch optimization (append to last)
        Opt::OPT_SO3_BiToBr | Opt::OPT_POS_BiInBr | Opt::OPT_TO_BiToBr | Opt::OPT_ACCE_BIAS |
            Opt::OPT_GYRO_BIAS};

    constexpr static std::array<Opt, 3> MultiLiDARIMU = {
        // first batch optimization
        Opt::OPT_SO3_SPLINE | Opt::OPT_SCALE_SPLINE | Opt::OPT_GRAVITY | Opt::OPT_SO3_LkToBr |
            Opt::OPT_POS_LkInBr | Opt::OPT_TO_LkToBr,
        // second batch optimization (append to last)
        Opt::OPT_SO3_BiToBr | Opt::OPT_POS_BiInBr | Opt::OPT_TO_BiToBr | Opt::OPT_ACCE_BIAS |
            Opt::OPT_GYRO_BIAS,
        // forth batch optimization (append to last)
        Opt::NONE,
    };

    constexpr static std::array<Opt, 2> MultiPosCameraIMU = {
        // first batch optimization
        Opt::OPT_SO3_SPLINE | Opt::OPT_SCALE_SPLINE | Opt::OPT_GRAVITY | Opt::OPT_SO3_CmToBr |
            Opt::OPT_POS_CmInBr | Opt::OPT_TO_CmToBr | Opt::OPT_VISUAL_DEPTH,
        // second batch optimization (append to last)
        Opt::OPT_SO3_BiToBr | Opt::OPT_POS_BiInBr | Opt::OPT_TO_BiToBr |
            Opt::OPT_RS_CAM_READOUT_TIME | Opt::OPT_CAM_FOCAL_LEN | Opt::OPT_CAM_PRINCIPAL_POINT |
            Opt::OPT_ACCE_BIAS | Opt::OPT_GYRO_BIAS};

    constexpr static std::array<Opt, 2> MultiVelCameraIMU = {
        // first batch optimization
        Opt::OPT_SO3_SPLINE | Opt::OPT_SCALE_SPLINE | Opt::OPT_GRAVITY | Opt::OPT_SO3_CmToBr |
            Opt::OPT_POS_CmInBr | Opt::OPT_TO_CmToBr |
            Opt::OPT_VISUAL_DEPTH,  // we always estimate the depth for vel-derived optical cameras
        // second batch optimization (append to last)
        Opt::OPT_SO3_BiToBr | Opt::OPT_POS_BiInBr | Opt::OPT_TO_BiToBr |
            Opt::OPT_RS_CAM_READOUT_TIME | Opt::OPT_CAM_FOCAL_LEN | Opt::OPT_CAM_PRINCIPAL_POINT |
            Opt::OPT_ACCE_BIAS | Opt::OPT_GYRO_BIAS};

    constexpr static std::array<Opt, 2> MultiRGBDIMU = {
        // first batch optimization
        Opt::OPT_SO3_SPLINE | Opt::OPT_SCALE_SPLINE | Opt::OPT_GRAVITY | Opt::OPT_SO3_DnToBr |
            Opt::OPT_POS_DnInBr | Opt::OPT_TO_DnToBr | Opt::OPT_VISUAL_DEPTH,
        // second batch optimization (append to last)
        Opt::OPT_SO3_BiToBr | Opt::OPT_POS_BiInBr | Opt::OPT_TO_BiToBr |
            Opt::OPT_RS_CAM_READOUT_TIME | Opt::OPT_CAM_FOCAL_LEN | Opt::OPT_CAM_PRINCIPAL_POINT |
            Opt::OPT_GYRO_BIAS | Opt::OPT_ACCE_BIAS};

public:
    static std::vector<Opt> GetOptions() {
        std::vector<Opt> options;
        if (!Configor::IsLiDARIntegrated() && !Configor::IsRadarIntegrated() &&
            !Configor::IsPosCameraIntegrated() && !Configor::IsVelCameraIntegrated() &&
            !Configor::IsRGBDIntegrated()) {
            // imu-only multi-imu calibration
            options = AryToVecWithAppend(MultiIMU);
        } else {
            if (Configor::IsLiDARIntegrated()) {
                options = MergeOptions(options, AryToVecWithAppend(MultiLiDARIMU));
            }
            if (Configor::IsPosCameraIntegrated()) {
                options = MergeOptions(options, AryToVecWithAppend(MultiPosCameraIMU));
            }
            if (Configor::IsVelCameraIntegrated()) {
                options = MergeOptions(options, AryToVecWithAppend(MultiVelCameraIMU));
            }
            if (Configor::IsRadarIntegrated()) {
                options = MergeOptions(options, AryToVecWithAppend(MultiRadarIMU));
            }
            if (Configor::IsRGBDIntegrated()) {
                options = MergeOptions(options, AryToVecWithAppend(MultiRGBDIMU));
            }
        }
        if (options.empty()) {
            throw Status(Status::CRITICAL, "unknown error happened! (unknown sensor suite)");
        }
        // if do not optimize temporal parameters, remove them
        if (!Configor::Prior::OptTemporalParams) {
            for (auto &opt : options) {
                RemoveOption(opt, Opt::OPT_TO_BiToBr);
                RemoveOption(opt, Opt::OPT_TO_RjToBr);
                RemoveOption(opt, Opt::OPT_TO_LkToBr);
                RemoveOption(opt, Opt::OPT_TO_CmToBr);
                RemoveOption(opt, Opt::OPT_TO_DnToBr);
                // we do not remove this optimization option, as the RS effect exists even if
                // sensors are hardware-synchronized. in other words, we still optimize this
                // parameter (readout time of RS camera) RemoveOption(opt,
                // Opt::OPT_RS_CAM_READOUT_TIME);
            }
        }
        return options;
    }

protected:
    template <std::size_t Dime>
    static std::vector<Opt> AryToVecWithAppend(const std::array<Opt, Dime> &ary) {
        std::vector<Opt> vec(ary.size(), Opt::NONE);
        for (int i = 0; i < static_cast<int>(ary.size()); ++i) {
            vec.at(i) = ary.at(i);
            // append
            if (i != 0) {
                vec.at(i) |= vec.at(i - 1);
            }
        }
        return vec;
    }

    static std::vector<Opt> MergeOptions(const std::vector<Opt> &optVec1,
                                         const std::vector<Opt> &optVec2) {
        if (optVec1.empty()) {
            return optVec2;
        }
        if (optVec2.empty()) {
            return optVec1;
        }
        std::vector<Opt> vec(std::max(optVec1.size(), optVec2.size()));
        int idx1 = 0, idx2 = 0;
        for (auto &item : vec) {
            Opt op1 = optVec1.at(idx1), op2 = optVec2.at(idx2);
            item = op1 | op2;
            if (idx1 < static_cast<int>(optVec1.size()) - 1) {
                ++idx1;
            }
            if (idx2 < static_cast<int>(optVec2.size()) - 1) {
                ++idx2;
            }
        }
        return vec;
    }

    static void RemoveOption(Opt &srcOpt, Opt optToRemove) {
        if (IsOptionWith(optToRemove, srcOpt)) {
            srcOpt ^= optToRemove;
        }
    }
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_BATCH_OPT_OPTION_HPP
