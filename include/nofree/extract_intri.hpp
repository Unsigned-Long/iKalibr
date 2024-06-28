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

#ifndef IKALIBR_EXTRACT_INTRI_HPP
#define IKALIBR_EXTRACT_INTRI_HPP

#include <utility>

#include "calib/calib_param_manager.h"

namespace ns_ikalibr {
struct ExtractIntriParam {
protected:
    CalibParamManager::Ptr _calibPar;

public:
    explicit ExtractIntriParam(CalibParamManager::Ptr calibPar)
        : _calibPar(std::move(calibPar)) {}

    void ToDisk(const std::string& dir,
                bool cam = true,
                bool imu = true,
                CerealArchiveType::Enum type = CerealArchiveType::Enum::YAML) {
        // cameras
        if (cam) {
            for (const auto& [topic, intri] : _calibPar->INTRI.Camera) {
                const std::string path = GetSavePath(topic, dir, type);
                CalibParamManager::ParIntri::SaveCameraIntri(intri, path, type);
            }
        }
        // imus
        if (imu) {
            for (const auto& [topic, intri] : _calibPar->INTRI.IMU) {
                const std::string path = GetSavePath(topic, dir, type);
                CalibParamManager::ParIntri::SaveIMUIntri(intri, path, type);
            }
        }
    }

    static void ToDisk(const std::vector<std::string>& wsVec,
                       bool cam = true,
                       bool imu = true,
                       CerealArchiveType::Enum type = CerealArchiveType::Enum::YAML) {
        for (const auto& ws : wsVec) {
            spdlog::info("extract intrinsics in workspace '{}'", ws);
            ExtractIntriParam(ns_ikalibr::CalibParamManager::Load(ws + "/ikalibr_param.yaml"))
                .ToDisk(ws, cam, imu, type);
        }
    }

protected:
    static std::string GetSavePath(const std::string& topic,
                                   const std::string& dir,
                                   CerealArchiveType::Enum type) {
        std::string filename = topic;
        std::replace(std::next(filename.begin()), filename.end(), '/', '-');
        std::string path =
            dir + '/' + filename + "-intri" + Configor::Preference::FileExtension.at(type);
        return path;
    }
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_EXTRACT_INTRI_HPP
