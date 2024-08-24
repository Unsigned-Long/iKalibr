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

#ifndef IKALIBR_IMU_INTRI_CALIB_H
#define IKALIBR_IMU_INTRI_CALIB_H

#include "util/utils.h"
#include "sensor/imu_intrinsic.hpp"
#include "cereal/types/utility.hpp"
#include "util/status.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
class IMUIntriCalibSolver {
public:
    using Ptr = std::shared_ptr<IMUIntriCalibSolver>;

public:
    struct Configor {
    public:
        struct Item {
            std::string bagPath;
            std::vector<std::pair<double, double>> staticPieces;

            Item() = default;

        public:
            template <class Archive>
            void serialize(Archive &ar) {
                ar(cereal::make_nvp("BagPath", bagPath),
                   cereal::make_nvp("StaticPieces", staticPieces));
            }
        };

    public:
        std::vector<Item> items;
        std::string IMUType;
        std::string IMUTopic;
        double gravityNorm{};
        std::string outputPath;

    public:
        Configor() = default;

        // load configure information from file
        template <class CerealArchiveType = CerealArchiveType::YAML>
        static Configor LoadConfigure(const std::string &filename) {
            std::ifstream file(filename);
            auto archive = GetInputArchive<CerealArchiveType>(file);
            auto configor = Configor();
            try {
                (*archive)(cereal::make_nvp("Configor", configor));
            } catch (const cereal::Exception &exception) {
                throw Status(Status::CRITICAL,
                             "The configuration file '{}' for 'IMUIntriCalibSolver::Configor' is "
                             "outdated or broken, and can not be loaded in iKalibr using cereal!!! "
                             "To make it right, please refer to our latest configuration file "
                             "template released at "
                             "https://github.com/Unsigned-Long/iKalibr/blob/master/config/tool/"
                             "config-imu-intri-calib.yaml, and then fix your custom configuration "
                             "file. Detailed cereal "
                             "exception information: \n'{}'",
                             filename, exception.what());
            }
            return configor;
        }

        // save configure information to file
        template <class CerealArchiveType = CerealArchiveType::YAML>
        bool SaveConfigure(const std::string &filename) {
            std::ofstream file(filename);
            auto archive = GetOutputArchive<CerealArchiveType>(file);
            (*archive)(cereal::make_nvp("Configor", *this));
            return true;
        }

    public:
        template <class Archive>
        void serialize(Archive &ar) {
            ar(cereal::make_nvp("IMUTopic", IMUTopic), CEREAL_NVP(IMUType),
               cereal::make_nvp("GravityNorm", gravityNorm),
               cereal::make_nvp("OutputPath", outputPath), cereal::make_nvp("ROSBags", items));
        }
    };

protected:
    Configor configor;
    std::vector<std::list<IMUFrame::Ptr>> data;
    std::vector<Eigen::Vector3d> gravity;
    IMUIntrinsics intrinsics;

public:
    explicit IMUIntriCalibSolver(Configor configor);

    static Ptr Create(const Configor &configor);

    void Process();

    [[nodiscard]] const IMUIntrinsics &GetIntrinsics() const;

protected:
    void LoadIMUData();

    static Eigen::Vector3d AverageAcce(const std::list<IMUFrame::Ptr> &frames);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_IMU_INTRI_CALIB_H
