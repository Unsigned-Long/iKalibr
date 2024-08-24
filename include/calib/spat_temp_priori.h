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

#ifndef IKALIBR_SPAT_TEMP_PRIORI_H
#define IKALIBR_SPAT_TEMP_PRIORI_H

#include "util/utils.h"
#include "util/cereal_archive_helper.hpp"
#include "cereal/types/utility.hpp"
#include "util/status.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct CalibParamManager;
struct Estimator;

class SpatialTemporalPriori {
public:
    using Ptr = std::shared_ptr<SpatialTemporalPriori>;
    using FromTo = std::pair<std::string, std::string>;
    constexpr static double PrioriWeight = 1E4;

public:
    // extrinsics
    std::map<FromTo, Sophus::SO3d> SO3_Sen1ToSen2;
    std::map<FromTo, Eigen::Vector3d> POS_Sen1InSen2;
    // time offsets
    std::map<FromTo, double> TO_Sen1ToSen2;
    // readout time of rs cameras
    std::map<std::string, double> RS_READOUT;

public:
    SpatialTemporalPriori() = default;

    static Ptr Create();

    [[nodiscard]] const std::map<FromTo, Sophus::SO3d> &GetExtriSO3() const;

    [[nodiscard]] const std::map<FromTo, Eigen::Vector3d> &GetExtriPOS() const;

    [[nodiscard]] const std::map<FromTo, double> &GetTimeOffset() const;

    [[nodiscard]] const std::map<std::string, double> &GetReadout() const;

    void CheckValidityWithConfigor() const;

    void AddSpatTempPrioriConstraint(Estimator &estimator, CalibParamManager &parMagr) const;

protected:
    template <class ValueType>
    static std::pair<bool, std::pair<std::string, std::string>> IsMapAmbiguous(
        const std::map<std::pair<std::string, std::string>, ValueType> &myMap) {
        for (const auto &entry : myMap) {
            std::pair<std::string, std::string> reversedKey = {entry.first.second,
                                                               entry.first.first};
            if (myMap.find(reversedKey) != myMap.end()) {
                return {true, entry.first};
            }
        }
        return {false, {}};
    }

public:
    // Serialization
    template <class Archive>
    void serialize(Archive &ar) {
        ar(CEREAL_NVP(SO3_Sen1ToSen2), CEREAL_NVP(POS_Sen1InSen2), CEREAL_NVP(TO_Sen1ToSen2),
           CEREAL_NVP(RS_READOUT));
    }

    // save the parameters to file using cereal library
    template <class CerealArchiveType = CerealArchiveType::YAML>
    void Save(const std::string &filename) const {
        std::ofstream file(filename, std::ios::out);
        auto ar = GetOutputArchive<CerealArchiveType>(file);

        (*ar)(cereal::make_nvp("SpatialTemporalPriori", *this));
    }

    // load the parameters from file using cereal library
    template <class CerealArchiveType = CerealArchiveType::YAML>
    static SpatialTemporalPriori::Ptr Load(const std::string &filename) {
        auto priori = SpatialTemporalPriori::Create();
        std::ifstream file(filename, std::ios::in);
        auto ar = GetInputArchive<CerealArchiveType>(file);
        try {
            (*ar)(cereal::make_nvp("SpatialTemporalPriori", *priori));
        } catch (const cereal::Exception &exception) {
            throw Status(Status::CRITICAL,
                         "The configuration file '{}' for 'SpatialTemporalPriori' is "
                         "outdated or broken, and can not be loaded in iKalibr using cereal!!! "
                         "To make it right, please refer to our latest configuration file "
                         "template released at "
                         "https://github.com/Unsigned-Long/iKalibr/blob/master/config/"
                         "spat-temp-priori.yaml, and then fix your custom configuration "
                         "file. Detailed cereal "
                         "exception information: \n'{}'",
                         filename, exception.what());
        }
        return priori;
    }
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_SPAT_TEMP_PRIORI_H
