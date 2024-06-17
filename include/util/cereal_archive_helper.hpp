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

#ifndef IKALIBR_CEREAL_ARCHIVE_HELPER_HPP
#define IKALIBR_CEREAL_ARCHIVE_HELPER_HPP

#include "util/cereal_yaml.hpp"
#include "variant"
#include "cereal/archives/json.hpp"
#include "cereal/archives/xml.hpp"
#include "cereal/cereal.hpp"
#include "cereal/archives/binary.hpp"
#include "memory"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct CerealArchiveType {
    enum class Enum { JSON, XML, YAML, BINARY };
    using InputArchiveVariant = std::variant<std::shared_ptr<cereal::JSONInputArchive>,
                                             std::shared_ptr<cereal::XMLInputArchive>,
                                             std::shared_ptr<cereal::YAMLInputArchive>,
                                             std::shared_ptr<cereal::BinaryInputArchive>>;

    using OutputArchiveVariant = std::variant<std::shared_ptr<cereal::JSONOutputArchive>,
                                              std::shared_ptr<cereal::XMLOutputArchive>,
                                              std::shared_ptr<cereal::YAMLOutputArchive>,
                                              std::shared_ptr<cereal::BinaryOutputArchive>>;

    struct JSON {};
    struct XML {};
    struct YAML {};
    struct BINARY {};
};

template <class Type>
struct CerealArchiveTypeExtractor {
    using InputArchive = CerealArchiveType;
    using OutputArchive = CerealArchiveType;
};

template <>
struct CerealArchiveTypeExtractor<CerealArchiveType::JSON> {
    using InputArchive = cereal::JSONInputArchive;
    using OutputArchive = cereal::JSONOutputArchive;
};

template <>
struct CerealArchiveTypeExtractor<CerealArchiveType::XML> {
    using InputArchive = cereal::XMLInputArchive;
    using OutputArchive = cereal::XMLOutputArchive;
};

template <>
struct CerealArchiveTypeExtractor<CerealArchiveType::YAML> {
    using InputArchive = cereal::YAMLInputArchive;
    using OutputArchive = cereal::YAMLOutputArchive;
};

template <>
struct CerealArchiveTypeExtractor<CerealArchiveType::BINARY> {
    using InputArchive = cereal::BinaryInputArchive;
    using OutputArchive = cereal::BinaryOutputArchive;
};

template <class ArchiveType>
static inline std::shared_ptr<typename CerealArchiveTypeExtractor<ArchiveType>::InputArchive>
GetInputArchive(std::ifstream &file) {
    return std::make_shared<typename CerealArchiveTypeExtractor<ArchiveType>::InputArchive>(file);
}

template <class ArchiveType>
static inline std::shared_ptr<typename CerealArchiveTypeExtractor<ArchiveType>::OutputArchive>
GetOutputArchive(std::ofstream &file) {
    return std::make_shared<typename CerealArchiveTypeExtractor<ArchiveType>::OutputArchive>(file);
}

static inline CerealArchiveType::InputArchiveVariant GetInputArchiveVariant(
    std::ifstream &file, CerealArchiveType::Enum archiveType) {
    switch (archiveType) {
        case CerealArchiveType::Enum::JSON:
            return std::make_shared<cereal::JSONInputArchive>(file);
        case CerealArchiveType::Enum::XML:
            return std::make_shared<cereal::XMLInputArchive>(file);
        case CerealArchiveType::Enum::YAML:
            return std::make_shared<cereal::YAMLInputArchive>(file);
        case CerealArchiveType::Enum::BINARY:
            return std::make_shared<cereal::BinaryInputArchive>(file);
        default:
            throw std::runtime_error("unknown cereal archive type!");
    }
}

static inline CerealArchiveType::OutputArchiveVariant GetOutputArchiveVariant(
    std::ofstream &file, CerealArchiveType::Enum archiveType) {
    switch (archiveType) {
        case CerealArchiveType::Enum::JSON:
            return std::make_shared<cereal::JSONOutputArchive>(file);
        case CerealArchiveType::Enum::XML:
            return std::make_shared<cereal::XMLOutputArchive>(file);
        case CerealArchiveType::Enum::YAML:
            return std::make_shared<cereal::YAMLOutputArchive>(file);
        case CerealArchiveType::Enum::BINARY:
            return std::make_shared<cereal::BinaryOutputArchive>(file);
        default:
            throw std::runtime_error("unknown cereal archive type!");
    }
}

template <class... Types>
static inline void SerializeByInputArchiveVariant(const CerealArchiveType::InputArchiveVariant &ar,
                                                  CerealArchiveType::Enum archiveType,
                                                  Types &&...args) {
    switch (archiveType) {
        case CerealArchiveType::Enum::JSON:
            (*std::get<std::shared_ptr<cereal::JSONInputArchive>>(ar))(args...);
            break;
        case CerealArchiveType::Enum::XML:
            (*std::get<std::shared_ptr<cereal::XMLInputArchive>>(ar))(args...);
            break;
        case CerealArchiveType::Enum::YAML:
            (*std::get<std::shared_ptr<cereal::YAMLInputArchive>>(ar))(args...);
            break;
        case CerealArchiveType::Enum::BINARY:
            (*std::get<std::shared_ptr<cereal::BinaryInputArchive>>(ar))(args...);
            break;
        default:
            throw std::runtime_error("unknown cereal archive type!");
    }
}

template <class... Types>
static inline void SerializeByOutputArchiveVariant(
    const CerealArchiveType::OutputArchiveVariant &ar,
    CerealArchiveType::Enum archiveType,
    Types &&...args) {
    switch (archiveType) {
        case CerealArchiveType::Enum::JSON:
            (*std::get<std::shared_ptr<cereal::JSONOutputArchive>>(ar))(args...);
            break;
        case CerealArchiveType::Enum::XML:
            (*std::get<std::shared_ptr<cereal::XMLOutputArchive>>(ar))(args...);
            break;
        case CerealArchiveType::Enum::YAML:
            (*std::get<std::shared_ptr<cereal::YAMLOutputArchive>>(ar))(args...);
            break;
        case CerealArchiveType::Enum::BINARY:
            (*std::get<std::shared_ptr<cereal::BinaryOutputArchive>>(ar))(args...);
            break;
        default:
            throw std::runtime_error("unknown cereal archive type!");
    }
}
}  // namespace ns_ikalibr

#endif  // IKALIBR_CEREAL_ARCHIVE_HELPER_HPP
