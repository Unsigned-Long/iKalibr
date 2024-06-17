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

#ifndef IKALIBR_ENUM_CAST_HPP
#define IKALIBR_ENUM_CAST_HPP

#include "exception"
#include "magic_enum.hpp"
#include "string"
#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct EnumCast {
    template <class EnumType>
    static constexpr auto stringToEnum(const std::string &enumStr) {
        if (auto color = magic_enum::enum_cast<EnumType>(enumStr); color.has_value()) {
            return color.value();
        } else {
            throw std::runtime_error("'EnumCast::stringToEnum' cast failed");
        }
    }

    template <class EnumType>
    static constexpr auto integerToEnum(int enumValue) {
        if (auto color = magic_enum::enum_cast<EnumType>(enumValue); color.has_value()) {
            return color.value();
        } else {
            throw std::runtime_error("'EnumCast::integerToEnum' cast failed");
        }
    }

    template <class EnumType>
    static constexpr auto enumToInteger(EnumType enumType) {
        return magic_enum::enum_integer(enumType);
    }

    template <class EnumType>
    static constexpr auto enumToString(EnumType enumType) {
        return magic_enum::enum_name(enumType);
    }

    template <class EnumType>
    static constexpr auto stringToInteger(const std::string &enumStr) {
        return enumToInteger(stringToEnum<EnumType>(enumStr));
    }

    template <class EnumType>
    static constexpr auto integerToString(int enumValue) {
        return enumToString(integerToEnum<EnumType>(enumValue));
    }
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_ENUM_CAST_HPP
