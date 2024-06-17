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

#include "nofree/logo_svg.h"
#include "fstream"
#include "spdlog/fmt/fmt.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void LoGoMaker::SaveToFile(const std::string &filename) {
    std::string str =
        "_|  _|    _|            _|  _|  _|                \n"
        "    _|  _|      _|_|_|  _|      _|_|_|    _|  _|_|\n"
        "_|  _|_|      _|    _|  _|  _|  _|    _|  _|_|    \n"
        "_|  _|  _|    _|    _|  _|  _|  _|    _|  _|      \n"
        "_|  _|    _|    _|_|_|  _|  _|  _|_|_|    _|      ";
    float x = 1.0, y = 1.0;
    constexpr float w = 10, h = 10;
    std::ofstream file(filename, std::ios::out);
    file << R"(<svg width="252" height="52" xmlns="http://www.w3.org/2000/svg">)" << std::endl;
    for (size_t i = 0; i < str.size() - 1; i++) {
        auto c = str.at(i);
        auto c2 = str.at(i + 1);
        if (c == '_' && c2 == '|') {
            file << DrawCube(x, y, w, h) << std::endl;
            x += w;
            i += 1;
        } else if (c == '\n') {
            y += h;
            x = 1.0;
        } else {
            x += w * 0.5;
        }
    }
    file << "</svg>" << std::endl;
    file.close();
}

std::string LoGoMaker::DrawCube(float x, float y, float w, float h) {
    static std::default_random_engine e;
    static std::uniform_real_distribution<float> u1(0.6, 1.0);
    static std::uniform_real_distribution<float> u2(0.0, 0.4);
    static std::uniform_real_distribution<float> s(1.0, 2.0);
    std::string color;
    if (x == 1) {
        color = "red";
    } else {
        color = "green";
    }
    return fmt::format(
        "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" stroke=\"black\" stroke-width=\"0.5\" "
        "fill=\"{}\" rx=\"3\" ry=\"3\">\n"
        "    <animate attributeName=\"opacity\" values=\"{};{};{}\" dur=\"{}s\" "
        "repeatCount=\"indefinite\" />\n"
        "</rect>",
        x, y, w, h, color, u1(e), u2(e), u1(e), s(e));
}
}  // namespace ns_ikalibr
