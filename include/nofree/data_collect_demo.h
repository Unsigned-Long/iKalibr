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

#ifndef IKALIBR_DATA_COLLECT_DEMO_H
#define IKALIBR_DATA_COLLECT_DEMO_H

#include "util/utils.h"
#include "tiny-viewer/core/pose.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct DataCollectionMotionDemo {
private:
    constexpr static float Radius = 2.0;
    constexpr static float Time = 30.0;
    constexpr static int Frequency = 30;
    constexpr static int KeepTime = 3;

    constexpr static float CoordSizeMax = 0.4f;
    constexpr static float CoordSizeMin = 0.0f;
    constexpr static float CoordSizeRange = CoordSizeMax - CoordSizeMin;

    constexpr static float EightWidth = 1.0f;
    constexpr static float EightHeight = 0.6f;
    constexpr static float EightTime = 3.0f;

    constexpr static float CameraZBias = 0.5f;
    constexpr static float CameraYBias = 2.0f;

    constexpr static int KeepNum = KeepTime * Frequency;
    constexpr static float Dt = 1.0 / Frequency;
    constexpr static int PoseSize = Time / Dt;
    constexpr static float DTheta = 2 * M_PI / PoseSize;

public:
    static void Run(const std::array<std::vector<ns_viewer::Posef>, 2>& poseSeq);

    static std::array<std::vector<ns_viewer::Posef>, 2> GeneratePoseSeq();

    static std::array<std::vector<ns_viewer::Posef>, 2> GeneratePoseSeqLivox();
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_DATA_COLLECT_DEMO_H
