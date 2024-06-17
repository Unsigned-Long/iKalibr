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

#ifndef IKALIBR_LIDAR_H
#define IKALIBR_LIDAR_H

#include "util/utils.h"
#include "util/cloud_define.hpp"
#include "ctraj/utils/macros.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct LiDARFrame {
public:
    using Ptr = std::shared_ptr<LiDARFrame>;

private:
    // the timestamp of this lidar scan
    double _timestamp;
    // the lidar scan [x, y, z, timestamp]
    IKalibrPointCloud::Ptr _scan;

public:
    // constructor
    explicit LiDARFrame(double timestamp = INVALID_TIME_STAMP,
                        IKalibrPointCloud::Ptr scan = boost::make_shared<IKalibrPointCloud>());

    // creator
    static LiDARFrame::Ptr Create(
        double timestamp = INVALID_TIME_STAMP,
        const IKalibrPointCloud::Ptr &scan = boost::make_shared<IKalibrPointCloud>());

    // access
    [[nodiscard]] IKalibrPointCloud::Ptr GetScan() const;

    [[nodiscard]] double GetTimestamp() const;

    void SetTimestamp(double timestamp);

    friend std::ostream &operator<<(std::ostream &os, const LiDARFrame &frame);
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_LIDAR_H
