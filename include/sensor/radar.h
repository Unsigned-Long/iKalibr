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

#ifndef IKALIBR_RADAR_H
#define IKALIBR_RADAR_H

#include "util/utils.h"
#include "ctraj/utils/macros.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct RadarTarget {
public:
    using Ptr = std::shared_ptr<RadarTarget>;

private:
    // the timestamp of this frame
    double _timestamp;

    Eigen::Vector3d _target;
    double _radialVel;

    double _range, _invRange;

public:
    /**
     * @attention rawMes: [ range | theta | phi | target radial vel with respect to radar in frame
     * {R} ]
     */
    explicit RadarTarget(double timestamp = INVALID_TIME_STAMP,
                         const Eigen::Vector4d &rawMes = Eigen::Vector4d::Zero());

    /**
     * @attention rawMes: [ xyz | target radial vel with respect to radar in frame {R} ]
     */
    explicit RadarTarget(double timestamp, Eigen::Vector3d target, double radialVel);

    /**
     * @attention rawMes: [ range | theta | phi | target radial vel with respect to radar in frame
     * {R} ]
     */
    static RadarTarget::Ptr Create(double timestamp = INVALID_TIME_STAMP,
                                   const Eigen::Vector4d &rawMes = Eigen::Vector4d::Zero());

    /**
     * @attention rawMes: [ xyz | target radial vel with respect to radar in frame {R} ]
     */
    static RadarTarget::Ptr Create(double timestamp,
                                   const Eigen::Vector3d &target,
                                   double radialVel);

    // access
    [[nodiscard]] double GetTimestamp() const;

    void SetTimestamp(double timestamp);

    [[nodiscard]] const Eigen::Vector3d &GetTargetXYZ() const;

    [[nodiscard]] Eigen::Vector3d GetTargetRTP() const;

    [[nodiscard]] double GetRadialVelocity() const;

    [[nodiscard]] double GetRange() const;

    [[nodiscard]] double GetInvRange() const;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    template <class Archive>
    void serialize(Archive &ar) {
        ar(cereal::make_nvp("timestamp", _timestamp), cereal::make_nvp("target", _target),
           cereal::make_nvp("radial_vel", _radialVel));
    }
};

struct RadarTargetArray {
public:
    using Ptr = std::shared_ptr<RadarTargetArray>;

private:
    // the timestamp of this array
    double _timestamp;
    std::vector<RadarTarget::Ptr> _targets;

public:
    explicit RadarTargetArray(double timestamp = INVALID_TIME_STAMP,
                              const std::vector<RadarTarget::Ptr> &targets = {});

    static RadarTargetArray::Ptr Create(double timestamp = INVALID_TIME_STAMP,
                                        const std::vector<RadarTarget::Ptr> &targets = {});

    [[nodiscard]] double GetTimestamp() const;

    void SetTimestamp(double timestamp);

    [[nodiscard]] const std::vector<RadarTarget::Ptr> &GetTargets() const;

    // save radar frames sequence to disk
    static bool SaveTargetArraysToDisk(const std::string &filename,
                                       const std::vector<RadarTargetArray::Ptr> &arrays,
                                       int precision = 10);

    Eigen::Vector3d RadarVelocityFromStaticTargetArray();

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    template <class Archive>
    void serialize(Archive &ar) {
        ar(cereal::make_nvp("timestamp", _timestamp), cereal::make_nvp("targets", _targets));
    }
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_RADAR_H
