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

#ifndef EVENT_H
#define EVENT_H

#include "util/utils.h"
#include "ctraj/utils/macros.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_veta {
struct PinholeIntrinsic;
using PinholeIntrinsicPtr = std::shared_ptr<PinholeIntrinsic>;
}  // namespace ns_veta

namespace ns_ikalibr {
class Event {
public:
    using Ptr = std::shared_ptr<Event>;

private:
    // the timestamp of this event
    double _timestamp;
    Eigen::Vector2d _pos;
    bool _polarity;

public:
    explicit Event(double timestamp = INVALID_TIME_STAMP,
                   Eigen::Vector2d pos = Eigen::Vector2d::Zero(),
                   bool polarity = {});

    static Ptr Create(double timestamp = INVALID_TIME_STAMP,
                      const Eigen::Vector2d& pos = Eigen::Vector2d::Zero(),
                      bool polarity = {});

    [[nodiscard]] double GetTimestamp() const;

    void SetTimestamp(double timestamp);

    [[nodiscard]] Eigen::Vector2d GetPos() const;

    [[nodiscard]] bool GetPolarity() const;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    template <class Archive>
    void serialize(Archive& ar) {
        ar(cereal::make_nvp("timestamp", _timestamp), cereal::make_nvp("pos", _pos),
           cereal::make_nvp("polarity", _polarity));
    }
};

class EventArray {
public:
    using Ptr = std::shared_ptr<EventArray>;

private:
    double _timestamp;
    std::vector<Event::Ptr> _events;

public:
    explicit EventArray(double timestamp = INVALID_TIME_STAMP,
                        const std::vector<Event::Ptr>& events = {});

    static Ptr Create(double timestamp = INVALID_TIME_STAMP,
                      const std::vector<Event::Ptr>& events = {});

    [[nodiscard]] double GetTimestamp() const;

    [[nodiscard]] std::vector<Event::Ptr> GetEvents() const;

    void SetTimestamp(double timestamp);

    [[nodiscard]] cv::Mat DrawRawEventFrame(const ns_veta::PinholeIntrinsicPtr& intri) const;

    static cv::Mat DrawRawEventFrame(const std::vector<Ptr>::const_iterator& sIter,
                                     const std::vector<Ptr>::const_iterator& eIter,
                                     const ns_veta::PinholeIntrinsicPtr& intri);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    template <class Archive>
    void serialize(Archive& ar) {
        ar(cereal::make_nvp("timestamp", _timestamp), cereal::make_nvp("events", _events));
    }
};
}  // namespace ns_ikalibr

#endif  // EVENT_H
