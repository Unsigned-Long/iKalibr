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

#ifndef EVENT_DATA_LOADER_H
#define EVENT_DATA_LOADER_H

#include "sensor/event.h"
#include "sensor/sensor_model.h"
#include "rosbag/message_instance.h"
#include "util/enum_cast.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
class EventDataLoader {
public:
    using Ptr = std::shared_ptr<EventDataLoader>;

protected:
    EventModelType _model;

public:
    explicit EventDataLoader(EventModelType model);

    virtual EventArray::Ptr UnpackData(const rosbag::MessageInstance &msgInstance) = 0;

    static EventDataLoader::Ptr GetLoader(const std::string &modelStr);

    [[nodiscard]] EventModelType GetEventModel() const;

    virtual ~EventDataLoader() = default;

protected:
    template <class MsgType>
    void CheckMessage(typename MsgType::ConstPtr msg) {
        if (msg == nullptr) {
            throw std::runtime_error(
                "Wrong sensor model: '" + std::string(EnumCast::enumToString(GetEventModel())) +
                "' for event cameras! It's incompatible with the type of ros message to load in!");
        }
    }
};

class PropheseeEventDataLoader : public EventDataLoader {
public:
    using Ptr = std::shared_ptr<PropheseeEventDataLoader>;

public:
    explicit PropheseeEventDataLoader(EventModelType model);

    static Ptr Create(EventModelType model);

    EventArray::Ptr UnpackData(const rosbag::MessageInstance &msgInstance) override;
};
}  // namespace ns_ikalibr

#endif  // EVENT_DATA_LOADER_H
