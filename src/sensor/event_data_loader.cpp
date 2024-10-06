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

#include "sensor/event_data_loader.h"
#include "ikalibr/PropheseeEventArray.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
EventDataLoader::EventDataLoader(EventModelType model)
    : _model(model) {}

EventDataLoader::Ptr EventDataLoader::GetLoader(const std::string& modelStr) {
    // try extract radar model
    EventModelType model;
    try {
        model = EnumCast::stringToEnum<EventModelType>(modelStr);
    } catch (...) {
        throw Status(Status::ERROR, EventModel::UnsupportedEventModelMsg(modelStr));
    }
    EventDataLoader::Ptr dataLoader;
    switch (model) {
        case EventModelType::PROPHESEE_EVENT:
            dataLoader = PropheseeEventDataLoader::Create(model);
            break;
        default:
            throw Status(Status::ERROR, EventModel::UnsupportedEventModelMsg(modelStr));
    }
    return dataLoader;
}

EventModelType EventDataLoader::GetEventModel() const { return _model; }

PropheseeEventDataLoader::PropheseeEventDataLoader(EventModelType model)
    : EventDataLoader(model) {}

PropheseeEventDataLoader::Ptr PropheseeEventDataLoader::Create(EventModelType model) {
    return std::make_shared<PropheseeEventDataLoader>(model);
}

EventArray::Ptr PropheseeEventDataLoader::UnpackData(const rosbag::MessageInstance& msgInstance) {
    ikalibr::PropheseeEventArrayPtr msg = msgInstance.instantiate<ikalibr::PropheseeEventArray>();

    CheckMessage<ikalibr::PropheseeEventArray>(msg);

    std::vector<Event::Ptr> events(msg->events.size());

    for (int i = 0; i < static_cast<int>(msg->events.size()); i++) {
        const auto& event = msg->events.at(i);
        events.at(i) =
            Event::Create(event.ts.toSec(), Eigen::Vector2d(event.x, event.y), event.polarity);
    }

    return EventArray::Create(msg->header.stamp.toSec(), events);
}
}  // namespace ns_ikalibr