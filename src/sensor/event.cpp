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

#include "sensor/event.h"
#include "veta/camera/pinhole.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
Event::Event(double timestamp, Eigen::Vector2d pos, bool polarity)
    : _timestamp(timestamp),
      _pos(std::move(pos)),
      _polarity(polarity) {}

Event::Ptr Event::Create(double timestamp, const Eigen::Vector2d& pos, bool polarity) {
    return std::make_shared<Event>(timestamp, pos, polarity);
}

double Event::GetTimestamp() const { return _timestamp; }

void Event::SetTimestamp(double timestamp) { _timestamp = timestamp; }

Eigen::Vector2d Event::GetPos() const { return _pos; }

bool Event::GetPolarity() const { return _polarity; }

EventArray::EventArray(double timestamp, const std::vector<Event::Ptr>& events)
    : _timestamp(timestamp),
      _events(events) {}

EventArray::Ptr EventArray::Create(double timestamp, const std::vector<Event::Ptr>& events) {
    return std::make_shared<EventArray>(timestamp, events);
}

double EventArray::GetTimestamp() const { return _timestamp; }

std::vector<Event::Ptr> EventArray::GetEvents() const { return _events; }

void EventArray::SetTimestamp(double timestamp) { _timestamp = timestamp; }

cv::Mat EventArray::DrawRawEventFrame(const ns_veta::PinholeIntrinsic::Ptr& intri) const {
    cv::Mat eventFrame =
        cv::Mat(static_cast<int>(intri->imgHeight), static_cast<int>(intri->imgWidth), CV_8UC3,
                cv::Scalar(255, 255, 255));

    for (const auto& event : _events) {
        cv::Vec3b color;
        if (event->GetPolarity()) {
            // red
            color = cv::Vec3b(0, 0, 255);
        } else {
            // green
            color = cv::Vec3b(255, 0, 0);
        }
        eventFrame.at<cv::Vec3b>(static_cast<int>(event->GetPos()(1)),
                                 static_cast<int>(event->GetPos()(0))) = color;
    }
    return eventFrame;
}

cv::Mat EventArray::DrawRawEventFrame(const std::vector<Ptr>::const_iterator& sIter,
                                      const std::vector<Ptr>::const_iterator& eIter,
                                      const ns_veta::PinholeIntrinsicPtr& intri) {
    cv::Mat eventFrame =
        cv::Mat(static_cast<int>(intri->imgHeight), static_cast<int>(intri->imgWidth), CV_8UC3,
                cv::Scalar(255, 255, 255));

    for (auto iter = sIter; iter != eIter; ++iter) {
        for (const auto& event : (*iter)->GetEvents()) {
            cv::Vec3b color;
            if (event->GetPolarity()) {
                // red
                color = cv::Vec3b(0, 0, 255);
            } else {
                // green
                color = cv::Vec3b(255, 0, 0);
            }
            eventFrame.at<cv::Vec3b>(static_cast<int>(event->GetPos()(1)),
                                     static_cast<int>(event->GetPos()(0))) = color;
        }
    }

    return eventFrame;
}

}  // namespace ns_ikalibr
