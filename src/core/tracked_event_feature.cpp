// iKalibr: Unified Targetless Spatiotemporal Calibration Framework
// Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China
// https://github.com/Unsigned-Long/iKalibr.git
// Author: Shuolong Chen (shlchen@whu.edu.cn)
// GitHub: https://github.com/Unsigned-Long
//  ORCID: 0000-0002-5283-9057
// Purpose: See .h/.hpp file.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * The names of its contributors can not be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
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

#include "core/tracked_event_feature.h"
#include "utility"
#include "core/event_trace_sac.h"

namespace ns_ikalibr {
EventFeature::EventFeature(double timestamp, Eigen::Vector2d pos)
    : timestamp(timestamp),
      pos(std::move(pos)) {}

void EventTrackingFilter::FilterByTrackingLength(EventFeatTrackingBatch &tracking,
                                                 double acceptedTrackedThdCompBest) {
    std::size_t maxLength = 0;
    for (const auto &[featId, trackingList] : tracking) {
        if (trackingList.size() > maxLength) {
            maxLength = trackingList.size();
        }
    }
    // auto oldSize = tracking.size();
    auto acceptedMinLength = static_cast<double>(maxLength) * acceptedTrackedThdCompBest;
    for (auto it = tracking.begin(); it != tracking.end();) {
        if (it->second.size() < static_cast<std::size_t>(acceptedMinLength)) {
            it = tracking.erase(it);
        } else {
            ++it;
        }
    }
}

void EventTrackingFilter::FilterByTraceFittingSAC(EventFeatTrackingBatch &tracking, double thd) {
    for (auto it = tracking.begin(); it != tracking.end();) {
        auto &trackingList = it->second;
        if (trackingList.size() < 3) {
            it = tracking.erase(it);
            continue;
        }
        auto res = EventTrackingTraceSacProblem::EventTrackingTraceSac(trackingList, thd);
        if (res.first != nullptr) {
            trackingList = res.second;
            std::sort(trackingList.begin(), trackingList.end(),
                      [](const std::shared_ptr<EventFeature> &f1,
                         const std::shared_ptr<EventFeature> &f2) {
                          return f1->timestamp < f2->timestamp;
                      });
            ++it;
        } else {
            it = tracking.erase(it);
        }
    }
}

void EventTrackingFilter::FilterByTrackingAge(EventFeatTrackingBatch &tracking,
                                              double acceptedTrackedThdCompBest) {
    double maxAge = 0.0;
    for (const auto &[featId, trackingList] : tracking) {
        const double age = trackingList.back()->timestamp - trackingList.front()->timestamp;
        if (age > maxAge) {
            maxAge = age;
        }
    }
    auto acceptedMinAge = maxAge * acceptedTrackedThdCompBest;
    for (auto it = tracking.begin(); it != tracking.end();) {
        const auto &trackingList = it->second;
        const double age = trackingList.back()->timestamp - trackingList.front()->timestamp;

        if (age < acceptedMinAge) {
            it = tracking.erase(it);
        } else {
            ++it;
        }
    }
}
void EventTrackingFilter::FilterByTrackingFreq(EventFeatTrackingBatch &tracking,
                                               double acceptedTrackedThdCompBest) {
    double maxFreq = 0.0;
    for (const auto &[featId, trackingList] : tracking) {
        const double age = trackingList.back()->timestamp - trackingList.front()->timestamp;
        const double freq = static_cast<double>(trackingList.size()) / age;
        if (freq > maxFreq) {
            maxFreq = freq;
        }
    }
    auto acceptedMinFreq = maxFreq * acceptedTrackedThdCompBest;
    for (auto it = tracking.begin(); it != tracking.end();) {
        const auto &trackingList = it->second;
        const double age = trackingList.back()->timestamp - trackingList.front()->timestamp;
        const double freq = static_cast<double>(trackingList.size()) / age;

        if (freq < acceptedMinFreq) {
            it = tracking.erase(it);
        } else {
            ++it;
        }
    }
}

}  // namespace ns_ikalibr