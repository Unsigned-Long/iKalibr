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

#ifndef EVENT_PREPROCESSING_H
#define EVENT_PREPROCESSING_H

#include "util/utils.h"
#include "sensor/event.h"
#include "opencv2/imgproc.hpp"
#include "core/visual_distortion.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_veta {
struct PinholeIntrinsic;
using PinholeIntrinsicPtr = std::shared_ptr<PinholeIntrinsic>;
}  // namespace ns_veta

namespace ns_ikalibr {

struct VisualUndistortionMap;
using VisualUndistortionMapPtr = std::shared_ptr<VisualUndistortionMap>;

class ActiveEventSurface {
public:
    using Ptr = std::shared_ptr<ActiveEventSurface>;

private:
    const double FILTER_THD;
    ns_veta::PinholeIntrinsicPtr _intri;
    VisualUndistortionMapPtr _undistoMap;

    Eigen::MatrixXd _sae[2];        // save sae
    Eigen::MatrixXd _saeLatest[2];  // save previous sae

    cv::Mat _eventImgMat;

public:
    explicit ActiveEventSurface(const ns_veta::PinholeIntrinsicPtr& intri, double filterThd = 0.01);

    static Ptr Create(const ns_veta::PinholeIntrinsicPtr& intri, double filterThd = 0.01);

    void GrabEvent(const Event::Ptr& event, bool drawEventMat = false);

    void GrabEvent(const EventArray::Ptr& events, bool drawEventMat = false);

    [[nodiscard]] cv::Mat GetEventImgMat(bool resetMat, bool undistoMat = false);

    cv::Mat TimeSurface(double curTime,
                        bool ignorePolarity = false,
                        bool undistoMat = false,
                        int medianBlurKernelSize = 0,
                        double decaySec = 0.02);

    cv::Mat RawTimeSurface(bool ignorePolarity = false, bool undistoMat = false);
};

}  // namespace ns_ikalibr

#endif  // EVENT_PREPROCESSING_H
