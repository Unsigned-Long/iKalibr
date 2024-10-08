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

#include "solver/calib_solver.h"
#include "util/tqdm.h"
#include "spdlog/spdlog.h"
#include "calib/calib_data_manager.h"
#include "calib/calib_param_manager.h"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
void CalibSolver::InitPrepEventInertialAlign() const {
    if (!Configor::IsEventIntegrated()) {
        return;
    }
    // const auto &so3Spline = _splines->GetSo3Spline(Configor::Preference::SO3_SPLINE);
    // const auto &scaleSpline = _splines->GetRdSpline(Configor::Preference::SCALE_SPLINE);
    /**
     * we throw the head and tail data as the rotations from the fitted SO3 Spline in that range are
     * poor
     */
    // const double st = std::max(so3Spline.MinTime(), scaleSpline.MinTime()) +  // the max as start
    //                   Configor::Prior::TimeOffsetPadding;
    // const double et = std::min(so3Spline.MaxTime(), scaleSpline.MaxTime()) -  // the min as end
    //                   Configor::Prior::TimeOffsetPadding;

    /**
     * we first perform event-based feature tracking.
     * We first dedistort the events, then export them to files and use third-party software for
     * feature tracking.
     */
    for (const auto &[topic, eventMes] : _dataMagr->GetEventMeasurements()) {
        const auto &intri = _parMagr->INTRI.Camera.at(topic);
        const std::size_t EVENT_FRAME_NUM_THD = intri->imgHeight * intri->imgWidth / 2;
        /**
         * |--> 'outputSIter1'
         * |            |<- BATCH_TIME_WIN_THD ->|<- BATCH_TIME_WIN_THD ->|
         * ------------------------------------------------------------------
         * |<- BATCH_TIME_WIN_THD ->|<- BATCH_TIME_WIN_THD ->|
         * | data in this windown would be output for event-based feature tracking
         * |--> 'outputSIter2'
         */
        constexpr double BATCH_TIME_WIN_THD = 1.0;
        constexpr double BATCH_TIME_WIN_THD_HALF = BATCH_TIME_WIN_THD * 0.5;

        // this queue maintain two iterators
        std::queue<std::vector<EventArray::Ptr>::const_iterator> outputSIter;
        outputSIter.push(eventMes.cbegin());
        outputSIter.push(eventMes.cbegin());

        auto matSIter = eventMes.cbegin();
        std::size_t accumulatedEventNum = 0;
        cv::Mat eventFrameMat;

        // create a workspace for event-based feature tracking
        const std::string outputDir =
            Configor::DataStream::OutputPath + "/events/" + topic + "/haste_ws";
        if (!std::filesystem::exists(outputDir)) {
            if (!std::filesystem::create_directories(outputDir)) {
                throw Status(Status::CRITICAL,
                             "can not create output directory '{}' for event camera '{}'!!!",
                             outputDir, topic);
            }
        }
        int subEventDataIdx = 0;
        std::shared_ptr<tqdm> bar = std::make_shared<tqdm>();
        auto totalSubSeqNum = (eventMes.back()->GetTimestamp() - eventMes.front()->GetTimestamp()) /
                              BATCH_TIME_WIN_THD * 2;
        spdlog::info("output sub event data sequence for feature tracking for camera '{}'.", topic);
        for (auto curIter = matSIter; curIter != eventMes.cend(); ++curIter) {
            accumulatedEventNum += (*curIter)->GetEvents().size();
            if (accumulatedEventNum > EVENT_FRAME_NUM_THD) {
                /**
                 * If the number of events accumulates to a certain number, we construct it into an
                 * event frame. The event frame is used for visualization, and for calculating rough
                 * texture positions for feature tracking (haste-based). Note that this even frame
                 * is a distorted one
                 */
                eventFrameMat = EventArray::DrawRawEventFrame(matSIter, curIter, intri);
                cv::imshow("Event Frame", eventFrameMat);
                cv::waitKey(1);

                matSIter = curIter;
                accumulatedEventNum = 0;
            }

            if ((*curIter)->GetTimestamp() - (*outputSIter.back())->GetTimestamp() >
                BATCH_TIME_WIN_THD_HALF) {
                if ((*curIter)->GetTimestamp() - (*outputSIter.front())->GetTimestamp() >
                    BATCH_TIME_WIN_THD) {
                    bar->progress(subEventDataIdx, static_cast<int>(totalSubSeqNum));

                    // calculate rough texture positions for feature tracking (haste-based)
                    auto vertex = FindTexturePoints(eventFrameMat);

                    // the directory to save sub event data
                    const std::string subOutputDir =
                        outputDir + "/" + std::to_string(subEventDataIdx);
                    if (!std::filesystem::exists(subOutputDir)) {
                        if (!std::filesystem::create_directories(subOutputDir)) {
                            throw Status(Status::CRITICAL,
                                         "can not create sub output directory '{}' for event "
                                         "camera '{}' sub event data sequence '{}'!!!",
                                         subOutputDir, topic, subEventDataIdx);
                        }
                    }

                    SaveEventDataForFeatureTracking(outputSIter.front(),  // from
                                                    curIter,              // to
                                                    intri,                // intrinsics
                                                    vertex,               // position
                                                    subOutputDir          // directory
                    );
                    ++subEventDataIdx;
                }

                outputSIter.push(curIter);
                outputSIter.pop();
            }
        }
        bar->progress(subEventDataIdx, subEventDataIdx);
        bar->finish();
    }
    cv::destroyAllWindows();

    spdlog::warn("developing!!!");
    std::cin.get();
}
}  // namespace ns_ikalibr