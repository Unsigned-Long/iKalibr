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

#include "core/haste_data_io.h"
#include "sensor/event.h"
#include "veta/camera/pinhole_brown.h"
#include "util/status.hpp"
#include "config/configor.h"
#include "filesystem"
#include "util/tqdm.h"
#include "core/event_trace_sac.h"
#include "core/feature_tracking.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

EventsInfo::SubBatch::SubBatch(int index,
                               double start_time,
                               double end_time,
                               std::size_t event_count)
    : index(index),
      start_time(start_time),
      end_time(end_time),
      event_count(event_count) {}

EventsInfo::EventsInfo(std::string topic,
                       std::string root_path,
                       double raw_start_time,
                       const std::vector<SubBatch> &batches)
    : topic(std::move(topic)),
      root_path(std::move(root_path)),
      raw_start_time(raw_start_time),
      batches(batches) {}

void EventsInfo::SaveToDisk(const std::string &filename,
                            CerealArchiveType::Enum archiveType) const {
    std::ofstream ofInfo(filename);
    auto ar = GetOutputArchiveVariant(ofInfo, archiveType);
    SerializeByOutputArchiveVariant(ar, archiveType, cereal::make_nvp("info", *this));
}

EventsInfo EventsInfo::LoadFromDisk(const std::string &filename,
                                    CerealArchiveType::Enum archiveType) {
    EventsInfo info;
    std::ifstream ifInfo(filename);
    auto ar = GetInputArchiveVariant(ifInfo, archiveType);
    SerializeByInputArchiveVariant(ar, archiveType, cereal::make_nvp("info", info));
    return info;
}

std::pair<std::string, EventsInfo::SubBatch> HASTEDataIO::SaveRawEventDataAsText(
    const std::vector<EventArray::Ptr>::const_iterator &fromIter,
    const std::vector<EventArray::Ptr>::const_iterator &toIter,
    const ns_veta::PinholeIntrinsic::Ptr &intri,
    const std::vector<Eigen::Vector2d> &seeds,
    double seedsTime,
    const std::string &subWS,
    int batchIdx) {
    // event.txt
    const std::string &eventsPath = subWS + "/events.txt";
    std::ofstream ofEvents(eventsPath, std::ios::out);
    std::stringstream buffer;
    std::size_t eventCount = 0;
    for (auto iter = fromIter; iter != toIter; ++iter) {
        const auto &events = (*iter)->GetEvents();
        for (const auto &event : events) {
            // todo: this is too too slow!!! modify haste to support binary data loading
            buffer << fmt::format("{:.9f} {} {} {}\n",                    // time, x, y, polarity
                                  event->GetTimestamp(),                  // time
                                  event->GetPos()(0),                     // x
                                  event->GetPos()(1),                     // y
                                  static_cast<int>(event->GetPolarity())  // polarity
            );
            ++eventCount;
        }
    }
    ofEvents << buffer.str();
    ofEvents.close();

    // calib.txt
    const std::string &calibPath = subWS + "/calib.txt";
    std::ofstream ofCalib(calibPath, std::ios::out);
    // since the events have been undistorted, the distortion parameters are all zeros
    // fx, fy, cx, cy, k1, k2, p1, p2, k3
    auto t2Intri = std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri);
    if (t2Intri == nullptr) {
        // the intrinsics of this camera is not 'ns_veta::PinholeIntrinsicBrownT2'
        throw Status(Status::CRITICAL,
                     "intrinsics of event camera is invalid, only the "
                     "'PinholeIntrinsicBrownT2' is supported currently!!!");
    }
    ofCalib << fmt::format(
        // fx, fy, cx, cy, k1, k2, p1, p2, k3
        "{:.3f} {:.3f} {:.3f} {:.3f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}\n",  // params
        t2Intri->GetParams().at(0), t2Intri->GetParams().at(1),              // fx, fy
        t2Intri->GetParams().at(2), t2Intri->GetParams().at(3),              // cx, cy
        t2Intri->GetParams().at(4), t2Intri->GetParams().at(5),              // k1, k2
        t2Intri->GetParams().at(7), t2Intri->GetParams().at(8),              // p1, p2
        t2Intri->GetParams().at(6)                                           // k3
    );
    ofCalib.close();

    // seeds.txt
    const std::string &seedsPath = subWS + "/seeds.txt";
    std::ofstream ofSeeds(seedsPath, std::ios::out);
    buffer = std::stringstream();
    for (int id = 0; id != static_cast<int>(seeds.size()); ++id) {
        const Eigen::Vector2d &seed = intri->GetUndistoPixel(seeds.at(id));
        // todo: this is too too slow!!! modify haste to support binary data loading
        buffer << fmt::format("{:.9f},{:.3f},{:.3f},0.0,{}\n",  // t, x, y, theta, id
                              seedsTime,                        // t
                              seed(0), seed(1),                 // x, y
                              id                                // id
        );
    }
    ofSeeds << buffer.str();
    ofSeeds.close();

    // command
    const std::string hasteProg = "/home/csl/Software/haste/build/tracking_app_file";
    auto command = fmt::format(
        "{} "
        "-events_file={} "  // Plain text file with events
        "-seeds_file={} "   // Plain text file with several initial tracking seeds
        "-tracker_type={} "  // correlation|haste_correlation|haste_correlation_star|haste_difference|haste_difference_star
        "-centered_initialization={} "  // Force tracker to be centered/non-centered initialized
        "-camera_params_file={} "       // Load pinhole camera calibration model
        "-camera_size={}x{} "           // Set image sensor resolution
        "-visualize={} "                // Visualize internal tracker state
        "-output_file={}",              // Write tracking results to file
        hasteProg, eventsPath, seedsPath, "haste_correlation_star", true, calibPath,
        intri->imgWidth, intri->imgHeight, false, subWS + "/haste_results.txt");

    auto batchInfo = EventsInfo::SubBatch(batchIdx, (*fromIter)->GetTimestamp(),
                                          (*toIter)->GetTimestamp(), eventCount);

    return {command, batchInfo};
}

std::pair<std::string, EventsInfo::SubBatch> HASTEDataIO::SaveRawEventDataAsBinary(
    const std::vector<EventArrayPtr>::const_iterator &fromIter,
    const std::vector<EventArrayPtr>::const_iterator &toIter,
    const ns_veta::PinholeIntrinsicPtr &intri,
    const std::vector<Eigen::Vector2d> &seeds,
    double seedsTime,
    const std::string &subWS,
    int batchIdx) {
    // event.bin
    const std::string &eventsPath = subWS + "/events.bin";
    std::ofstream ofEvents(eventsPath, std::ios::binary);
    std::size_t eventCount = 0;
    for (auto iter = fromIter; iter != toIter; ++iter) {
        const auto &events = (*iter)->GetEvents();
        for (const auto &event : events) {
            // time (float), x (uint16_t), y (uint16_t), polarity (boolean)
            auto time = static_cast<float>(event->GetTimestamp());
            std::uint16_t x = event->GetPos()(0);
            std::uint16_t y = event->GetPos()(1);
            bool polarity = event->GetPolarity();

            ofEvents.write(reinterpret_cast<const char *>(&time), sizeof(time));
            ofEvents.write(reinterpret_cast<const char *>(&x), sizeof(x));
            ofEvents.write(reinterpret_cast<const char *>(&y), sizeof(y));
            ofEvents.write(reinterpret_cast<const char *>(&polarity), sizeof(polarity));

            ++eventCount;
        }
    }
    ofEvents.close();

    // calib.bin
    const std::string &calibPath = subWS + "/calib.bin";
    std::ofstream ofCalib(calibPath, std::ios::binary);
    // since the events have been undistorted, the distortion parameters are all zeros
    // fx, fy, cx, cy, k1, k2, p1, p2, k3
    auto t2Intri = std::dynamic_pointer_cast<ns_veta::PinholeIntrinsicBrownT2>(intri);
    if (t2Intri == nullptr) {
        // the intrinsics of this camera is not 'ns_veta::PinholeIntrinsicBrownT2'
        throw Status(Status::CRITICAL,
                     "intrinsics of event camera is invalid, only the "
                     "'PinholeIntrinsicBrownT2' is supported currently!!!");
    }
    float fx = static_cast<float>(t2Intri->GetParams().at(0)),
          fy = static_cast<float>(t2Intri->GetParams().at(1));
    float cx = static_cast<float>(t2Intri->GetParams().at(2)),
          cy = static_cast<float>(t2Intri->GetParams().at(3));
    float k1 = static_cast<float>(t2Intri->GetParams().at(4)),
          k2 = static_cast<float>(t2Intri->GetParams().at(5));
    float p1 = static_cast<float>(t2Intri->GetParams().at(7)),
          p2 = static_cast<float>(t2Intri->GetParams().at(8));
    float k3 = static_cast<float>(t2Intri->GetParams().at(6));
    // fx, fy, cx, cy, k1, k2, p1, p2, k3 (all are float format)
    ofCalib.write(reinterpret_cast<const char *>(&fx), sizeof(fx));
    ofCalib.write(reinterpret_cast<const char *>(&fy), sizeof(fy));
    ofCalib.write(reinterpret_cast<const char *>(&cx), sizeof(cx));
    ofCalib.write(reinterpret_cast<const char *>(&cy), sizeof(cy));
    ofCalib.write(reinterpret_cast<const char *>(&k1), sizeof(k1));
    ofCalib.write(reinterpret_cast<const char *>(&k2), sizeof(k2));
    ofCalib.write(reinterpret_cast<const char *>(&p1), sizeof(p1));
    ofCalib.write(reinterpret_cast<const char *>(&p2), sizeof(p2));
    ofCalib.write(reinterpret_cast<const char *>(&k3), sizeof(k3));
    ofCalib.close();

    // seeds.bin
    const std::string &seedsPath = subWS + "/seeds.bin";
    std::ofstream ofSeeds(seedsPath, std::ios::binary);
    for (int i = 0; i != static_cast<int>(seeds.size()); ++i) {
        const Eigen::Vector2d &seed = intri->GetUndistoPixel(seeds.at(i));
        auto time = static_cast<float>(seedsTime);
        auto x = static_cast<float>(seed(0));
        auto y = static_cast<float>(seed(1));
        float theta = 0.0f;
        std::uint64_t id = i;

        // t (float), x (float), y (float), theta (float), id (uint64_t)
        ofSeeds.write(reinterpret_cast<const char *>(&time), sizeof(time));
        ofSeeds.write(reinterpret_cast<const char *>(&x), sizeof(x));
        ofSeeds.write(reinterpret_cast<const char *>(&y), sizeof(y));
        ofSeeds.write(reinterpret_cast<const char *>(&theta), sizeof(theta));
        ofSeeds.write(reinterpret_cast<const char *>(&id), sizeof(id));
    }
    ofSeeds.close();

    // command
    const std::string hasteProg = "/home/csl/Software/haste/build/tracking_app_binary";
    auto command = fmt::format(
        "{} "
        "-events_file={} "  // Plain text file with events
        "-seeds_file={} "   // Plain text file with several initial tracking seeds
        "-tracker_type={} "  // correlation|haste_correlation|haste_correlation_star|haste_difference|haste_difference_star
        "-centered_initialization={} "  // Force tracker to be centered/non-centered initialized
        "-camera_params_file={} "       // Load pinhole camera calibration model
        "-camera_size={}x{} "           // Set image sensor resolution
        "-visualize={} "                // Visualize internal tracker state
        "-output_file={}",              // Write tracking results to file
        hasteProg, eventsPath, seedsPath, "haste_correlation_star", true, calibPath,
        intri->imgWidth, intri->imgHeight, false, subWS + "/haste_results.bin");

    auto batchInfo = EventsInfo::SubBatch(batchIdx, (*fromIter)->GetTimestamp(),
                                          (*toIter)->GetTimestamp(), eventCount);

    return {command, batchInfo};
}

std::optional<HASTEDataIO::TrackingResultsType> HASTEDataIO::TryLoadHASTEResultsFromTXT(
    const EventsInfo &info, const ns_veta::PinholeIntrinsicPtr &intri, double newRawStartTime) {
    if (info.batches.empty()) {
        spdlog::warn("there is no any batch in '{}'!!!", info.root_path);
        return {};
    }
    std::shared_ptr<tqdm> bar = std::make_shared<tqdm>();
    // batch index, feature id, tracking list
    TrackingResultsType trackResults;
    std::size_t trackedFeatCount = 0;
    double minTime = std::numeric_limits<double>::max();
    double maxTime = std::numeric_limits<double>::min();
    for (const auto &batch : info.batches) {
        bar->progress(batch.index, static_cast<int>(info.batches.size()));

        const std::string resultsPath =
            fmt::format("{}/{}/haste_results.txt", info.root_path, batch.index);
        if (!std::filesystem::exists(resultsPath)) {
            spdlog::warn("haste result file of batch indexed as {} is not found!!!", batch.index);
            continue;
        }

        std::ifstream ifResults(resultsPath, std::ios::in);
        std::string strLine;
        auto &tracking = trackResults[batch.index];

        while (std::getline(ifResults, strLine)) {
            auto strVec = SplitString(strLine, ',');
            if (strVec.size() != 5) {
                spdlog::warn("the feature tracking result item '{}' may broken!!!", strLine);
                continue;
            }
            Feature feature;
            feature.timestamp = std::stod(strVec.at(0)) + info.raw_start_time - newRawStartTime;
            // results from haste are features that have been undistorted
            feature.undistorted = cv::Point2f(std::stof(strVec.at(1)), std::stof(strVec.at(2)));
            Eigen::Vector2d rawFeat = intri->GetDistoPixel(feature.Undistorted().cast<double>());
            feature.raw = cv::Point2f(rawFeat(0), rawFeat(1));

            int id = std::stoi(strVec.at(4));
            tracking[id].push_back(std::make_shared<Feature>(feature));

            ++trackedFeatCount;
            if (feature.timestamp < minTime) {
                minTime = feature.timestamp;
            }
            if (feature.timestamp > maxTime) {
                maxTime = feature.timestamp;
            }
        }
        ifResults.close();
    }
    bar->finish();
    spdlog::info(
        "load tracking results finished, batch count: {}, total tracked feature count: {}, time "
        "span from '{:.5f}' to '{:.5f}', time range: '{:.5f}'",
        trackResults.size(), trackedFeatCount, minTime, maxTime, maxTime - minTime);
    return trackResults;
}

std::optional<HASTEDataIO::TrackingResultsType> HASTEDataIO::TryLoadHASTEResultsFromBinary(
    const EventsInfo &info, const ns_veta::PinholeIntrinsicPtr &intri, double newRawStartTime) {
    if (info.batches.empty()) {
        spdlog::warn("there is no any batch in '{}'!!!", info.root_path);
        return {};
    }
    std::shared_ptr<tqdm> bar = std::make_shared<tqdm>();
    // batch index, feature id, tracking list
    TrackingResultsType trackResults;
    std::size_t trackedFeatCount = 0;
    double minTime = std::numeric_limits<double>::max();
    double maxTime = std::numeric_limits<double>::min();
    for (const auto &batch : info.batches) {
        bar->progress(batch.index, static_cast<int>(info.batches.size()));

        const std::string resultsPath =
            fmt::format("{}/{}/haste_results.bin", info.root_path, batch.index);
        if (!std::filesystem::exists(resultsPath)) {
            spdlog::warn("haste result file of batch indexed as {} is not found!!!", batch.index);
            continue;
        }

        std::ifstream ifResults(resultsPath, std::ios::binary);

        float time;
        float x;
        float y;
        float theta;
        std::uint64_t id;

        auto &tracking = trackResults[batch.index];

        while (ifResults.read(reinterpret_cast<char *>(&time), sizeof(time))) {
            ifResults.read(reinterpret_cast<char *>(&x), sizeof(x));
            ifResults.read(reinterpret_cast<char *>(&y), sizeof(y));
            ifResults.read(reinterpret_cast<char *>(&theta), sizeof(theta));
            ifResults.read(reinterpret_cast<char *>(&id), sizeof(id));

            Feature feature;
            feature.timestamp = time + info.raw_start_time - newRawStartTime;
            // results from haste are features that have been undistorted
            feature.undistorted = cv::Point2f{x, y};
            Eigen::Vector2d rawFeat = intri->GetDistoPixel(feature.Undistorted().cast<double>());
            feature.raw = cv::Point2f(rawFeat(0), rawFeat(1));

            tracking[static_cast<int>(id)].push_back(std::make_shared<Feature>(feature));

            ++trackedFeatCount;
            if (feature.timestamp < minTime) {
                minTime = feature.timestamp;
            }
            if (feature.timestamp > maxTime) {
                maxTime = feature.timestamp;
            }
        }
        ifResults.close();
    }
    bar->finish();
    spdlog::info(
        "load tracking results finished, batch count: {}, total tracked feature count: {}, time "
        "span from '{:.5f}' to '{:.5f}', time range: '{:.5f}'",
        trackResults.size(), trackedFeatCount, minTime, maxTime, maxTime - minTime);
    return trackResults;
}

void HASTEDataIO::SaveEventsInfo(const EventsInfo &info, const std::string &ws) {
    info.SaveToDisk(ws + "/info" + Configor::GetFormatExtension(),
                    Configor::Preference::OutputDataFormat);
}

std::optional<EventsInfo> HASTEDataIO::TryLoadEventsInfo(const std::string &ws) {
    const std::string filepath = ws + "/info" + Configor::GetFormatExtension();
    if (!std::filesystem::exists(filepath)) {
        return {};
    }
    return EventsInfo::LoadFromDisk(filepath, Configor::Preference::OutputDataFormat);
}

void EventTrackingFilter::FilterByTrackingLength(FeatureVecMap &tracking,
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

void EventTrackingFilter::FilterByTraceFittingSAC(FeatureVecMap &tracking, double thd) {
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
                      [](const std::shared_ptr<Feature> &f1, const std::shared_ptr<Feature> &f2) {
                          return f1->timestamp < f2->timestamp;
                      });
            ++it;
        } else {
            it = tracking.erase(it);
        }
    }
}

void EventTrackingFilter::FilterByTrackingAge(FeatureVecMap &tracking,
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
void EventTrackingFilter::FilterByTrackingFreq(FeatureVecMap &tracking,
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