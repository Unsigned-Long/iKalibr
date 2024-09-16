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

#ifndef IKALIBR_CALIB_DATA_MANAGER_H
#define IKALIBR_CALIB_DATA_MANAGER_H

#include "config/configor.h"
#include "sensor/camera.h"
#include "sensor/imu.h"
#include "sensor/lidar.h"
#include "sensor/radar.h"
#include "sensor/rgbd.h"
#include "util/status.hpp"
#include "veta/veta.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct OpticalFlowTripleTrace;
using OpticalFlowTripleTracePtr = std::shared_ptr<OpticalFlowTripleTrace>;

class CalibDataManager {
public:
    using Ptr = std::shared_ptr<CalibDataManager>;

private:
    std::map<std::string, std::vector<IMUFrame::Ptr>> _imuMes;
    std::map<std::string, std::vector<RadarTargetArray::Ptr>> _radarMes;
    std::map<std::string, std::vector<LiDARFrame::Ptr>> _lidarMes;
    std::map<std::string, std::vector<CameraFrame::Ptr>> _camMes;
    std::map<std::string, std::vector<RGBDFrame::Ptr>> _rgbdMes;

    std::map<std::string, ns_veta::Veta::Ptr> _sfmData;

    std::map<std::string, std::vector<OpticalFlowTripleTracePtr>> _visualOpticalFlowTrace;

    double _rawStartTimestamp{};
    double _rawEndTimestamp{};

    double _alignedStartTimestamp{};
    double _alignedEndTimestamp{};

public:
    // using config information to load and adjust data in this constructor
    CalibDataManager();

    // the creator
    static CalibDataManager::Ptr Create();

    // get raw imu measurements
    [[nodiscard]] const std::map<std::string, std::vector<IMUFrame::Ptr>> &GetIMUMeasurements()
        const;

    [[nodiscard]] const std::vector<IMUFrame::Ptr> &GetIMUMeasurements(
        const std::string &imuTopic) const;

    // get raw radar measurements
    [[nodiscard]] const std::map<std::string, std::vector<RadarTargetArray::Ptr>> &
    GetRadarMeasurements() const;

    [[nodiscard]] const std::vector<RadarTargetArray::Ptr> &GetRadarMeasurements(
        const std::string &radarTopic) const;

    // get raw lidar measurements
    [[nodiscard]] const std::map<std::string, std::vector<LiDARFrame::Ptr>> &GetLiDARMeasurements()
        const;

    [[nodiscard]] const std::vector<LiDARFrame::Ptr> &GetLiDARMeasurements(
        const std::string &lidarTopic) const;

    // get raw camera measurements
    // [[nodiscard]] const std::map<std::string, std::vector<CameraFrame::Ptr>> &
    // GetCameraMeasurements() const;

    [[nodiscard]] const std::vector<CameraFrame::Ptr> &GetCameraMeasurements(
        const std::string &camTopic) const;

    // get raw rgbd measurements
    [[nodiscard]] const std::map<std::string, std::vector<RGBDFrame::Ptr>> &GetRGBDMeasurements()
        const;

    [[nodiscard]] const std::vector<RGBDFrame::Ptr> &GetRGBDMeasurements(
        const std::string &rgbdTopic) const;

    // get raw SfM data
    [[nodiscard]] const std::map<std::string, ns_veta::Veta::Ptr> &GetSfMData() const;

    [[nodiscard]] const ns_veta::Veta::Ptr &GetSfMData(const std::string &camTopic) const;

    void SetSfMData(const std::string &camTopic, const ns_veta::Veta::Ptr &veta);

    // [[nodiscard]] const std::map<std::string, std::vector<OpticalFlowTripleTracePtr>> &
    // GetVisualOpticalFlowTrace() const;

    [[nodiscard]] const std::vector<OpticalFlowTripleTracePtr> &GetVisualOpticalFlowTrace(
        const std::string &visualTopic) const;

    void SetVisualOpticalFlowTrace(const std::string &visualTopic,
                                   const std::vector<OpticalFlowTripleTracePtr> &dynamics);

    [[nodiscard]] double GetRawStartTimestamp() const;

    [[nodiscard]] double GetRawEndTimestamp() const;

    [[nodiscard]] double GetAlignedStartTimestamp() const;

    [[nodiscard]] double GetAlignedEndTimestamp() const;

    [[nodiscard]] double GetAlignedTimeRange() const;

    [[nodiscard]] double GetCalibStartTimestamp() const;

    [[nodiscard]] double GetCalibEndTimestamp() const;

    [[nodiscard]] double GetCalibTimeRange() const;

    // sensor measurement frequency

    [[nodiscard]] double GetLiDARAvgFrequency() const;

    [[nodiscard]] double GetCameraAvgFrequency() const;

    [[nodiscard]] double GetRGBDAvgFrequency() const;

    [[nodiscard]] double GetRadarAvgFrequency() const;

    [[nodiscard]] double GetLiDARAvgFrequency(const std::string &topic) const;

    [[nodiscard]] double GetCameraAvgFrequency(const std::string &topic) const;

    [[nodiscard]] double GetRGBDAvgFrequency(const std::string &topic) const;

    [[nodiscard]] double GetRadarAvgFrequency(const std::string &topic) const;

    static auto ExtractIMUDataPiece(const std::vector<IMUFrame::Ptr> &data, double st, double et) {
        auto sIter = std::find_if(data.begin(), data.end(), [st](const IMUFrame::Ptr &frame) {
            return frame->GetTimestamp() > st;
        });
        auto eIter = std::find_if(data.rbegin(), data.rend(), [et](const IMUFrame::Ptr &frame) {
                         return frame->GetTimestamp() < et;
                     }).base();
        return std::pair(sIter, eIter);
    }

    auto ExtractIMUDataPiece(const std::string &topic, double st, double et) {
        return ExtractIMUDataPiece(_imuMes.at(topic), st, et);
    }

    // load camera, lidar, imu data from the ros bag [according to the config file]
    void LoadCalibData();

protected:
    // make sure the first imu frame is before camera and lidar data
    // assign the '_alignedStartTimestamp' and '_alignedEndTimestamp'
    void AdjustCalibDataSequence();

    // align the timestamp to zero
    void AlignTimestamp();

    // remove the head data according to the pred
    template <typename ElemType, typename Pred>
    void EraseSeqHeadData(std::vector<ElemType> &seq, Pred pred, const std::string &errorMsg) {
        auto iter = std::find_if(seq.begin(), seq.end(), pred);
        if (iter == seq.end()) {
            // find failed
            throw Status(Status::ERROR, errorMsg);
        } else {
            // adjust
            seq.erase(seq.begin(), iter);
        }
    }

    // remove the tail data according to the pred
    template <typename ElemType, typename Pred>
    void EraseSeqTailData(std::vector<ElemType> &seq, Pred pred, const std::string &errorMsg) {
        auto iter = std::find_if(seq.rbegin(), seq.rend(), pred);
        if (iter == seq.rend()) {
            // find failed
            throw Status(Status::ERROR, errorMsg);
        } else {
            // adjust
            seq.erase(iter.base(), seq.end());
        }
    }

    // output the data status
    void OutputDataStatus() const;

    template <typename MesSeqType>
    static void CheckTopicExists(const std::string &topic,
                                 const std::map<std::string, MesSeqType> &mesSeq) {
        if (mesSeq.count(topic) == 0) {
            throw Status(Status::CRITICAL,
                         "there is no data in topic '{}'! "
                         "check your configure file and rosbag!",
                         topic);
        }
    }

    template <typename MesType>
    static double GetSensorAvgFrequency(const std::vector<typename MesType::Ptr> &mes) {
        if (mes.empty()) {
            return -1.0;
        } else {
            return static_cast<double>(mes.size()) /
                   (mes.back()->GetTimestamp() - mes.front()->GetTimestamp());
        }
    }

    template <typename MesType>
    static double GetSensorAvgFrequency(
        const std::map<std::string, std::vector<typename MesType::Ptr>> &mesMap) {
        if (mesMap.empty()) {
            return -1.0;
        } else {
            double hz = 0.0;
            for (const auto &[topic, mes] : mesMap) {
                if (double v = GetSensorAvgFrequency<MesType>(mes); v > 0.0) hz += v;
            }
            return hz / static_cast<double>(mesMap.size());
        }
    }
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_CALIB_DATA_MANAGER_H
