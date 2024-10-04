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

#include "calib/calib_data_manager.h"
#include "core/optical_flow_trace.h"
#include "opencv4/opencv2/imgcodecs.hpp"
#include "rosbag/view.h"
#include "sensor/camera_data_loader.h"
#include "sensor/depth_data_loader.h"
#include "sensor/imu_data_loader.h"
#include "sensor/lidar_data_loader.h"
#include "sensor/radar_data_loader.h"
#include "spdlog/spdlog.h"
#include "util/tqdm.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

// ----------------
// CalibDataManager
// ----------------

CalibDataManager::CalibDataManager() = default;

CalibDataManager::Ptr CalibDataManager::Create() { return std::make_shared<CalibDataManager>(); }

void CalibDataManager::LoadCalibData() {
    spdlog::info("loading calibration data...");

    // open the ros bag
    auto bag = std::make_unique<rosbag::Bag>();
    if (!std::filesystem::exists(Configor::DataStream::BagPath)) {
        spdlog::error("the ros bag path '{}' is invalid!", Configor::DataStream::BagPath);
    } else {
        bag->open(Configor::DataStream::BagPath, rosbag::BagMode::Read);
    }

    auto view = rosbag::View();

    // using a temp view to check the time range of the source ros bag
    auto viewTemp = rosbag::View();

    std::vector<std::string> topicsToQuery;
    // add topics to vector
    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        topicsToQuery.push_back(topic);
    }
    for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
        topicsToQuery.push_back(topic);
    }
    for (const auto &[topic, _] : Configor::DataStream::LiDARTopics) {
        topicsToQuery.push_back(topic);
    }
    for (const auto &[topic, _] : Configor::DataStream::CameraTopics) {
        topicsToQuery.push_back(topic);
    }
    for (const auto &[topic, info] : Configor::DataStream::RGBDTopics) {
        topicsToQuery.push_back(topic);
        topicsToQuery.push_back(info.DepthTopic);
    }

    viewTemp.addQuery(*bag, rosbag::TopicQuery(topicsToQuery));
    auto begTime = viewTemp.getBeginTime();
    auto endTime = viewTemp.getEndTime();
    spdlog::info("source data duration: from '{:.5f}' to '{:.5f}'.", begTime.toSec(),
                 endTime.toSec());

    // adjust the data time range
    if (Configor::DataStream::BeginTime > 0.0) {
        begTime += ros::Duration(Configor::DataStream::BeginTime);
        if (begTime > endTime) {
            spdlog::warn(
                "begin time '{:.5f}' is out of the bag's data range, set begin time to '{:.5f}'.",
                begTime.toSec(), viewTemp.getBeginTime().toSec());
            begTime = viewTemp.getBeginTime();
        }
    }
    if (Configor::DataStream::Duration > 0.0) {
        endTime = begTime + ros::Duration(Configor::DataStream::Duration);
        if (endTime > viewTemp.getEndTime()) {
            spdlog::warn(
                "end time '{:.5f}' is out of the bag's data range, set end time to '{:.5f}'.",
                endTime.toSec(), viewTemp.getEndTime().toSec());
            endTime = viewTemp.getEndTime();
        }
    }
    spdlog::info("expect data duration: from '{:.5f}' to '{:.5f}'.", begTime.toSec(),
                 endTime.toSec());

    view.addQuery(*bag, rosbag::TopicQuery(topicsToQuery), begTime, endTime);

    // create IMU data loader
    std::map<std::string, IMUDataLoader::Ptr> imuDataLoaders;
    std::map<std::string, RadarDataLoader::Ptr> radarDataLoaders;
    std::map<std::string, LiDARDataLoader::Ptr> lidarDataLoaders;
    std::map<std::string, CameraDataLoader::Ptr> cameraDataLoaders;
    // for rgbd cameras
    std::map<std::string, CameraDataLoader::Ptr> rgbdColorDataLoaders;
    std::map<std::string, DepthDataLoader::Ptr> rgbdDepthDataLoaders;
    // temporal data containers ('list' containers)
    std::map<std::string, std::list<CameraFrame::Ptr>> rgbdColorMesTemp;
    std::map<std::string, std::list<DepthFrame::Ptr>> rgbdDepthMesTemp;

    // get type enum from the string
    for (const auto &[topic, config] : Configor::DataStream::IMUTopics) {
        imuDataLoaders.insert({topic, IMUDataLoader::GetLoader(config.Type)});
    }
    for (const auto &[topic, config] : Configor::DataStream::RadarTopics) {
        radarDataLoaders.insert({topic, RadarDataLoader::GetLoader(config.Type)});
    }
    for (const auto &[topic, config] : Configor::DataStream::LiDARTopics) {
        lidarDataLoaders.insert({topic, LiDARDataLoader::GetLoader(config.Type)});
    }
    for (const auto &[topic, config] : Configor::DataStream::CameraTopics) {
        cameraDataLoaders.insert({topic, CameraDataLoader::GetLoader(config.Type)});
    }
    for (const auto &[topic, config] : Configor::DataStream::RGBDTopics) {
        rgbdColorDataLoaders.insert({topic, CameraDataLoader::GetLoader(config.Type)});
        bool isInverse = config.DepthFactor < 0.0f;
        rgbdDepthDataLoaders.insert(
            {config.DepthTopic, DepthDataLoader::GetLoader(config.Type, isInverse)});
    }

    // read raw data
    auto bar = std::make_shared<tqdm>();
    int idx = 0;
    for (auto iter = view.begin(); iter != view.end(); ++iter, ++idx) {
        bar->progress(idx, static_cast<int>(view.size()));
        const auto &item = *iter;
        const std::string &topic = item.getTopic();
        if (imuDataLoaders.cend() != imuDataLoaders.find(topic)) {
            // is an inertial frame
            auto mes = imuDataLoaders.at(topic)->UnpackFrame(item);
            if (mes != nullptr) {
                _imuMes[topic].push_back(mes);
            }
        } else if (radarDataLoaders.cend() != radarDataLoaders.find(topic)) {
            // is a radar frame
            auto mes = radarDataLoaders.at(topic)->UnpackScan(item);
            if (mes != nullptr) {
                _radarMes[topic].push_back(mes);
            }
        } else if (lidarDataLoaders.cend() != lidarDataLoaders.find(topic)) {
            // is a lidar frame
            auto mes = lidarDataLoaders.at(topic)->UnpackScan(item);
            if (mes != nullptr) {
                _lidarMes[topic].push_back(mes);
            }
        } else if (rgbdColorDataLoaders.cend() != rgbdColorDataLoaders.find(topic)) {
            // is a rgbd color frame
            auto mes = rgbdColorDataLoaders.at(topic)->UnpackFrame(item);
            if (mes != nullptr) {
                // id: uint64_t from timestamp (raw, millisecond)
                mes->SetId(static_cast<ns_veta::IndexT>(mes->GetTimestamp() * 1E3));
                rgbdColorMesTemp[topic].push_back(mes);
            }
        } else if (rgbdDepthDataLoaders.cend() != rgbdDepthDataLoaders.find(topic)) {
            // is a rgbd depth frame
            auto mes = rgbdDepthDataLoaders.at(topic)->UnpackFrame(item);
            if (mes != nullptr) {
                // id: uint64_t from timestamp (raw, millisecond)
                mes->SetId(static_cast<ns_veta::IndexT>(mes->GetTimestamp() * 1E3));
                rgbdDepthMesTemp[topic].push_back(mes);
            }
        } else if (cameraDataLoaders.cend() != cameraDataLoaders.find(topic)) {
            // is a camera frame
            auto mes = cameraDataLoaders.at(topic)->UnpackFrame(item);
            if (mes != nullptr) {
                // id: uint64_t from timestamp (raw, millisecond)
                mes->SetId(static_cast<ns_veta::IndexT>(mes->GetTimestamp() * 1E3));
                _camMes[topic].push_back(mes);
            }
        }
    }
    bar->finish();
    bag->close();

    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        CheckTopicExists(topic, _imuMes);
    }
    for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
        CheckTopicExists(topic, _radarMes);
    }
    for (const auto &[topic, _] : Configor::DataStream::LiDARTopics) {
        CheckTopicExists(topic, _lidarMes);
    }
    for (const auto &[topic, _] : Configor::DataStream::CameraTopics) {
        CheckTopicExists(topic, _camMes);
    }
    for (const auto &[topic, info] : Configor::DataStream::RGBDTopics) {
        // measurements are stored in 'rgbdColorMesTemp' and 'rgbdDepthMesTemp' temporally
        CheckTopicExists(topic, rgbdColorMesTemp);
        CheckTopicExists(info.DepthTopic, rgbdDepthMesTemp);
    }

    // if the radar is AWR1843BOOST, data should be reorganized,
    // i.e., merge multiple radar target measurements to radar array measurements
    // note that although radar targets are wrapped as scans here (just for unification and
    // convenience), they are still fused separately in batch optimizations (a tightly-coupled
    // optimization framework)
    std::map<std::string, std::vector<RadarTargetArray::Ptr>> oldRadarMes = _radarMes;
    for (const auto &[topic, loader] : radarDataLoaders) {
        if (loader->GetRadarModel() == RadarModelType::AWR1843BOOST_RAW ||
            loader->GetRadarModel() == RadarModelType::AWR1843BOOST_CUSTOM) {
            const auto &mes = oldRadarMes.at(topic);
            std::vector<RadarTarget::Ptr> targets;
            std::vector<RadarTargetArray::Ptr> arrays;
            for (const auto &item : mes) {
                // merge measurements by 10 HZ (0.1 s)
                if (targets.empty() ||
                    std::abs(targets.front()->GetTimestamp() - item->GetTimestamp()) < 0.1) {
                    targets.push_back(item->GetTargets().front());
                } else {
                    // compute average time as the timestamp of radar target array
                    double t = 0.0;
                    for (const auto &target : targets) {
                        t += target->GetTimestamp() / static_cast<double>(targets.size());
                    }
                    arrays.push_back(RadarTargetArray::Create(t, targets));
                    targets.clear();
                    targets.push_back(item->GetTargets().front());
                }
            }
            _radarMes.at(topic) = arrays;
        }
    }

    // match color and depth images for rgbd cameras
    // 'rgbdColorMesTemp' + 'rgbdDepthMesTemp' -->  '_rgbdMes'
    for (const auto &[colorTopic, colorFrames] : rgbdColorMesTemp) {
        auto depthTopic = Configor::DataStream::RGBDTopics.at(colorTopic).DepthTopic;
        auto &curRGBDMes = _rgbdMes[colorTopic];
        auto &depthImgs = rgbdDepthMesTemp.at(depthTopic);
        for (const auto &colorFrame : colorFrames) {
            const double timestamp = colorFrame->GetTimestamp();
            // find a matched depth image for current color image
            // the matched depth and color images should be close enough to each other in time
            // domain
            auto iter = std::find_if(
                depthImgs.begin(), depthImgs.end(), [timestamp](const DepthFrame::Ptr &depthFrame) {
                    // time distance: 0.001 s, i.e., 1 ms
                    return std::abs(depthFrame->GetTimestamp() - timestamp) < 1E-3;
                });
            if (iter != depthImgs.cend()) {
                curRGBDMes.push_back(RGBDFrame::Create(timestamp,                    // timestamp
                                                       colorFrame->GetImage(),       // grey image
                                                       colorFrame->GetColorImage(),  // color image
                                                       (*iter)->GetDepthImage(),     // depth image
                                                       colorFrame->GetId()           // image index
                                                       ));
                // remove this depth image
                depthImgs.erase(iter);
            } else {
                spdlog::warn(
                    "can not find a matched depth image from '{}' for color image from '{}' at "
                    "time '{:.5f}'",
                    depthTopic, colorTopic, timestamp);
            }
        }
    }
    rgbdColorMesTemp.clear(), rgbdDepthMesTemp.clear();
    for (const auto &[topic, info] : Configor::DataStream::RGBDTopics) {
        CheckTopicExists(topic, _rgbdMes);
    }

    OutputDataStatus();

    AdjustCalibDataSequence();
    AlignTimestamp();

    /**
     * to calibrate velocity-spline-derived cameras, high sampling frequency is required (larger
     * than 30 Hz), to perform high-precision optical flow velocity recovery
     */
    for (const auto &[topic, _] : Configor::DataStream::VelCameraTopics()) {
        auto freq = GetCameraAvgFrequency(topic);
        spdlog::info("sampling frequency for camera '{}': {:.3f}", topic, freq);
        if (freq < 29.0) {
            throw Status(
                Status::WARNING,
                "Sampling frequency of vel camera '{}' (freq: {:.3f}) is too small!!! "
                "Frequency larger than 30 Hz is required!!! Please change 'ScaleSplineType' of "
                "this camera to 'LIN_POS_SPLINE' which would perform a SfM-based calibration!!! Do "
                "not forget to change its weight!",
                topic, freq);
        }
    }

    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        auto freq = GetRGBDAvgFrequency(topic);
        spdlog::info("sampling frequency for rgbd camera '{}': {:.3f}", topic, freq);
        if (freq < 29.0) {
            throw Status(Status::WARNING,
                         "Sampling frequency of rgbd camera '{}' (freq: {:.3f}) is too small!!! "
                         "Frequency larger than 30 Hz is required!!! Please throw the depth "
                         "information and treat it an optical camera, and perform "
                         "'LIN_POS_SPLINE'-based calibration.",
                         topic, freq);
        }
    }
}

void CalibDataManager::AdjustCalibDataSequence() {
    spdlog::info("adjust calibration data sequence...");

    // data sequence adjustment pattern (first step)
    /**
     *       |erased|                           |erased |
     * --------------------------------------------------
     * IMU1: |o o o |o o o o o o o o o o o o o o|o      |
     * IMU2: |      |o o o o o o o o o o o o o o|o o o  |
     * RAD1: |   o o|o o o o o o o o o o o o o o|       |
     * RAD2: | o o o|o o o o o o o o o o o o o o|o o o o|
     * CAM1: | o o o|o o o o o o o o o o o o o o|o      |
     * CAM2: | o o o|o o o o o o o o o o o o o o|o o    |
     * LID1: | o o o|o o o o o o o o o o o o o o|o      |
     * LID2: |   o o|o o o o o o o o o o o o o o|o      |
     * --------------------------------------------------
     *              |--> imuMinTime             |--> imuMaxTime
     */
    auto imuMinTime = std::max_element(_imuMes.begin(), _imuMes.end(),
                                       [](const auto &p1, const auto &p2) {
                                           return p1.second.front()->GetTimestamp() <
                                                  p2.second.front()->GetTimestamp();
                                       })
                          ->second.front()
                          ->GetTimestamp();
    auto imuMaxTime = std::min_element(_imuMes.begin(), _imuMes.end(),
                                       [](const auto &p1, const auto &p2) {
                                           return p1.second.back()->GetTimestamp() <
                                                  p2.second.back()->GetTimestamp();
                                       })
                          ->second.back()
                          ->GetTimestamp();

    _rawStartTimestamp = imuMinTime;
    _rawEndTimestamp = imuMaxTime;

    // data sequence adjustment pattern (second step)
    /**
     *              |--> imuMinTime             |--> imuMaxTime
     * --------------------------------------------------
     * IMU1:        |o|o o o o o o o o o o o o|o|
     * IMU2:        |o|o o o o o o o o o o o o|o|
     * RAD1:        |o|o o o o o o o o o o o o|o|
     * RAD2:        |o|o o o o o o o o o o o o|o|
     * CAM1:        | | |o o o o o o o o o o| | |
     * CAM2:        | | |o o o o o o o o o o| | |
     * LID1:        | | |o o o o o o o o o o| | |
     * LID2:        | | |o o o o o o o o o o| | |
     * --------------------------------------------------
     *                |--> calibSTime         |--> calibETime
     */

    if (Configor::IsRadarIntegrated()) {
        auto radarMinTime = std::max_element(_radarMes.begin(), _radarMes.end(),
                                             [](const auto &p1, const auto &p2) {
                                                 return p1.second.front()->GetTimestamp() <
                                                        p2.second.front()->GetTimestamp();
                                             })
                                ->second.front()
                                ->GetTimestamp();
        auto radarMaxTime = std::min_element(_radarMes.begin(), _radarMes.end(),
                                             [](const auto &p1, const auto &p2) {
                                                 return p1.second.back()->GetTimestamp() <
                                                        p2.second.back()->GetTimestamp();
                                             })
                                ->second.back()
                                ->GetTimestamp();
        _rawStartTimestamp = std::max(_rawStartTimestamp, radarMinTime);
        _rawEndTimestamp = std::min(_rawEndTimestamp, radarMaxTime);
    }

    if (Configor::IsLiDARIntegrated()) {
        auto lidarMinTime = std::max_element(_lidarMes.begin(), _lidarMes.end(),
                                             [](const auto &p1, const auto &p2) {
                                                 return p1.second.front()->GetTimestamp() <
                                                        p2.second.front()->GetTimestamp();
                                             })
                                ->second.front()
                                ->GetTimestamp();
        auto lidarMaxTime = std::min_element(_lidarMes.begin(), _lidarMes.end(),
                                             [](const auto &p1, const auto &p2) {
                                                 return p1.second.back()->GetTimestamp() <
                                                        p2.second.back()->GetTimestamp();
                                             })
                                ->second.back()
                                ->GetTimestamp();
        _rawStartTimestamp = std::max({_rawStartTimestamp, lidarMinTime});
        _rawEndTimestamp = std::min({_rawEndTimestamp, lidarMaxTime});
    }

    if (Configor::IsPosCameraIntegrated() || Configor::IsVelCameraIntegrated()) {
        auto camMinTime = std::max_element(_camMes.begin(), _camMes.end(),
                                           [](const auto &p1, const auto &p2) {
                                               return p1.second.front()->GetTimestamp() <
                                                      p2.second.front()->GetTimestamp();
                                           })
                              ->second.front()
                              ->GetTimestamp();
        auto camMaxTime = std::min_element(_camMes.begin(), _camMes.end(),
                                           [](const auto &p1, const auto &p2) {
                                               return p1.second.back()->GetTimestamp() <
                                                      p2.second.back()->GetTimestamp();
                                           })
                              ->second.back()
                              ->GetTimestamp();
        _rawStartTimestamp = std::max({_rawStartTimestamp, camMinTime});
        _rawEndTimestamp = std::min({_rawEndTimestamp, camMaxTime});
    }

    if (Configor::IsRGBDIntegrated()) {
        auto rgbdMinTime = std::max_element(_rgbdMes.begin(), _rgbdMes.end(),
                                            [](const auto &p1, const auto &p2) {
                                                return p1.second.front()->GetTimestamp() <
                                                       p2.second.front()->GetTimestamp();
                                            })
                               ->second.front()
                               ->GetTimestamp();
        auto rgbdMaxTime = std::min_element(_rgbdMes.begin(), _rgbdMes.end(),
                                            [](const auto &p1, const auto &p2) {
                                                return p1.second.back()->GetTimestamp() <
                                                       p2.second.back()->GetTimestamp();
                                            })
                               ->second.back()
                               ->GetTimestamp();
        _rawStartTimestamp = std::max({_rawStartTimestamp, rgbdMinTime});
        _rawEndTimestamp = std::min({_rawEndTimestamp, rgbdMaxTime});
    }

    for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
        // remove imu frames that are before the start time stamp
        EraseSeqHeadData(
            _imuMes.at(topic),
            [this](const IMUFrame::Ptr &frame) {
                return frame->GetTimestamp() > _rawStartTimestamp;
            },
            "the imu data is invalid, there is no intersection.");

        // remove imu frames that are after the end time stamp
        EraseSeqTailData(
            _imuMes.at(topic),
            [this](const IMUFrame::Ptr &frame) { return frame->GetTimestamp() < _rawEndTimestamp; },
            "the imu data is invalid, there is no intersection.");
    }

    for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
        // remove radar frames that are before the start time stamp
        EraseSeqHeadData(
            _radarMes.at(topic),
            [this](const RadarTargetArray::Ptr &frame) {
                return frame->GetTimestamp() >
                       _rawStartTimestamp + 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the radar data is invalid, there is no intersection between imu data and radar data.");

        // remove radar frames that are after the end time stamp
        EraseSeqTailData(
            _radarMes.at(topic),
            [this](const RadarTargetArray::Ptr &frame) {
                return frame->GetTimestamp() <
                       _rawEndTimestamp - 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the radar data is invalid, there is no intersection between imu data and radar data.");
    }

    for (const auto &[topic, _] : Configor::DataStream::LiDARTopics) {
        // remove lidar frames that are before the start time stamp
        EraseSeqHeadData(
            _lidarMes.at(topic),
            [this](const LiDARFrame::Ptr &frame) {
                // different from other sensor, a time offset padding is used here
                // to ensure that the time of first lidar frame (map time) is in spline time range
                return frame->GetTimestamp() >
                       _rawStartTimestamp + 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the lidar data is invalid, there is no intersection between imu data and lidar data.");

        // remove lidar frames that are after the end time stamp
        EraseSeqTailData(
            _lidarMes.at(topic),
            [this](const LiDARFrame::Ptr &frame) {
                return frame->GetTimestamp() <
                       _rawEndTimestamp - 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the lidar data is invalid, there is no intersection between imu data and lidar data.");
    }

    for (const auto &[topic, _] : Configor::DataStream::CameraTopics) {
        // remove camera frames that are before the start time stamp
        EraseSeqHeadData(
            _camMes.at(topic),
            [this](const CameraFrame::Ptr &frame) {
                // different from other sensor, a time offset padding is used here
                // to ensure that the time of first lidar frame (map time) is in spline time range
                return frame->GetTimestamp() >
                       _rawStartTimestamp + 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the camera data is invalid, there is no intersection between imu data and lidar "
            "data.");

        // remove camera frames that are after the end time stamp
        EraseSeqTailData(
            _camMes.at(topic),
            [this](const CameraFrame::Ptr &frame) {
                return frame->GetTimestamp() <
                       _rawEndTimestamp - 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the camera data is invalid, there is no intersection between imu data and lidar "
            "data.");
    }

    for (const auto &[topic, _] : Configor::DataStream::RGBDTopics) {
        // remove rgbd frames that are before the start time stamp
        EraseSeqHeadData(
            _rgbdMes.at(topic),
            [this](const RGBDFrame::Ptr &frame) {
                return frame->GetTimestamp() >
                       _rawStartTimestamp + 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the rgbd data is invalid, there is no intersection between imu data and rgbd data.");

        // remove rgbd frames that are after the end time stamp
        EraseSeqTailData(
            _rgbdMes.at(topic),
            [this](const RGBDFrame::Ptr &frame) {
                return frame->GetTimestamp() <
                       _rawEndTimestamp - 2 * Configor::Prior::TimeOffsetPadding;
            },
            "the rgbd data is invalid, there is no intersection between imu data and rgbd data.");
    }
    OutputDataStatus();
}

void CalibDataManager::AlignTimestamp() {
    spdlog::info("align calibration data timestamp...");

    // all time stamps minus  '_rawStartTimestamp'
    _alignedStartTimestamp = 0.0;
    _alignedEndTimestamp = _rawEndTimestamp - _rawStartTimestamp;
    for (auto &[imuTopic, mes] : _imuMes) {
        for (const auto &frame : mes) {
            frame->SetTimestamp(frame->GetTimestamp() - _rawStartTimestamp);
        }
    }
    for (const auto &[radarTopic, mes] : _radarMes) {
        for (const auto &array : mes) {
            // array
            array->SetTimestamp(array->GetTimestamp() - _rawStartTimestamp);
            // targets
            for (auto &tar : array->GetTargets()) {
                tar->SetTimestamp(tar->GetTimestamp() - _rawStartTimestamp);
            }
        }
    }
    for (const auto &[topic, data] : _lidarMes) {
        for (const auto &item : data) {
            item->SetTimestamp(item->GetTimestamp() - _rawStartTimestamp);
            for (auto &p : *item->GetScan()) {
                p.timestamp -= _rawStartTimestamp;
            }
        }
    }
    for (auto &[camTopic, mes] : _camMes) {
        for (const auto &frame : mes) {
            frame->SetTimestamp(frame->GetTimestamp() - _rawStartTimestamp);
        }
    }
    for (auto &[rgbdTopic, mes] : _rgbdMes) {
        for (const auto &frame : mes) {
            frame->SetTimestamp(frame->GetTimestamp() - _rawStartTimestamp);
        }
    }
    OutputDataStatus();
}

void CalibDataManager::OutputDataStatus() const {
    spdlog::info("calibration data info:");
    for (const auto &[topic, mes] : _imuMes) {
        spdlog::info(
            "IMU topic: '{}', data size: '{:06}', time span: from '{:+010.5f}' to '{:+010.5f}' (s)",
            topic, mes.size(), mes.front()->GetTimestamp(), mes.back()->GetTimestamp());
    }
    for (const auto &[topic, mes] : _radarMes) {
        spdlog::info(
            "Radar topic: '{}', data size: '{:06}', time span: from '{:+010.5f}' to '{:+010.5f}' "
            "(s)",
            topic, mes.size(), mes.front()->GetTimestamp(), mes.back()->GetTimestamp());
    }
    for (const auto &[topic, mes] : _lidarMes) {
        spdlog::info(
            "LiDAR topic: '{}', data size: '{:06}', time span: from '{:+010.5f}' to '{:+010.5f}' "
            "(s)",
            topic, mes.size(), mes.front()->GetTimestamp(), mes.back()->GetTimestamp());
    }
    for (const auto &[topic, mes] : _camMes) {
        spdlog::info(
            "Camera topic: '{}', data size: '{:06}', time span: from '{:+010.5f}' to '{:+010.5f}' "
            "(s)",
            topic, mes.size(), mes.front()->GetTimestamp(), mes.back()->GetTimestamp());
    }
    for (const auto &[topic, mes] : _rgbdMes) {
        spdlog::info(
            "RGBD topic: '{}', data size: '{:06}', time span: from '{:+010.5f}' to '{:+010.5f}' "
            "(s)",
            topic, mes.size(), mes.front()->GetTimestamp(), mes.back()->GetTimestamp());
    }

    spdlog::info("raw start time: '{:+010.5f}' (s), raw end time: '{:+010.5f}' (s)",
                 GetRawStartTimestamp(), GetRawEndTimestamp());
    spdlog::info("aligned start time: '{:+010.5f}' (s), aligned end time: '{:+010.5f}' (s)",
                 GetAlignedStartTimestamp(), GetAlignedEndTimestamp());
    spdlog::info("calib start time: '{:+010.5f}' (s), calib end time: '{:+010.5f}' (s)\n",
                 GetCalibStartTimestamp(), GetCalibEndTimestamp());
}

// -----------
// time access
// -----------

double CalibDataManager::GetRawStartTimestamp() const { return _rawStartTimestamp; }

double CalibDataManager::GetRawEndTimestamp() const { return _rawEndTimestamp; }

double CalibDataManager::GetAlignedStartTimestamp() const {
    // Function 'GetAlignedStartTimestamp' always returns 0.0
    return _alignedStartTimestamp;
}

double CalibDataManager::GetAlignedEndTimestamp() const { return _alignedEndTimestamp; }

double CalibDataManager::GetAlignedTimeRange() const {
    return _alignedEndTimestamp - _alignedStartTimestamp;
}

double CalibDataManager::GetCalibStartTimestamp() const {
    return _alignedStartTimestamp + Configor::Prior::TimeOffsetPadding;
}

double CalibDataManager::GetCalibEndTimestamp() const {
    return _alignedEndTimestamp - Configor::Prior::TimeOffsetPadding;
}

double CalibDataManager::GetCalibTimeRange() const {
    return GetCalibEndTimestamp() - GetCalibStartTimestamp();
}

double CalibDataManager::GetLiDARAvgFrequency() const {
    return GetSensorAvgFrequency<LiDARFrame>(_lidarMes);
}

double CalibDataManager::GetCameraAvgFrequency() const {
    return GetSensorAvgFrequency<CameraFrame>(_camMes);
}

double CalibDataManager::GetRGBDAvgFrequency() const {
    return GetSensorAvgFrequency<RGBDFrame>(_rgbdMes);
}

double CalibDataManager::GetRadarAvgFrequency() const {
    return GetSensorAvgFrequency<RadarTargetArray>(_radarMes);
}

double CalibDataManager::GetLiDARAvgFrequency(const std::string &topic) const {
    return GetSensorAvgFrequency<LiDARFrame>(_lidarMes.at(topic));
}

double CalibDataManager::GetCameraAvgFrequency(const std::string &topic) const {
    return GetSensorAvgFrequency<CameraFrame>(_camMes.at(topic));
}

double CalibDataManager::GetRGBDAvgFrequency(const std::string &topic) const {
    return GetSensorAvgFrequency<RGBDFrame>(_rgbdMes.at(topic));
}

double CalibDataManager::GetRadarAvgFrequency(const std::string &topic) const {
    return GetSensorAvgFrequency<RadarTargetArray>(_radarMes.at(topic));
}

// -----------
// data access
// -----------

const std::map<std::string, std::vector<IMUFrame::Ptr>> &CalibDataManager::GetIMUMeasurements()
    const {
    return _imuMes;
}

const std::vector<IMUFrame::Ptr> &CalibDataManager::GetIMUMeasurements(
    const std::string &imuTopic) const {
    return _imuMes.at(imuTopic);
}

const std::map<std::string, std::vector<RadarTargetArray::Ptr>> &
CalibDataManager::GetRadarMeasurements() const {
    return _radarMes;
}

const std::vector<RadarTargetArray::Ptr> &CalibDataManager::GetRadarMeasurements(
    const std::string &radarTopic) const {
    return _radarMes.at(radarTopic);
}

const std::map<std::string, std::vector<LiDARFrame::Ptr>> &CalibDataManager::GetLiDARMeasurements()
    const {
    return _lidarMes;
}

const std::vector<LiDARFrame::Ptr> &CalibDataManager::GetLiDARMeasurements(
    const std::string &lidarTopic) const {
    return _lidarMes.at(lidarTopic);
}

// const std::map<std::string, std::vector<CameraFrame::Ptr>> &
// CalibDataManager::GetCameraMeasurements() const {
//     return _camMes;
// }

const std::vector<CameraFrame::Ptr> &CalibDataManager::GetCameraMeasurements(
    const std::string &camTopic) const {
    return _camMes.at(camTopic);
}

const std::map<std::string, std::vector<RGBDFrame::Ptr>> &CalibDataManager::GetRGBDMeasurements()
    const {
    return _rgbdMes;
}

const std::vector<RGBDFrame::Ptr> &CalibDataManager::GetRGBDMeasurements(
    const std::string &rgbdTopic) const {
    return _rgbdMes.at(rgbdTopic);
}

const std::map<std::string, ns_veta::Veta::Ptr> &CalibDataManager::GetSfMData() const {
    return _sfmData;
}

const ns_veta::Veta::Ptr &CalibDataManager::GetSfMData(const std::string &camTopic) const {
    return _sfmData.at(camTopic);
}

void CalibDataManager::SetSfMData(const std::string &camTopic, const ns_veta::Veta::Ptr &veta) {
    _sfmData[camTopic] = veta;
}
void CalibDataManager::SetVisualOpticalFlowTrace(
    const std::string &visualTopic, const std::vector<OpticalFlowTripleTrace::Ptr> &dynamics) {
    _visualOpticalFlowTrace[visualTopic] = dynamics;
}

// const std::map<std::string, std::vector<OpticalFlowTripleTrace::Ptr>> &
// CalibDataManager::GetVisualOpticalFlowTrace() const {
//     return _visualOpticalFlowTrace;
// }

const std::vector<OpticalFlowTripleTrace::Ptr> &CalibDataManager::GetVisualOpticalFlowTrace(
    const std::string &visualTopic) const {
    return _visualOpticalFlowTrace.at(visualTopic);
}
}  // namespace ns_ikalibr