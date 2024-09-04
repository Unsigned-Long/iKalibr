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

#include "core/scan_undistortion.h"
#include "calib/calib_param_manager.h"
#include "sensor/lidar.h"
#include "util/tqdm.h"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

ScanUndistortion::ScanUndistortion(const SplineBundleType::Ptr &splines,
                                   CalibParamManager::Ptr calibParamManager)
    : _so3Spline(splines->GetSo3Spline(Configor::Preference::SO3_SPLINE)),
      _posSpline(splines->GetRdSpline(Configor::Preference::SCALE_SPLINE)),
      _parMagr(std::move(calibParamManager)) {}

ScanUndistortion::Ptr ScanUndistortion::Create(const SplineBundleType::Ptr &splines,
                                               const CalibParamManager::Ptr &calibParamManager) {
    return std::make_shared<ScanUndistortion>(splines, calibParamManager);
}

// ---------------
// UndistortToScan
// ---------------

std::vector<LiDARFrame::Ptr> ScanUndistortion::UndistortToScan(
    const std::vector<LiDARFrame::Ptr> &data, const std::string &topic, Option option) {
    std::vector<LiDARFrame::Ptr> lidarUndistFrames;

    auto bar = std::make_shared<tqdm>();
    bool correctPos = IsOptionWith(Option::UNDIST_POS, option);
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        bar->progress(i, static_cast<int>(data.size()));

        if (auto undistLidarFrame = UndistortToScan(data.at(i), topic, correctPos)) {
            lidarUndistFrames.push_back(*undistLidarFrame);
        } else {
            lidarUndistFrames.push_back(nullptr);
        }
    }
    bar->finish();
    return lidarUndistFrames;
}

std::optional<LiDARFrame::Ptr> ScanUndistortion::UndistortToScan(const LiDARFrame::Ptr &lidarFrame,
                                                                 const std::string &topic,
                                                                 bool correctPos) {
    double scanTimeByBr = lidarFrame->GetTimestamp() + _parMagr->TEMPORAL.TO_LkToBr.at(topic);
    // id this time stamp is invalid, return
    if (!_so3Spline.TimeStampInRange(scanTimeByBr) || !_posSpline.TimeStampInRange(scanTimeByBr)) {
        return {};
    }

    Sophus::SE3d scanRefIMUToW(_so3Spline.Evaluate(scanTimeByBr),
                               _posSpline.Evaluate(scanTimeByBr));
    auto scanToRef = scanRefIMUToW * _parMagr->EXTRI.SE3_LkToBr(topic);

    Sophus::SE3d refToScan = scanToRef.inverse();

    // prepare
    IKalibrPointCloud::Ptr undistScan(new IKalibrPointCloud);
    auto rawScan = lidarFrame->GetScan();
    // assign
    undistScan->header = rawScan->header;
    undistScan->height = rawScan->height;
    undistScan->width = rawScan->width;
    undistScan->resize(rawScan->height * rawScan->width);
    undistScan->is_dense = rawScan->is_dense;

    for (int h = 0; h < static_cast<int>(rawScan->height); h++) {
        for (int w = 0; w < static_cast<int>(rawScan->width); w++) {
            const auto &rawPoint = rawScan->points[h * rawScan->width + w];
            IKalibrPoint undistPoint;

            if (IS_POS_NAN(rawPoint)) {
                SET_POS_NAN(undistPoint)
            } else {
                double pTimeByBr = rawPoint.timestamp + _parMagr->TEMPORAL.TO_LkToBr.at(topic);

                if (_so3Spline.TimeStampInRange(pTimeByBr) &&
                    _posSpline.TimeStampInRange(pTimeByBr)) {
                    Sophus::SE3d pRefIMUToW(_so3Spline.Evaluate(pTimeByBr),
                                            _posSpline.Evaluate(pTimeByBr));
                    auto pointToRef = pRefIMUToW * _parMagr->EXTRI.SE3_LkToBr(topic);

                    Sophus::SE3d pointToScan = refToScan * pointToRef;

                    Eigen::Vector3d rp(rawPoint.x, rawPoint.y, rawPoint.z), up;
                    if (correctPos) {
                        up = pointToScan * rp;
                    } else {
                        up = pointToScan.so3() * rp;
                    }

                    undistPoint.x = static_cast<float>(up(0));
                    undistPoint.y = static_cast<float>(up(1));
                    undistPoint.z = static_cast<float>(up(2));
                    undistPoint.timestamp = rawPoint.timestamp;
                    // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
                    // undistPoint.intensity = rawPoint.intensity;
                } else {
                    // we can't undistort it
                    SET_POS_NAN(undistPoint)
                }
            }
            // undistScan->at(w, h) = undistPoint;
            undistScan->points[h * rawScan->width + w] = undistPoint;
        }
    }

    return LiDARFrame::Create(lidarFrame->GetTimestamp(), undistScan);
}

// --------------
// UndistortToRef
// --------------

std::vector<LiDARFrame::Ptr> ScanUndistortion::UndistortToRef(
    const std::vector<LiDARFrame::Ptr> &data, const std::string &topic, Option option) {
    std::vector<LiDARFrame::Ptr> lidarUndistFramesInRef;

    auto bar = std::make_shared<tqdm>();
    bool correctPos = IsOptionWith(Option::UNDIST_POS, option);
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        bar->progress(i, static_cast<int>(data.size()));

        if (auto undistLidarFrame = UndistortToRef(data.at(i), topic, correctPos)) {
            lidarUndistFramesInRef.push_back(*undistLidarFrame);
        } else {
            lidarUndistFramesInRef.push_back(nullptr);
        }
    }
    bar->finish();
    return lidarUndistFramesInRef;
}

std::optional<LiDARFrame::Ptr> ScanUndistortion::UndistortToRef(const LiDARFrame::Ptr &lidarFrame,
                                                                const std::string &topic,
                                                                bool correctPos) {
    double scanTimeByBr = lidarFrame->GetTimestamp() + _parMagr->TEMPORAL.TO_LkToBr.at(topic);
    if (!_so3Spline.TimeStampInRange(scanTimeByBr) || !_posSpline.TimeStampInRange(scanTimeByBr)) {
        return {};
    }
    // prepare
    IKalibrPointCloud::Ptr undistScan(new IKalibrPointCloud);
    auto rawScan = lidarFrame->GetScan();
    // assign
    undistScan->header = rawScan->header;
    undistScan->height = rawScan->height;
    undistScan->width = rawScan->width;
    undistScan->resize(rawScan->height * rawScan->width);
    undistScan->is_dense = rawScan->is_dense;

    for (int h = 0; h < static_cast<int>(rawScan->height); h++) {
        for (int w = 0; w < static_cast<int>(rawScan->width); w++) {
            const auto &rawPoint = rawScan->points[h * rawScan->width + w];
            IKalibrPoint undistPoint;

            if (IS_POS_NAN(rawPoint)) {
                SET_POS_NAN(undistPoint)
            } else {
                double pTimeByBr = rawPoint.timestamp + _parMagr->TEMPORAL.TO_LkToBr.at(topic);

                if (_so3Spline.TimeStampInRange(pTimeByBr) &&
                    _posSpline.TimeStampInRange(pTimeByBr)) {
                    Sophus::SE3d pBrToW(_so3Spline.Evaluate(pTimeByBr),
                                        _posSpline.Evaluate(pTimeByBr));
                    auto pointToW = pBrToW * _parMagr->EXTRI.SE3_LkToBr(topic);

                    Eigen::Vector3d rp(rawPoint.x, rawPoint.y, rawPoint.z), up;
                    if (correctPos) {
                        up = pointToW * rp;
                    } else {
                        up = pointToW.so3() * rp;
                    }

                    undistPoint.x = static_cast<float>(up(0));
                    undistPoint.y = static_cast<float>(up(1));
                    undistPoint.z = static_cast<float>(up(2));
                    undistPoint.timestamp = rawPoint.timestamp;
                    // attention: use 'PointXYZT' as 'IKalibrPoint' rather than 'PointXYZIT' here
                    // undistPoint.intensity = rawPoint.intensity;
                } else {
                    // we can't undistort it (no condition)
                    SET_POS_NAN(undistPoint)
                }
            }
            // undistScan->at(w, h) = undistPoint;
            undistScan->points[h * rawScan->width + w] = undistPoint;
        }
    }

    return LiDARFrame::Create(lidarFrame->GetTimestamp(), undistScan);
}
}  // namespace ns_ikalibr