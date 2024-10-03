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

#ifndef IKALIBR_CALIB_SOLVER_TPL_HPP
#define IKALIBR_CALIB_SOLVER_TPL_HPP

#include "calib_solver.h"
#include "calib/estimator_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

template <TimeDeriv::ScaleSplineType type>
void CalibSolver::AddRadarFactor(Estimator::Ptr &estimator,
                                 const std::string &radarTopic,
                                 Estimator::Opt option) const {
    double weight = Configor::DataStream::RadarTopics.at(radarTopic).Weight;

    for (const auto &targetAry : _dataMagr->GetRadarMeasurements(radarTopic)) {
        for (const auto &tar : targetAry->GetTargets()) {
            estimator->AddRadarMeasurement<type>(tar, radarTopic, option, weight);
        }
    }
}

template <TimeDeriv::ScaleSplineType type>
void CalibSolver::AddAcceFactor(Estimator::Ptr &estimator,
                                const std::string &imuTopic,
                                Estimator::Opt option) const {
    double weight = Configor::DataStream::IMUTopics.at(imuTopic).AcceWeight;

    for (const auto &item : _dataMagr->GetIMUMeasurements(imuTopic)) {
        estimator->AddIMUAcceMeasurement<type>(item, imuTopic, option, weight);
    }
}

template <TimeDeriv::ScaleSplineType type>
void CalibSolver::AddLiDARPointToSurfelFactor(Estimator::Ptr &estimator,
                                              const std::string &lidarTopic,
                                              const std::vector<PointToSurfelCorrPtr> &corrs,
                                              Estimator::Opt option) {
    double weight = Configor::DataStream::LiDARTopics.at(lidarTopic).Weight;

    for (const auto &corr : corrs) {
        estimator->AddLiDARPointTiSurfelConstraint<type>(corr, lidarTopic, option,
                                                         weight * corr->weight);
    }
}

template <TimeDeriv::ScaleSplineType type>
void CalibSolver::AddRGBDPointToSurfelFactor(Estimator::Ptr &estimator,
                                             const std::string &rgbdTopic,
                                             const std::vector<PointToSurfelCorrPtr> &corrs,
                                             Estimator::Opt option) {
    double weight = Configor::DataStream::RGBDTopics.at(rgbdTopic).Weight;

    for (const auto &corr : corrs) {
        estimator->AddRGBDPointTiSurfelConstraint<type>(corr, rgbdTopic, option,
                                                        weight * corr->weight);
    }
}

template <TimeDeriv::ScaleSplineType type>
void CalibSolver::AddVisualReprojectionFactor(Estimator::Ptr &estimator,
                                              const std::string &camTopic,
                                              const std::vector<VisualReProjCorrSeq::Ptr> &corrs,
                                              double *globalScale,
                                              Estimator::Opt option) {
    double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;

    for (const auto &corr : corrs) {
        for (const auto &c : corr->corrs) {
            estimator->AddVisualReprojection<type>(
                c, camTopic, globalScale, corr->invDepthFir.get(), option, weight * c->weight);
        }
    }
}

template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void CalibSolver::AddRGBDOpticalFlowFactor(Estimator::Ptr &estimator,
                                           const std::string &rgbdTopic,
                                           const std::vector<OpticalFlowCorr::Ptr> &corrs,
                                           Estimator::Opt option) {
    double weight = Configor::DataStream::RGBDTopics.at(rgbdTopic).Weight;
    for (const auto &corr : corrs) {
        estimator->AddRGBDOpticalFlowConstraint<type, IsInvDepth>(corr, rgbdTopic, option,
                                                                  weight * corr->weight);
    }
}

template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void CalibSolver::AddVisualOpticalFlowFactor(Estimator::Ptr &estimator,
                                             const std::string &camTopic,
                                             const std::vector<OpticalFlowCorr::Ptr> &corrs,
                                             Estimator::Opt option) {
    double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;
    for (const auto &corr : corrs) {
        estimator->AddVisualOpticalFlowConstraint<type, IsInvDepth>(corr, camTopic, option,
                                                                    weight * corr->weight);
    }
}

template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void CalibSolver::AddVisualOpticalFlowReprojFactor(Estimator::Ptr &estimator,
                                                   const std::string &camTopic,
                                                   const std::vector<OpticalFlowCorr::Ptr> &corrs,
                                                   Estimator::Opt option) {
    if constexpr (type != TimeDeriv::LIN_POS_SPLINE) {
        throw Status(Status::CRITICAL,
                     "the 'AddVisualOpticalFlowReprojConstraint' only works for pos spline!!!");
    }
    double weight = Configor::DataStream::CameraTopics.at(camTopic).Weight;
    for (const auto &corr : corrs) {
        /**
         * given a optical flow tracking correspondence (triple tracking, three points), we throw
         * the middle feature to the camera frame and reproject it to the first and last camera
         * image plane, just like this:
         *                         +-<---<---(*)--->--->-+
         *                         |          ^          |
         *                         v          |          v
         *                      [ fir        mid        last ] -> a optical flow tracking (triple)
         */
        estimator->AddVisualOpticalFlowReprojConstraint<type, IsInvDepth>(corr, camTopic, option,
                                                                          weight * corr->weight);
    }
}

template <TimeDeriv::ScaleSplineType type, bool IsInvDepth>
void CalibSolver::AddRGBDOpticalFlowReprojFactor(Estimator::Ptr &estimator,
                                                 const std::string &camTopic,
                                                 const std::vector<OpticalFlowCorr::Ptr> &corrs,
                                                 Estimator::Opt option) {
    if constexpr (type != TimeDeriv::LIN_POS_SPLINE) {
        throw Status(Status::CRITICAL,
                     "the 'AddRGBDOpticalFlowReprojConstraint' only works for pos spline!!!");
    }
    double weight = Configor::DataStream::RGBDTopics.at(camTopic).Weight;
    for (const auto &corr : corrs) {
        /**
         * given a optical flow tracking correspondence (triple tracking, three points), we throw
         * the middle feature to the camera frame and reproject it to the first and last camera
         * image plane, just like this:
         *                         --<---<---(*)--->--->--
         *                         |          ^          |
         *                         v          |          v
         *                      [ fir        mid        last ] -> a optical flow tracking (triple)
         */
        estimator->AddRGBDOpticalFlowReprojConstraint<type, IsInvDepth>(corr, camTopic, option,
                                                                        weight * corr->weight);
    }
}
}  // namespace ns_ikalibr

#endif  // IKALIBR_CALIB_SOLVER_TPL_HPP
