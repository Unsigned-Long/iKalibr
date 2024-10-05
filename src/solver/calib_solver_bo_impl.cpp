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

#include "solver/calib_solver_tpl.hpp"
#include "magic_enum_flags.hpp"
#include "util/utils_tpl.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

CalibSolver::BackUp::Ptr CalibSolver::BatchOptimization(
    OptOption optOption,
    const std::map<std::string, std::vector<PointToSurfelCorr::Ptr>> &lidarPtsCorrs,
    const std::map<std::string, std::vector<VisualReProjCorrSeq::Ptr>> &visualReprojCorrs,
    const std::map<std::string, std::vector<OpticalFlowCorr::Ptr>> &rgbdCorrs,
    const std::map<std::string, std::vector<OpticalFlowCorr::Ptr>> &visualVelCorrs,
    const std::optional<std::map<std::string, std::vector<PointToSurfelCorrPtr>>> &rgbdPtsCorrs)
    const {
    // a lambda function to obtain the string of current optimization option
    auto GetOptString = [](OptOption opt) -> std::string {
        std::stringstream stringStream;
        stringStream << magic_enum::enum_flags_name(opt);
        return stringStream.str();
    };
    // a lambda function to refine optimization option for gs/rs cameras
    auto RefineReadoutTimeOptForCameras = [](const std::string &topic,
                                             const OptOption opt) -> OptOption {
        OptOption visualOpt = opt;
        /*
         * if this camera is a rolling shutter camera and optimization option is with
         * 'OPT_RS_CAM_READOUT_TIME' option, we optimize the rs readout time.
         * if this camera is a global shutter camera and optimization option is with
         * 'OPT_RS_CAM_READOUT_TIME' option, we remove this option and do not optimize the
         * rs readout time.
         */
        if (IsOptionWith(OptOption::OPT_RS_CAM_READOUT_TIME, visualOpt)) {
            if (IsRSCamera(topic)) {
                spdlog::info(
                    "camera '{}' is a rolling shutter (RS) camera, use optimization option "
                    "'OPT_RS_CAM_READOUT_TIME'",
                    topic);
            } else {
                visualOpt ^= OptOption::OPT_RS_CAM_READOUT_TIME;
                spdlog::info(
                    "camera '{}' is a global shutter (GS) camera, remove optimization option "
                    "'OPT_RS_CAM_READOUT_TIME'",
                    topic);
            }
        }
        return visualOpt;
    };

    spdlog::info("Optimization option: {}", GetOptString(optOption));

    auto estimator = Estimator::Create(_splines, _parMagr);
    auto visualGlobalScale = std::make_shared<double>(1.0);
    constexpr bool OPTICAL_FLOW_EST_INV_DEPTH = true;

    switch (GetScaleType()) {
        case TimeDeriv::LIN_ACCE_SPLINE: {
            /**
             * only when imu-only multi-imu calibration is required, the linear acceleration spline
             * would be maintained
             */
            for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
                this->AddAcceFactor<TimeDeriv::LIN_ACCE_SPLINE>(estimator, topic, optOption);
                this->AddGyroFactor(estimator, topic, optOption);
            }
        } break;
        case TimeDeriv::LIN_VEL_SPLINE: {
            /**
             * when rgbds or radars are involved in the calibration, a linear velocity spline would
             * be maintained in the estimator
             */
            for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
                this->AddRadarFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
            }
            for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
                this->AddAcceFactor<TimeDeriv::LIN_VEL_SPLINE>(estimator, topic, optOption);
                this->AddGyroFactor(estimator, topic, optOption);
            }
            for (const auto &[topic, corrs] : rgbdCorrs) {
                this->AddRGBDOpticalFlowFactor<TimeDeriv::LIN_VEL_SPLINE,
                                               OPTICAL_FLOW_EST_INV_DEPTH>(
                    estimator, topic, corrs, RefineReadoutTimeOptForCameras(topic, optOption));
            }
            for (const auto &[topic, corrs] : visualVelCorrs) {
                this->AddVisualOpticalFlowFactor<TimeDeriv::LIN_VEL_SPLINE,
                                                 OPTICAL_FLOW_EST_INV_DEPTH>(
                    estimator, topic, corrs, RefineReadoutTimeOptForCameras(topic, optOption));
            }
        } break;
        case TimeDeriv::LIN_POS_SPLINE: {
            /*
             * when lidars or optical cameras are involved in the calibration, a translation spline
             * would be maintained in the estimator
             */
            for (const auto &[topic, corrs] : lidarPtsCorrs) {
                this->AddLiDARPointToSurfelFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic,
                                                                             corrs, optOption);
            }
            for (const auto &[topic, corrs] : visualReprojCorrs) {
                this->AddVisualReprojectionFactor<TimeDeriv::LIN_POS_SPLINE>(
                    estimator, topic, corrs, visualGlobalScale.get(),
                    RefineReadoutTimeOptForCameras(topic, optOption));
            }
            for (const auto &[topic, _] : Configor::DataStream::RadarTopics) {
                this->AddRadarFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic, optOption);
            }
            for (const auto &[topic, _] : Configor::DataStream::IMUTopics) {
                this->AddAcceFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic, optOption);
                this->AddGyroFactor(estimator, topic, optOption);
            }
            for (const auto &[topic, corrs] : rgbdCorrs) {
                this->AddRGBDOpticalFlowFactor<TimeDeriv::LIN_POS_SPLINE,
                                               OPTICAL_FLOW_EST_INV_DEPTH>(
                    estimator, topic, corrs, RefineReadoutTimeOptForCameras(topic, optOption));
                /**
                 * when pos spline is maintained, we add additional reprojection constraints for
                 * optical flow tracking correspondence
                 */
                this->AddRGBDOpticalFlowReprojFactor<TimeDeriv::LIN_POS_SPLINE,
                                                     OPTICAL_FLOW_EST_INV_DEPTH>(
                    estimator, topic, corrs, RefineReadoutTimeOptForCameras(topic, optOption));
            }
            for (const auto &[topic, corrs] : visualVelCorrs) {
                this->AddVisualOpticalFlowFactor<TimeDeriv::LIN_POS_SPLINE,
                                                 OPTICAL_FLOW_EST_INV_DEPTH>(
                    estimator, topic, corrs, RefineReadoutTimeOptForCameras(topic, optOption));
                /**
                 * when pos spline is maintained, we add additional reprojection constraints for
                 * optical flow tracking correspondence
                 */
                this->AddVisualOpticalFlowReprojFactor<TimeDeriv::LIN_POS_SPLINE,
                                                       OPTICAL_FLOW_EST_INV_DEPTH>(
                    estimator, topic, corrs, RefineReadoutTimeOptForCameras(topic, optOption));
            }
            if (rgbdPtsCorrs != std::nullopt) {
                /**
                 * such point-to-surfel data association is not necessary, if it exists, we add to
                 * the estimator
                 */
                for (const auto &[topic, corrs] : *rgbdPtsCorrs) {
                    this->AddRGBDPointToSurfelFactor<TimeDeriv::LIN_POS_SPLINE>(estimator, topic,
                                                                                corrs, optOption);
                }
            }
        } break;
    }

    // make this problem full rank
    estimator->SetRefIMUParamsConstant();

    auto sum = estimator->Solve(_ceresOption, this->_priori);
    spdlog::info("here is the summary:\n{}\n", sum.BriefReport());

    // align states to the gravity after the batch optimization is finished
    AlignStatesToGravity();

    // for better map consistency in visualization, we update the veta every time
    for (const auto &[topic, reprojCorrVec] : visualReprojCorrs) {
        auto &veta = _dataMagr->GetSfMData(topic);
        auto &intri = _parMagr->INTRI.Camera.at(topic);
        // compute the pose sequence based on the estimated bsplines and extrinsics
        for (const auto &[viewId, view] : veta->views) {
            auto SE3_CurCmToW = CurCmToW(view->timestamp, topic);
            if (SE3_CurCmToW == std::nullopt) {
                throw Status(Status::CRITICAL,
                             "can not find pose from B-splines for camera '{}'!!!", topic);
            }
            veta->poses.at(view->poseId) =
                ns_veta::Posed(SE3_CurCmToW->so3(), SE3_CurCmToW->translation());
        }

        for (const auto &reprojCorrSeq : reprojCorrVec) {
            // recover point in camera frame
            Eigen::Vector2d pInCamPlane = intri->ImgToCam(reprojCorrSeq->firObv.x);
            double depth = *visualGlobalScale * 1.0 / *reprojCorrSeq->invDepthFir;
            Eigen::Vector3d pInCam(pInCamPlane(0) * depth, pInCamPlane(1) * depth, depth);
            // transform point to world frame (we do not consider the RS effect here, which only
            // affects the visualization)
            auto pose = veta->poses.at(veta->views.at(reprojCorrSeq->firObvViewId)->poseId);
            Eigen::Vector3d pInW = pose.Rotation() * pInCam + pose.Translation();
            veta->structure.at(reprojCorrSeq->lmId).X = pInW;
        }
    }

    // update depth information for rgbds
    for (const auto &[topic, corrs] : rgbdCorrs) {
        for (const auto &corr : corrs) {
            if constexpr (OPTICAL_FLOW_EST_INV_DEPTH) {
                // the inverse depth is estimated, we update the depth
                corr->depth = corr->invDepth > 1E-3 ? 1.0 / corr->invDepth : -1.0;
            } else {
                // the depth is estimated, we update the inverse depth
                corr->invDepth = corr->depth > 1E-3 ? 1.0 / corr->depth : -1.0;
            }
        }
    }

    // update depth information for vel cameras
    for (const auto &[topic, corrs] : visualVelCorrs) {
        for (const auto &corr : corrs) {
            if constexpr (OPTICAL_FLOW_EST_INV_DEPTH) {
                // the inverse depth is estimated, we update the depth
                corr->depth = corr->invDepth > 1E-3 ? 1.0 / corr->invDepth : -1.0;
            } else {
                // the depth is estimated, we update the inverse depth
                corr->invDepth = corr->depth > 1E-3 ? 1.0 / corr->depth : -1.0;
            }
        }
    }

    // these quantities need to be backup for Hessian matrix finding in ceres
    auto backUp = std::make_shared<BackUp>();
    backUp->estimator = estimator;
    backUp->visualGlobalScale = visualGlobalScale;
    backUp->lidarCorrs = lidarPtsCorrs;
    backUp->visualCorrs = visualReprojCorrs;
    backUp->ofCorrs = rgbdCorrs;
    for (const auto &[topic, corrs] : visualVelCorrs) {
        backUp->ofCorrs.insert({topic, corrs});
    }
    // do not back up the maps
    backUp->lidarMap = nullptr;
    backUp->radarMap = nullptr;
    return backUp;
}
}  // namespace ns_ikalibr