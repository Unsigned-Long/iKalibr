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

#include "core/visual_reproj_association.h"
#include "factor/data_correspondence.h"
#include "veta/veta.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

VisualReProjAssociator::VisualReProjAssociator(const CameraModelType &type)
    : ExposureFactor(CameraModel::RSCameraExposureFactor(type)) {}

VisualReProjAssociator::Ptr VisualReProjAssociator::Create(const CameraModelType &type) {
    return std::make_shared<VisualReProjAssociator>(type);
}

std::vector<VisualReProjCorrSeq::Ptr> VisualReProjAssociator::Association(
    const ns_veta::Veta &veta, const ns_veta::PinholeIntrinsic::Ptr &intri) const {
    // scale weight from image pixel to real scale
    const double weight = intri->ImagePlaneToCameraPlaneError(1.0);

    std::vector<VisualReProjCorrSeq::Ptr> corrVec;
    corrVec.reserve(veta.structure.size());

    for (const auto &[lmId, lm] : veta.structure) {
        auto begIter = lm.obs.cbegin();
        const auto &[viewIdFir, featFir] = *begIter;
        const auto &viewFir = veta.views.find(viewIdFir)->second;
        // row / image height - ExposureFactor
        // attention: computed based on raw pixel rather undistorted pixel
        const double lFir =
            intri->GetDistoPixel(featFir.x)(1) / static_cast<double>(viewFir->imgHeight) -
            ExposureFactor;

        auto corrSeq = std::make_shared<VisualReProjCorrSeq>();

        // bring landmark from world frame to the first camera frame which first obverses this
        // landmark
        Eigen::Vector3d lmInFir = veta.poses.at(viewFir->poseId).Inverse().operator()(lm.X);
        // inverse depth
        corrSeq->invDepthFir = std::make_unique<double>(1.0 / lmInFir(2));
        corrSeq->lmId = lmId;
        corrSeq->corrs.reserve(lm.obs.size() - 1);
        corrSeq->firObvViewId = viewIdFir;
        corrSeq->firObv = featFir;

        for (auto curIter = std::next(begIter); curIter != lm.obs.cend(); ++curIter) {
            const auto &[viewIdCur, featCur] = *curIter;
            const auto &viewCur = veta.views.find(viewIdCur)->second;
            // row / image height - ExposureFactor
            // attention: computed based on raw pixel rather undistorted pixel
            const double lCur =
                intri->GetDistoPixel(featCur.x)(1) / static_cast<double>(viewCur->imgHeight) -
                ExposureFactor;

            corrSeq->corrs.push_back(VisualReProjCorr::Create(
                // timestamps
                viewFir->timestamp, viewCur->timestamp,
                // feature location in image plane (has been undistorted)
                featFir.x, featCur.x,
                // row / image height - ExposureFactor: v/h - ExposureFactor
                lFir, lCur,
                // rough weight
                weight));
        }
        corrVec.push_back(corrSeq);
    }

    return corrVec;
}
}  // namespace ns_ikalibr