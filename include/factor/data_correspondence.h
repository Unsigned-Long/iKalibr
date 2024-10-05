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

#ifndef DATA_CORRESPONDENCE_H
#define DATA_CORRESPONDENCE_H

#include "util/utils_tpl.hpp"
#include "ufo/map/octree/node.h"
#include "veta/landmark.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
class CameraFrame;
using CameraFramePtr = std::shared_ptr<CameraFrame>;

struct PointToSurfelCorr {
    using Ptr = std::shared_ptr<PointToSurfelCorr>;

public:
    double timestamp;

    Eigen::Vector3d pInScan;
    double weight;
    // [norm dir, dist]
    Eigen::Vector4d surfelInW;

    // just for visualization
    Eigen::Vector3d pInMap;
    // for ufomap associator
    ufo::map::Node node;

    PointToSurfelCorr(double timestamp,
                      Eigen::Vector3d pInScan,
                      double weight,
                      Eigen::Vector4d surfelInW);

    static Ptr Create(double timestamp,
                      const Eigen::Vector3d &pInScan,
                      double weight,
                      const Eigen::Vector4d &surfelInW);
};

struct VisualReProjCorr {
public:
    using Ptr = std::shared_ptr<VisualReProjCorr>;

public:
    double ti{}, tj{};

    Eigen::Vector2d fi, fj;

    // row / image height - 'ExposureFactor'
    double li{}, lj{};

    double weight{};

public:
    VisualReProjCorr(double ti,
                     double tj,
                     Eigen::Vector2d fi,
                     Eigen::Vector2d fj,
                     double li,
                     double lj,
                     double weight);

    static Ptr Create(double ti,
                      double tj,
                      const Eigen::Vector2d &fi,
                      const Eigen::Vector2d &fj,
                      double li,
                      double lj,
                      double weight);

    VisualReProjCorr();

    template <class T>
    static void TransformImgToCam(const T *FX_INV,
                                  const T *FY_INV,
                                  const T *CX,
                                  const T *CY,
                                  const Eigen::Vector2<T> &feat,
                                  Eigen::Vector3<T> *P) {
        P->operator()(0) = (feat(0) - *CX) * *FX_INV;
        P->operator()(1) = (feat(1) - *CY) * *FY_INV;
        P->operator()(2) = (T)1.0;
    }

    template <class T>
    static void TransformCamToImg(const T *FX,
                                  const T *FY,
                                  const T *CX,
                                  const T *CY,
                                  const Eigen::Vector3<T> &P,
                                  Eigen::Vector2<T> *feat) {
        feat->operator()(0) = *FX * P(0) + *CX;
        feat->operator()(1) = *FY * P(1) + *CY;
    }
};

struct VisualReProjCorrSeq {
public:
    using Ptr = std::shared_ptr<VisualReProjCorrSeq>;

public:
    std::vector<VisualReProjCorr::Ptr> corrs;

    std::unique_ptr<double> invDepthFir;

    ns_veta::IndexT lmId{};

    ns_veta::IndexT firObvViewId{};
    ns_veta::Observation firObv;

public:
    VisualReProjCorrSeq();
};

struct OpticalFlowCorr {
public:
    using Ptr = std::shared_ptr<OpticalFlowCorr>;
    static constexpr int FIR = 0;
    static constexpr int MID = 1;
    static constexpr int LAST = 2;

public:
    std::array<double, 3> timeAry;
    std::array<double, 3> xTraceAry;
    std::array<double, 3> yTraceAry;
    // row / image height - rsExpFactor
    std::array<double, 3> rdFactorAry;
    double depth;
    double invDepth;
    CameraFramePtr frame;
    // if this dynamic is with depth observability
    bool withDepthObservability;
    double weight = 1.0;

public:
    OpticalFlowCorr(const std::array<double, 3> &timeAry,
                    const std::array<double, 3> &xTraceAry,
                    const std::array<double, 3> &yTraceAry,
                    double depth,
                    const CameraFramePtr &frame,
                    double rsExpFactor);

    static Ptr Create(const std::array<double, 3> &timeAry,
                      const std::array<double, 3> &xDynamicAry,
                      const std::array<double, 3> &yDynamicAry,
                      double depth,
                      const CameraFramePtr &frame,
                      double rsExpFactor);

    [[nodiscard]] Eigen::Vector2d FirPoint() const;

    [[nodiscard]] Eigen::Vector2d MidPoint() const;

    [[nodiscard]] Eigen::Vector2d LastPoint() const;

    template <class Type>
    [[nodiscard]] Type MidPointTime(Type readout) const {
        return timeAry[MID] + rdFactorAry[MID] * readout;
    }

    [[nodiscard]] double MidReadoutFactor() const;

    template <class Type>
    [[nodiscard]] Eigen::Vector2<Type> MidPointVel(Type readout) const {
        if constexpr (std::is_same<Type, double>::value) {
            std::array<Type, 3> newTimeAry{};
            for (int i = 0; i < 3; ++i) {
                newTimeAry[i] = timeAry[i] + rdFactorAry[i] * readout;
            }
            return {ns_ikalibr::LagrangePolynomialTripleMidFOD<Type>(newTimeAry, xTraceAry),
                    ns_ikalibr::LagrangePolynomialTripleMidFOD<Type>(newTimeAry, yTraceAry)};
        } else {
            std::array<Type, 3> newTimeAry{};
            std::array<Type, 3> newXAry{};
            std::array<Type, 3> newYAry{};
            for (int i = 0; i < 3; ++i) {
                newTimeAry[i] = timeAry[i] + rdFactorAry[i] * readout;
                newXAry[i] = (Type)xTraceAry[i];
                newYAry[i] = (Type)yTraceAry[i];
            }
            return {ns_ikalibr::LagrangePolynomialTripleMidFOD<Type>(newTimeAry, newXAry),
                    ns_ikalibr::LagrangePolynomialTripleMidFOD<Type>(newTimeAry, newYAry)};
        }
    }

    template <class T>
    static void SubAMat(
        const T *fx, const T *fy, const T &up, const T &vp, Eigen::Matrix<T, 2, 3> *aMat) {
        *aMat = Eigen::Matrix<T, 2, 3>::Zero();
        (*aMat)(0, 0) = -*fx;
        (*aMat)(1, 1) = -*fy;
        (*aMat)(0, 2) = up;
        (*aMat)(1, 2) = vp;
    }

    template <class T>
    static void SubBMat(
        const T *fx, const T *fy, const T &up, const T &vp, Eigen::Matrix<T, 2, 3> *bMat) {
        *bMat = Eigen::Matrix<T, 2, 3>::Zero();
        (*bMat)(0, 0) = up * vp / *fy;
        (*bMat)(0, 1) = -*fx - up * up / *fx;
        (*bMat)(0, 2) = *fx * vp / *fy;
        (*bMat)(1, 0) = *fy + vp * vp / *fy;
        (*bMat)(1, 1) = -up * vp / *fx;
        (*bMat)(1, 2) = -*fy * up / *fx;
    }

    template <class T>
    static void SubMats(const T *fx,
                        const T *fy,
                        const T *cx,
                        const T *cy,
                        const Eigen::Vector2<T> feat,
                        Eigen::Matrix<T, 2, 3> *aMat,
                        Eigen::Matrix<T, 2, 3> *bMat) {
        const T up = feat(0) - *cx, vp = feat(1) - *cy;
        SubAMat(fx, fy, up, vp, aMat);
        SubBMat(fx, fy, up, vp, bMat);
    }
};

}  // namespace ns_ikalibr

#endif  // DATA_CORRESPONDENCE_H
