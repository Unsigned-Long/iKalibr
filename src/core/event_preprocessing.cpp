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

#include "core/event_preprocessing.h"
#include "veta/camera/pinhole.h"
#include "opencv2/imgproc.hpp"
#include "veta/camera/pinhole_brown.h"
#include "veta/camera/pinhole_fisheye.h"
#include "util/status.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui.hpp"
#include "opengv/sac/Ransac.hpp"

#include <solver/batch_opt_option.hpp>

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
/**
 * ActiveEventSurface
 */
ActiveEventSurface::ActiveEventSurface(const ns_veta::PinholeIntrinsic::Ptr &intri,
                                       double filterThd)
    : FILTER_THD(filterThd),
      _intri(intri),
      _undistoMap(VisualUndistortionMap::Create(_intri)),
      _eventImgMat(cv::Size(_intri->imgWidth, _intri->imgHeight), CV_8UC3, cv::Scalar(0, 0, 0)) {
    _sae[0] = Eigen::MatrixXd::Zero(_intri->imgWidth, _intri->imgHeight);
    _sae[1] = Eigen::MatrixXd::Zero(_intri->imgWidth, _intri->imgHeight);
    _saeLatest[0] = Eigen::MatrixXd::Zero(_intri->imgWidth, _intri->imgHeight);
    _saeLatest[1] = Eigen::MatrixXd::Zero(_intri->imgWidth, _intri->imgHeight);
    _timeLatest = 0.0;
}

ActiveEventSurface::Ptr ActiveEventSurface::Create(const ns_veta::PinholeIntrinsicPtr &intri,
                                                   double filterThd) {
    return std::make_shared<ActiveEventSurface>(intri, filterThd);
}

void ActiveEventSurface::GrabEvent(const Event::Ptr &event, bool drawEventMat) {
    const bool ep = event->GetPolarity();
    const std::uint16_t ex = event->GetPos()(0), ey = event->GetPos()(1);
    const double et = event->GetTimestamp();

    // update Surface of Active Events
    const int pol = ep ? 1 : 0;
    const int polInv = !ep ? 1 : 0;
    double &tLast = _saeLatest[pol](ex, ey);
    double &tLastInv = _saeLatest[polInv](ex, ey);

    if ((et > tLast + FILTER_THD) || (tLastInv > tLast)) {
        tLast = et;
        _sae[pol](ex, ey) = et;
    } else {
        tLast = et;
    }
    _timeLatest = et;

    if (drawEventMat) {
        // draw image
        _eventImgMat.at<cv::Vec3b>(cv::Point2d(ex, ey)) =
            ep ? cv::Vec3b(255, 0, 0) : cv::Vec3b(0, 0, 255);
    }
}

void ActiveEventSurface::GrabEvent(const EventArray::Ptr &events, bool drawEventMat) {
    for (const auto &event : events->GetEvents()) {
        GrabEvent(event, drawEventMat);
    }
}

cv::Mat ActiveEventSurface::GetEventImgMat(bool resetMat, bool undistoMat) {
    auto mat = _eventImgMat.clone();
    if (resetMat) {
        _eventImgMat =
            cv::Mat(cv::Size(_intri->imgWidth, _intri->imgHeight), CV_8UC3, cv::Scalar(0, 0, 0));
    }
    if (undistoMat) {
        return _undistoMap->RemoveDistortion(mat);
    } else {
        return mat;
    }
}

cv::Mat ActiveEventSurface::TimeSurface(bool ignorePolarity,
                                        bool undistoMat,
                                        int medianBlurKernelSize,
                                        double decaySec) {
    // create exponential-decayed Time Surface map
    const auto imgSize = cv::Size(_intri->imgWidth, _intri->imgHeight);
    cv::Mat timeSurfaceMap = cv::Mat::zeros(imgSize, CV_64F);

    for (int y = 0; y < imgSize.height; ++y) {
        for (int x = 0; x < imgSize.width; ++x) {
            const double mostRecentStampAtCoord = std::max(_sae[1](x, y), _sae[0](x, y));

            const double dt = _timeLatest - mostRecentStampAtCoord;
            double expVal = std::exp(-dt / decaySec);

            if (!ignorePolarity) {
                double polarity = _sae[1](x, y) > _sae[0](x, y) ? 1.0 : -1.0;
                expVal *= polarity;
            }
            timeSurfaceMap.at<double>(y, x) = expVal;
        }
    }

    if (!ignorePolarity) {
        timeSurfaceMap = 255.0 * (timeSurfaceMap + 1.0) / 2.0;
    } else {
        timeSurfaceMap = 255.0 * timeSurfaceMap;
    }
    timeSurfaceMap.convertTo(timeSurfaceMap, CV_8U);

    if (medianBlurKernelSize > 0) {
        cv::medianBlur(timeSurfaceMap, timeSurfaceMap, 2 * medianBlurKernelSize + 1);
    }

    if (undistoMat) {
        return _undistoMap->RemoveDistortion(timeSurfaceMap);
    } else {
        return timeSurfaceMap;
    }
}

std::pair<cv::Mat, cv::Mat> ActiveEventSurface::RawTimeSurface(bool ignorePolarity,
                                                               bool undistoMat) {
    // create exponential-decayed Time Surface map
    const auto imgSize = cv::Size(_intri->imgWidth, _intri->imgHeight);
    cv::Mat timeSurfaceMap = cv::Mat::zeros(imgSize, CV_64FC1);
    cv::Mat polarityMap = cv::Mat::zeros(imgSize, CV_8UC1);

    for (int y = 0; y < imgSize.height; ++y) {
        for (int x = 0; x < imgSize.width; ++x) {
            double mostRecentStampAtCoord = std::max(_sae[1](x, y), _sae[0](x, y));
            double polarity = _sae[1](x, y) > _sae[0](x, y) ? 1.0 : -1.0;
            polarityMap.at<uchar>(y, x) = polarity;
            if (!ignorePolarity) {
                mostRecentStampAtCoord *= polarity;
            }
            timeSurfaceMap.at<double>(y, x) = mostRecentStampAtCoord;
        }
    }
    if (undistoMat) {
        return {_undistoMap->RemoveDistortion(timeSurfaceMap),
                _undistoMap->RemoveDistortion(polarityMap)};
    } else {
        return {timeSurfaceMap, polarityMap};
    }
}

double ActiveEventSurface::GetTimeLatest() const { return _timeLatest; }

EventNormFlow::NormFlow::NormFlow(double timestamp,
                                  const Eigen::Vector2i &p,
                                  const Eigen::Vector2d &nf)
    : timestamp(timestamp),
      p(p),
      nf(nf) {}

EventNormFlow::NormFlow::Ptr EventNormFlow::NormFlow::Create(double timestamp,
                                                             const Eigen::Vector2i &p,
                                                             const Eigen::Vector2d &nf) {
    return std::make_shared<NormFlow>(timestamp, p, nf);
}

EventArray::Ptr EventNormFlow::NormFlowPack::ActiveEvents(double dt) const {
    std::list<Event::Ptr> events;
    const int rows = rawTimeSurfaceMap.rows;
    const int cols = rawTimeSurfaceMap.cols;
    for (int ey = 0; ey < rows; ey++) {
        for (int ex = 0; ex < cols; ex++) {
            const auto &et = rawTimeSurfaceMap.at<double>(ey, ex);
            // whether this pixel is assigned
            if (et < 1E-3 || timestamp - et > dt) {
                continue;
            }
            const auto &ep = polarityMap.at<uchar>(ey, ex);
            events.push_back(Event::Create(et, {ex, ey}, ep));
        }
    }
    if (events.size() > 0) {
        std::vector<Event::Ptr> eventsVec(events.cbegin(), events.cend());
        return EventArray::Create(eventsVec.back()->GetTimestamp(), eventsVec);
    } else {
        return nullptr;
    }
}

EventArray::Ptr EventNormFlow::NormFlowPack::NormFlowEvents() const {
    std::list<Event::Ptr> events;
    const int rows = rawTimeSurfaceMap.rows;
    const int cols = rawTimeSurfaceMap.cols;
    for (int ey = 0; ey < rows; ey++) {
        for (int ex = 0; ex < cols; ex++) {
            const auto &et = rawTimeSurfaceMap.at<double>(ey, ex);
            // whether this pixel is assigned
            if (et < 1E-3) {
                continue;
            }
            if (inliersOccupy.at<uchar>(ey, ex) == 0) {
                continue;
            }
            const auto &ep = polarityMap.at<uchar>(ey, ex);
            events.push_back(Event::Create(et, {ex, ey}, ep));
        }
    }
    if (events.size() > 0) {
        std::vector<Event::Ptr> eventsVec(events.cbegin(), events.cend());
        return EventArray::Create(eventsVec.back()->GetTimestamp(), eventsVec);
    } else {
        return nullptr;
    }
}

cv::Mat EventNormFlow::NormFlowPack::Visualization(double dt) const {
    cv::Mat m1;
    cv::hconcat(nfSeedsImg, nfsImg, m1);

    cv::Mat actEventMat(nfSeedsImg.size(), CV_8UC3, cv::Scalar(0, 0, 0));
    for (const auto &event : this->ActiveEvents(dt)->GetEvents()) {
        auto ex = event->GetPos()(0), ey = event->GetPos()(1);
        auto ep = event->GetPolarity();
        actEventMat.at<cv::Vec3b>(cv::Point2d(ex, ey)) =
            ep ? cv::Vec3b(255, 0, 0) : cv::Vec3b(0, 0, 255);
    }

    cv::Mat nfEventMat(nfSeedsImg.size(), CV_8UC3, cv::Scalar(0, 0, 0));
    for (const auto &event : this->NormFlowEvents()->GetEvents()) {
        auto ex = event->GetPos()(0), ey = event->GetPos()(1);
        auto ep = event->GetPolarity();
        nfEventMat.at<cv::Vec3b>(cv::Point2d(ex, ey)) =
            ep ? cv::Vec3b(255, 0, 0) : cv::Vec3b(0, 0, 255);
    }

    cv::Mat m2;
    cv::hconcat(actEventMat, nfEventMat, m2);

    cv::Mat m3;
    cv::vconcat(m1, m2, m3);

    return m3;
}

EventNormFlow::NormFlowPack EventNormFlow::ExtractNormFlows(double decaySec,
                                                            int winSize,
                                                            int neighborDist,
                                                            double goodRatioThd,
                                                            double timeDistEventToPlaneThd,
                                                            int ransacMaxIter) const {
    // CV_64FC1
    auto [rtsMat, pMat] = _sea->RawTimeSurface(true, true /*todo: undistorted*/);
    // CV_8UC1
    auto tsImg = _sea->TimeSurface(true, true /*todo: undistorted*/, 0, decaySec);

    cv::Mat mask;
    cv::inRange(tsImg, 50, 230, mask);

    cv::cvtColor(tsImg, tsImg, cv::COLOR_GRAY2BGR);
    auto tsImgNfs = tsImg.clone();

    const int ws = winSize;
    const int subTravSize = std::max(winSize, neighborDist);
    const int winSampleCount = (2 * ws + 1) * (2 * ws + 1);
    const int winSampleCountThd = winSampleCount * goodRatioThd;
    const int rows = mask.rows;
    const int cols = mask.cols;
    cv::Mat occupy = cv::Mat::zeros(rows, cols, CV_8UC1);
    cv::Mat inliersOccupy = cv::Mat::zeros(rows, cols, CV_8UC1);
    std::list<NormFlow::Ptr> nfs;
    for (int y = subTravSize; y < mask.rows - subTravSize; y++) {
        for (int x = subTravSize; x < mask.cols - subTravSize; x++) {
            if (mask.at<uchar>(y /*row*/, x /*col*/) != 255) {
                continue;
            }
            // for this window, obtain the values [x, y, timestamp]
            std::vector<std::tuple<int, int, double>> inRangeData;
            double timeCen = 0.0;
            inRangeData.reserve(winSampleCount);
            bool jumpCurPixel = false;
            for (int dy = -subTravSize; dy <= subTravSize; ++dy) {
                for (int dx = -subTravSize; dx <= subTravSize; ++dx) {
                    int nx = x + dx;
                    int ny = y + dy;

                    // this pixel is in neighbor range
                    if (std::abs(dx) <= neighborDist && std::abs(dy) <= neighborDist) {
                        if (occupy.at<uchar>(ny /*row*/, nx /*col*/) == 255) {
                            // this pixl has been occupied, thus the current pixel would not be
                            // considered in norm flow estimation
                            jumpCurPixel = true;
                            break;
                        }
                    }

                    // this pixel is not considered in the window
                    if (std::abs(dx) > ws || std::abs(dy) > ws) {
                        continue;
                    }

                    // in window but not involved in norm flow estimation
                    if (mask.at<uchar>(ny /*row*/, nx /*col*/) != 255) {
                        continue;
                    }

                    double timestamp = rtsMat.at<double>(ny /*row*/, nx /*col*/);
                    inRangeData.emplace_back(nx, ny, timestamp);
                    if (nx == x && ny == y) {
                        timeCen = timestamp;
                    }
                }
                if (jumpCurPixel) {
                    break;
                }
            }
            // data in this window is sufficient
            if (jumpCurPixel || static_cast<int>(inRangeData.size()) < winSampleCountThd) {
                continue;
            }
            /**
             * drawing
             */
            tsImg.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 255);  // selected but not verified
            occupy.at<uchar>(y /*row*/, x /*col*/) = 255;

            // try fit planes using ransac
            opengv::sac::Ransac<EventLocalPlaneSacProblem> ransac;
            std::shared_ptr<EventLocalPlaneSacProblem> probPtr(
                new EventLocalPlaneSacProblem(inRangeData));
            ransac.sac_model_ = probPtr;
            // the point to plane threshold in temporal domain
            ransac.threshold_ = timeDistEventToPlaneThd;
            ransac.max_iterations_ = ransacMaxIter;
            auto res = ransac.computeModel();

            if (!res || ransac.inliers_.size() / (double)inRangeData.size() < goodRatioThd) {
                continue;
            }
            // success
            Eigen::Vector3d abd;
            probPtr->optimizeModelCoefficients(ransac.inliers_, ransac.model_coefficients_, abd);

            // 'abd' is the params we are interested in
            const double dtdx = -abd(0), dtdy = -abd(1);
            Eigen::Vector2d nf = 1.0 / (dtdx * dtdx + dtdy * dtdy) * Eigen::Vector2d(dtdx, dtdy);
            nfs.push_back(NormFlow::Create(timeCen, Eigen::Vector2i{x, y}, nf));

            for (int idx : ransac.inliers_) {
                const auto &[ex, ey, et] = inRangeData.at(idx);
                inliersOccupy.at<uchar>(ey /*row*/, ex /*col*/) = 255;
            }

            /**
             * drawing
             */
            tsImg.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 255, 0);  // selected and verified
            DrawLineOnCVMat(tsImgNfs, Eigen::Vector2d{x, y} + 0.01 * nf, {x, y});
        }
    }
    NormFlowPack pack;
    pack.nfs = nfs;
    pack.inliersOccupy = inliersOccupy;
    pack.polarityMap = pMat;
    pack.rawTimeSurfaceMap = rtsMat;
    pack.timestamp = _sea->GetTimeLatest();
    // for visualization
    pack.nfsImg = tsImgNfs;
    pack.nfSeedsImg = tsImg;
    return pack;
}

/**
 * EventLocalPlaneSacProblem
 */
int EventLocalPlaneSacProblem::getSampleSize() const { return 3; }

double EventLocalPlaneSacProblem::PointToPlaneDistance(
    double x, double y, double t, double A, double B, double C) {
    // double numerator = std::abs(A * x + B * y + t + C);
    // double denominator = std::sqrt(A * A + B * B + 1);
    // return numerator / denominator;
    double tPred = -(A * x + B * y + C);
    return std::abs(t - tPred);
}

bool EventLocalPlaneSacProblem::computeModelCoefficients(const std::vector<int> &indices,
                                                         model_t &outModel) const {
    Eigen::MatrixXd M(indices.size(), 3);
    Eigen::VectorXd b(indices.size());

    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        const auto &[x, y, t] = _data.at(indices.at(i));
        M(i, 0) = x;
        M(i, 1) = y;
        M(i, 2) = 1;
        b(i) = -t;
    }
    outModel = (M.transpose() * M).ldlt().solve(M.transpose() * b);
    return true;
}

void EventLocalPlaneSacProblem::getSelectedDistancesToModel(const model_t &model,
                                                            const std::vector<int> &indices,
                                                            std::vector<double> &scores) const {
    scores.resize(indices.size());
    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        const auto &[x, y, t] = _data.at(indices.at(i));
        scores.at(i) = PointToPlaneDistance(x, y, t, model(0), model(1), model(2));
    }
}

void EventLocalPlaneSacProblem::optimizeModelCoefficients(const std::vector<int> &inliers,
                                                          const model_t &model,
                                                          model_t &optimized_model) {
    computeModelCoefficients(inliers, optimized_model);
}
}  // namespace ns_ikalibr