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

#include "core/ev_circle_tracking.h"
#include "factor/data_correspondence.h"
#include <opencv2/highgui.hpp>
#include <tiny-viewer/entity/entity.h>

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
/**
 * EventCircleTracking::CircleClusterInfo
 */
EventCircleTracking::CircleClusterInfo::CircleClusterInfo(const std::list<NormFlowPtr>& nf_cluster,
                                                          bool polarity,
                                                          const Eigen::Vector2d& center,
                                                          const Eigen::Vector2d& dir)
    : nfCluster(nf_cluster),
      polarity(polarity),
      center(center),
      dir(dir) {}

EventCircleTracking::CircleClusterInfo::Ptr EventCircleTracking::CircleClusterInfo::Create(
    const std::list<NormFlowPtr>& nf_cluster,
    bool polarity,
    const Eigen::Vector2d& center,
    const Eigen::Vector2d& dir) {
    return std::make_shared<CircleClusterInfo>(nf_cluster, polarity, center, dir);
}

/**
 * EventCircleTracking
 */
void EventCircleTracking::ExtractCircles(const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    auto [pNormFlowCluster, nNormFlowCluster] = this->ClusterNormFlowEvents(nfPack);
    auto pCenDir = ComputeCenterDir(pNormFlowCluster, nfPack);
    auto nCenDir = ComputeCenterDir(nNormFlowCluster, nfPack);

    auto pClusterType = IdentifyCategory(pNormFlowCluster, pCenDir, nfPack);
    auto nClusterType = IdentifyCategory(nNormFlowCluster, nCenDir, nfPack);

    // clean 'NormFlowCluster' and 'CenDir' using 'ClusterType'
    RemoveClusterTypes(pClusterType, pNormFlowCluster, pCenDir, ClusterType::OTHER);
    RemoveClusterTypes(nClusterType, nNormFlowCluster, nCenDir, ClusterType::OTHER);

    std::map<ClusterType, std::vector<CircleClusterInfo::Ptr>> clusters;
    for (int i = 0; i < static_cast<int>(pClusterType.size()); i++) {
        const auto& cenDir = pCenDir[i];
        clusters[pClusterType[i]].push_back(std::make_shared<CircleClusterInfo>(
            pNormFlowCluster[i], true, cenDir.first, cenDir.second));
    }
    for (int i = 0; i < static_cast<int>(nClusterType.size()); i++) {
        const auto& cenDir = nCenDir[i];
        clusters[nClusterType[i]].push_back(std::make_shared<CircleClusterInfo>(
            nNormFlowCluster[i], false, cenDir.first, cenDir.second));
    }

    auto mat = nfPack->tsImg.clone();
    DrawCircleCluster(mat, clusters, nfPack, 10);
    cv::imshow("Circle-Oriented Cluster Results", mat);
}

std::vector<EventCircleTracking::ClusterType> EventCircleTracking::IdentifyCategory(
    const std::vector<std::list<NormFlowPtr>>& clusters,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& cenDirs,
    const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    assert(clusters.size() == cenDirs.size());
    const auto size = static_cast<int>(clusters.size());
    std::vector<ClusterType> types(size, ClusterType::OTHER);
    for (int i = 0; i < size; ++i) {
        types.at(i) = IdentifyCategory(clusters[i], cenDirs[i], nfPack);
    }
    return types;
}

EventCircleTracking::ClusterType EventCircleTracking::IdentifyCategory(
    const std::list<NormFlowPtr>& clusters,
    const std::pair<Eigen::Vector2d, Eigen::Vector2d>& cenDir,
    const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    const Eigen::Vector2d &cen = cenDir.first, dir = cenDir.second;

    enum Side : int { LEFT = 0, RIGHT = 1, ON_LINE = 2 };

    auto CheckPointPosition = [](const Eigen::Vector2d& cen, const Eigen::Vector2d& dir,
                                 const Eigen::Vector2d& p) {
        Eigen::Vector2d A = p - cen;
        Eigen::Vector2d B = dir;

        double crossProduct = A.x() * B.y() - A.y() * B.x();

        if (crossProduct > 0) {
            return LEFT;
        } else if (crossProduct < 0) {
            return RIGHT;
        } else {
            return ON_LINE;
        }
    };
    auto CheckDirectionPosition = [](const Eigen::Vector2d& dir, const Eigen::Vector2d& v) {
        double crossProduct = v.x() * dir.y() - v.y() * dir.x();

        if (crossProduct > 0) {
            return LEFT;
        } else if (crossProduct < 0) {
            return RIGHT;
        } else {
            return ON_LINE;
        }
    };

    // row: (avgDir, nfDir), col: (line, point)
    Eigen::Matrix2i typeCount = Eigen::Matrix2i::Zero();
    for (const auto& nf : clusters) {
        // row: (avgDir, nfDir)
        auto dirPosition = CheckDirectionPosition(dir, nf->nfDir);
        if (dirPosition == ON_LINE) {
            continue;
        }
        for (const auto& [ex, ey, et] : nfPack->nfs.at(nf)) {
            // col: (line, point)
            auto ptPosition = CheckPointPosition(cen, dir, Eigen::Vector2d(ex, ey));
            if (ptPosition == ON_LINE) {
                continue;
            }
            typeCount(dirPosition, ptPosition) += 1;
        }
    }

    Eigen::Matrix2d featMatrix = typeCount.cast<double>().normalized();
    Eigen::EigenSolver<Eigen::Matrix2d> solver(featMatrix);
    Eigen::Vector2d eigenvalues = solver.eigenvalues().real();

    double eigenValAbsDiff = std::abs(std::abs(eigenvalues(0)) - std::abs(eigenvalues(1)));
    if (eigenValAbsDiff < 0.5) {
        return (eigenvalues(0) * eigenvalues(1)) > 0 ? ClusterType::RUN : ClusterType::CHASE;
    } else {
        return ClusterType::OTHER;
    }
}

std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> EventCircleTracking::ComputeCenterDir(
    const std::vector<std::list<NormFlowPtr>>& clusters,
    const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> results;
    results.reserve(clusters.size());
    for (const auto& cluster : clusters) {
        results.emplace_back(ComputeCenterDir(cluster, nfPack));
    }
    return results;
}

std::pair<Eigen::Vector2d, Eigen::Vector2d> EventCircleTracking::ComputeCenterDir(
    const std::list<NormFlowPtr>& cluster, const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    Eigen::Vector2d center(0.0, 0.0), dir(0.0, 0.0);
    int size = 0;
    for (const auto& nf : cluster) {
        auto inliers = nfPack->nfs.at(nf);
        // There is spatial overlap between inliers of norm flows, but it is not a big problem
        for (const auto& [ex, ey, et] : inliers) {
            center(0) += ex, center(1) += ey;
        }
        size += inliers.size();
        dir += inliers.size() * nf->nfDir.normalized();
    }
    center /= size;
    dir /= size;
    dir.normalize();
    return {center, dir};
}

std::pair<std::vector<std::list<NormFlow::Ptr>>, std::vector<std::list<NormFlow::Ptr>>>
EventCircleTracking::ClusterNormFlowEvents(const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    cv::Mat pMat(nfPack->rawTimeSurfaceMap.size(), CV_8UC1, cv::Scalar(0));
    cv::Mat nMat(nfPack->rawTimeSurfaceMap.size(), CV_8UC1, cv::Scalar(0));
    for (const auto& event : nfPack->NormFlowInlierEvents()) {
        Event::PosType::Scalar ex = event->GetPos()(0), ey = event->GetPos()(1);
        if (event->GetPolarity()) {
            pMat.at<uchar>(ey, ex) = 255;
        } else {
            nMat.at<uchar>(ey, ex) = 255;
        }
    }

    auto pContours = FindContours(pMat), nContours = FindContours(nMat);
    FilterContoursUsingArea(pContours, CLUSTER_AREA_THD);
    FilterContoursUsingArea(nContours, CLUSTER_AREA_THD);

    std::vector<std::list<NormFlow::Ptr>> pNormFlowCluster, ncNormFlowCluster;
    pNormFlowCluster.resize(pContours.size()), ncNormFlowCluster.resize(nContours.size());
    for (const auto& [nf, inliers] : nfPack->nfs) {
        const auto &x = nf->p(0), y = nf->p(1);
        std::vector<std::list<NormFlow::Ptr>>* nfCluster;
        std::vector<std::vector<cv::Point>>* contours;
        if (nfPack->polarityMap.at<uchar>(y, x)) {
            nfCluster = &pNormFlowCluster;
            contours = &pContours;
        } else {
            nfCluster = &ncNormFlowCluster;
            contours = &nContours;
        }
        for (int i = 0; i != static_cast<int>(contours->size()); ++i) {
            if (cv::pointPolygonTest(contours->at(i), cv::Point(x, y), false) < 0) {
                // out of contours
                continue;
            }
            nfCluster->at(i).push_back(nf);
        }
    }

    // auto mat = nfPack->NormFlowInlierEventMat();
    // DrawContours(mat, pContours, {255, 255, 255}), DrawContours(mat, nContours, {0, 255, 0});
    // cv::imshow("Cluster Contours", mat);
    // cv::imshow("pMat", pMat);
    // cv::imshow("nMat", nMat);

    return {pNormFlowCluster, ncNormFlowCluster};
}

void EventCircleTracking::FilterContoursUsingArea(std::vector<std::vector<cv::Point>>& contours,
                                                  int areaThd) {
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [areaThd](const std::vector<cv::Point>& contour) {
                                      return cv::contourArea(contour) < areaThd;
                                  }),
                   contours.end());
}

std::vector<std::vector<cv::Point>> EventCircleTracking::FindContours(const cv::Mat& binaryImg) {
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(binaryImg, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE);
    return contours;
}

void EventCircleTracking::DrawCircleCluster(
    cv::Mat& mat,
    const std::map<ClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
    const EventNormFlow::NormFlowPack::Ptr& nfPack,
    double scale) {
    for (const auto& [type, cluster] : clusters) {
        cv::Scalar color;
        switch (type) {
            case ClusterType::CHASE:
                color = cv::Scalar(127, 127, 127);
                break;
            case ClusterType::RUN:
                color = cv::Scalar(255, 255, 255);
                break;
            case ClusterType::OTHER:
                color = cv::Scalar(0, 0, 0);
                break;
        }

        for (const auto& c : cluster) {
            DrawCluster(mat, c->nfCluster, nfPack);

            DrawLineOnCVMat(mat, c->center, c->center + scale * c->dir, color);

            const double f = c->polarity ? 1 : -1;
            Eigen::Vector2d p = c->center + f * scale * Eigen::Vector2d(-c->dir(1), c->dir(0));
            DrawLineOnCVMat(mat, c->center, p, color);

            DrawKeypointOnCVMat(mat, c->center, false, color);
        }
    }
}

void EventCircleTracking::DrawCenterDir(
    cv::Mat& mat,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& cenDirVec,
    double scale) {
    for (const auto& [c, d] : cenDirVec) {
        DrawKeypointOnCVMat(mat, c, false, cv::Scalar(255, 255, 255));
        DrawLineOnCVMat(mat, c, c + scale * d, cv::Scalar(255, 255, 255));
    }
}

void EventCircleTracking::DrawCenterDir(
    cv::Mat& mat,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& cenDirVec,
    const std::vector<ClusterType>& types,
    double scale) {
    assert(cenDirVec.size() == types.size());
    const int size = static_cast<int>(cenDirVec.size());
    for (int i = 0; i < size; ++i) {
        const auto& [c, d] = cenDirVec[i];
        cv::Scalar color;
        switch (types[i]) {
            case ClusterType::CHASE:
                color = cv::Scalar(127, 127, 127);
                break;
            case ClusterType::RUN:
                color = cv::Scalar(255, 255, 255);
                break;
            case ClusterType::OTHER:
                color = cv::Scalar(0, 0, 0);
                break;
        }
        DrawLineOnCVMat(mat, c, c + scale * d, color);
        DrawKeypointOnCVMat(mat, c, false, color);
    }
}

void EventCircleTracking::DrawCluster(cv::Mat& mat,
                                      const std::vector<std::list<NormFlow::Ptr>>& clusters,
                                      const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    for (const auto& cluster : clusters) {
        DrawCluster(mat, cluster, nfPack);
    }
}

void EventCircleTracking::DrawCluster(cv::Mat& mat,
                                      const std::list<NormFlowPtr>& clusters,
                                      const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    auto color = ns_viewer::Entity::GetUniqueColour();
    cv::Vec3b c(color.b * 255, color.g * 255, color.r * 255);
    for (const auto& nf : clusters) {
        for (const auto& [ex, ey, et] : nfPack->nfs.at(nf)) {
            mat.at<cv::Vec3b>(ey, ex) = c;
        }
    }
}

void EventCircleTracking::DrawContours(cv::Mat& mat,
                                       const std::vector<std::vector<cv::Point>>& contours,
                                       const cv::Vec3b& color) {
    for (const auto& cVec : contours) {
        for (const auto& p : cVec) {
            mat.at<cv::Vec3b>(p) = color;
        }
    }
}

}  // namespace ns_ikalibr
