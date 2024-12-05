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
    // todo: Two adjacent circular triggered event clusters may stick together, leading to
    // clustering errors!!!
    auto [pNormFlowCluster, nNormFlowCluster] =
        this->ClusterNormFlowEvents(nfPack, CLUSTER_AREA_THD);
    auto pCenDir = ComputeCenterDir(pNormFlowCluster, nfPack);
    auto nCenDir = ComputeCenterDir(nNormFlowCluster, nfPack);

    auto mat0 = nfPack->tsImg.clone();
    DrawCluster(mat0, pNormFlowCluster, nfPack);
    DrawCluster(mat0, nNormFlowCluster, nfPack);

    auto pClusterType = IdentifyCategory(pNormFlowCluster, pCenDir, nfPack);
    auto nClusterType = IdentifyCategory(nNormFlowCluster, nCenDir, nfPack);

    std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>> clusters;
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

    auto mat1 = nfPack->tsImg.clone();
    DrawCircleCluster(mat1, clusters, nfPack, 10);
    // cv::imshow("Circle Cluster Results", mat1);

    /**
     * First, we search for potential matching pairs within the clusters that are already registered
     * as 'RUN' or 'CHASE' types.
     */
    auto pairs = SearchMatchesInRunChasePair(clusters, std::cos(30.0 /*degree*/ * DEG2RAD));
    RemovingAmbiguousMatches(pairs);

    /**
     * Next, we re-pair those clusters that are already registered as 'RUN' or 'CHASE' types but
     * were not paired in the previous step, attempting to match them with the unregistered
     * clusters.
     */
    std::set<CircleClusterInfo::Ptr> alreadyMatched;
    for (const auto& [c1, c2] : pairs) {
        alreadyMatched.insert(c1), alreadyMatched.insert(c2);
    }
    auto newPairs = ReSearchMatchesCirclesOtherPair(clusters, alreadyMatched,
                                                    std::cos(30.0 /*degree*/ * DEG2RAD));
    RemovingAmbiguousMatches(newPairs);
    pairs.insert(newPairs.begin(), newPairs.end());

    auto mat2 = nfPack->tsImg.clone();
    DrawCircleClusterPair(mat2, pairs, nfPack);
    cv::Mat mat3;
    cv::hconcat(mat0, mat1, mat3);
    cv::hconcat(mat3, mat2, mat3);
    cv::imshow("Initial Cluster | Identification | Matching Results", mat3);
}

std::map<EventCircleTracking::CircleClusterInfo::Ptr, EventCircleTracking::CircleClusterInfo::Ptr>
EventCircleTracking::SearchMatchesInRunChasePair(
    const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
    double DIR_DIFF_COS_THD) {
    if (clusters.find(CircleClusterType::CHASE) == clusters.end() ||
        clusters.find(CircleClusterType::RUN) == clusters.end()) {
        return {};
    }
    // chase, run
    std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> pairs;
    for (const auto& cCluster : clusters.at(CircleClusterType::CHASE)) {
        CircleClusterInfo::Ptr bestRunCluster = nullptr;
        double bestDist = std::numeric_limits<double>::max();

        for (const auto& rCluster : clusters.at(CircleClusterType::RUN)) {
            double dist = TryMatchRunChaseCircleClusterPair(rCluster, cCluster, DIR_DIFF_COS_THD);
            if (dist < 0.0) {
                continue;
            }
            if (dist < bestDist && !ClusterExistsInCurCircle(clusters, rCluster, cCluster)) {
                bestRunCluster = rCluster;
                bestDist = dist;
            }
        }

        if (bestRunCluster != nullptr) {
            pairs[cCluster] = bestRunCluster;
        }
    }
    return pairs;
}

std::map<EventCircleTracking::CircleClusterInfo::Ptr, EventCircleTracking::CircleClusterInfo::Ptr>
EventCircleTracking::ReSearchMatchesCirclesOtherPair(
    const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
    const std::set<CircleClusterInfo::Ptr>& alreadyMatched,
    double DIR_DIFF_COS_THD) {
    if (clusters.find(CircleClusterType::OTHER) == clusters.end() ||
        clusters.find(CircleClusterType::CHASE) == clusters.end() ||
        clusters.find(CircleClusterType::RUN) == clusters.end()) {
        return {};
    }
    // chase, run
    std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> pairs;
    for (const auto& [type, clusterVec] : clusters) {
        if (type == CircleClusterType::OTHER) {
            continue;
        }
        for (const auto& cluster : clusterVec) {
            if (alreadyMatched.find(cluster) != alreadyMatched.end()) {
                // already been matched
                continue;
            }
            /**
             * For those clusters that have been assigned a category but haven't been matched yet,
             * we search for possible matching candidates in the unassigned clusters.
             */
            CircleClusterInfo::Ptr bestCluster = nullptr;
            double bestDist = std::numeric_limits<double>::max();

            for (const auto& oCluster : clusters.at(CircleClusterType::OTHER)) {
                double dist = -1.0;
                if (type == CircleClusterType::CHASE) {
                    // leftCluster as 'CircleClusterType::CHASE' one
                    dist = TryMatchRunChaseCircleClusterPair(oCluster, cluster, DIR_DIFF_COS_THD);
                } else if (type == CircleClusterType::RUN) {
                    // leftCluster as 'CircleClusterType::RUN' one
                    dist = TryMatchRunChaseCircleClusterPair(cluster, oCluster, DIR_DIFF_COS_THD);
                }
                if (dist < 0.0) {
                    continue;
                }
                if (dist < bestDist && !ClusterExistsInCurCircle(clusters, oCluster, cluster)) {
                    bestCluster = oCluster;
                    bestDist = dist;
                }
            }
            if (bestCluster != nullptr) {
                if (type == CircleClusterType::CHASE) {
                    pairs[cluster] = bestCluster;
                } else if (type == CircleClusterType::RUN) {
                    pairs[bestCluster] = cluster;
                }
            }
        }
    }
    return pairs;
}

double EventCircleTracking::TryMatchRunChaseCircleClusterPair(
    const CircleClusterInfo::Ptr& rCluster,
    const CircleClusterInfo::Ptr& cCluster,
    double DIR_DIFF_COS_THD) {
    // different polarity
    if (cCluster->polarity == rCluster->polarity) {
        return -1.0;
    }
    // small direction difference
    const double dirDiffCos = cCluster->dir.dot(rCluster->dir);
    if (dirDiffCos < DIR_DIFF_COS_THD) {
        return -1.0;
    }
    Eigen::Vector2d c2r = (rCluster->center - cCluster->center).normalized();
    if (c2r.dot(cCluster->dir) < DIR_DIFF_COS_THD || c2r.dot(rCluster->dir) < DIR_DIFF_COS_THD) {
        return -1.0;
    }
    double dist = (cCluster->center - rCluster->center).norm();
    return dist;
}

bool EventCircleTracking::ClusterExistsInCurCircle(
    const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
    const CircleClusterInfo::Ptr& p1,
    const CircleClusterInfo::Ptr& p2) {
    Eigen::Vector2d c = (p1->center + p2->center) * 0.5;
    double r2 = (0.5 * (p1->center - p2->center)).squaredNorm();

    for (const auto& [type, cluster] : clusters) {
        for (const auto& p : cluster) {
            if (p == p2 || p == p1) {
                continue;
            }
            if ((p->center - c).squaredNorm() < r2) {
                return true;
            }
        }
    }
    return false;
}

void EventCircleTracking::RemovingAmbiguousMatches(
    std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr>& pairs) {
    std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> reversedPairs;

    for (auto it = pairs.begin(); it != pairs.end(); ++it) {
        CircleClusterInfo::Ptr type1 = it->first;
        CircleClusterInfo::Ptr type2 = it->second;

        if (reversedPairs.find(type2) == reversedPairs.end()) {
            reversedPairs[type2] = type1;
        } else {
            CircleClusterInfo::Ptr prevType1 = reversedPairs[type2];
            if ((type1->center - type2->center).norm() <
                (prevType1->center - type2->center).norm()) {
                reversedPairs[type2] = type1;
            }
        }
    }

    pairs.clear();
    for (auto& entry : reversedPairs) {
        pairs[entry.second] = entry.first;
    }
}

std::vector<EventCircleTracking::CircleClusterType> EventCircleTracking::IdentifyCategory(
    const std::vector<std::list<NormFlowPtr>>& clusters,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& cenDirs,
    const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    assert(clusters.size() == cenDirs.size());
    const auto size = static_cast<int>(clusters.size());
    std::vector<CircleClusterType> types(size, CircleClusterType::OTHER);
    for (int i = 0; i < size; ++i) {
        types.at(i) = IdentifyCategory(clusters[i], cenDirs[i], nfPack);
    }
    return types;
}

EventCircleTracking::CircleClusterType EventCircleTracking::IdentifyCategory(
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
    Eigen::Matrix2d weightTypeCount = Eigen::Matrix2d::Zero();
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
            // const double weight = std::exp(-nf->nfNorm * std::abs(nfPack->timestamp - et));
            const double weight = 1.0;
            weightTypeCount(dirPosition, ptPosition) += weight;
        }
    }

    /**
     * CircleClusterType::RUN
     * [ 0.7, 0.;
     *   0.,  0.7 ]
     *
     * CircleClusterType::CHASE
     * [ 0.,  0.7;
     *   0.7, 0.  ]
     */

    // todo: How to better identify categories
    Eigen::Matrix2d featMatrix = weightTypeCount.normalized();

    // const double t1 = featMatrix(0, 0) + featMatrix(1, 1);
    // const double t2 = featMatrix(0, 1) + featMatrix(1, 0);
    // if (t1 > 2.0 * t2) {
    //     return CircleClusterType::RUN;
    // } else if (t2 > 2.0 * t1) {
    //     return CircleClusterType::CHASE;
    // } else {
    //     return CircleClusterType::OTHER;
    // }

    Eigen::EigenSolver<Eigen::Matrix2d> solver(featMatrix);
    Eigen::Vector2d eigenvalues = solver.eigenvalues().real();
    double eigenValAbsDiff = std::abs(std::abs(eigenvalues(0)) - std::abs(eigenvalues(1)));
    if (eigenValAbsDiff < 0.5) {
        return (eigenvalues(0) * eigenvalues(1)) > 0 ? CircleClusterType::RUN
                                                     : CircleClusterType::CHASE;
    } else {
        return CircleClusterType::OTHER;
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
EventCircleTracking::ClusterNormFlowEvents(const EventNormFlow::NormFlowPack::Ptr& nfPack,
                                           double clusterAreaThd) {
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
    InterruptionInTimeDomain(pMat, nfPack->rawTimeSurfaceMap, 7E-3);
    InterruptionInTimeDomain(nMat, nfPack->rawTimeSurfaceMap, 7E-3);

    auto pContours = FindContours(pMat), nContours = FindContours(nMat);
    FilterContoursUsingArea(pContours, clusterAreaThd);
    FilterContoursUsingArea(nContours, clusterAreaThd);

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

void EventCircleTracking::InterruptionInTimeDomain(cv::Mat& pMat, const cv::Mat& tMat, double thd) {
    assert(pMat.type() == CV_8UC1);
    assert(tMat.type() == CV_64FC1);
    assert(pMat.size() == tMat.size());

    int rows = pMat.rows;
    int cols = pMat.cols;
    auto LargerThanThd = [&pMat, &tMat, &thd](const double ct, int r, int c) {
        if (pMat.at<uchar>(r, c) != 255) {
            return false;
        }
        const double t = tMat.at<double>(r, c);
        if (std::abs(ct - t) < thd) {
            return false;
        }
        return true;
    };
    cv::Mat pMatNew = pMat.clone();
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            if (pMat.at<uchar>(i, j) != 255) {
                continue;
            }
            const double ct = tMat.at<double>(i, j);

            if (LargerThanThd(ct, i - 1, j) || LargerThanThd(ct, i + 1, j) ||  // top and bottom
                LargerThanThd(ct, i, j - 1) || LargerThanThd(ct, i, j + 1) ||  // left and right
                LargerThanThd(ct, i - 1, j - 1) || LargerThanThd(ct, i + 1, j + 1) ||
                LargerThanThd(ct, i - 1, j + 1) || LargerThanThd(ct, i + 1, j - 1)) {
                pMatNew.at<uchar>(i, j) = 0;
            }
        }
    }
    pMat = pMatNew;
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
    const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
    const EventNormFlow::NormFlowPack::Ptr& nfPack,
    double scale) {
    for (const auto& [type, cluster] : clusters) {
        DrawCircleCluster(mat, cluster, type, nfPack, scale);
    }
}

void EventCircleTracking::DrawCircleClusterPair(
    cv::Mat& mat,
    const std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr>& pair,
    const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    for (const auto& [chase, run] : pair) {
        DrawCluster(mat, chase->nfCluster, nfPack);
        DrawCluster(mat, run->nfCluster, nfPack);
    }
    for (const auto& [chase, run] : pair) {
        DrawLineOnCVMat(mat, chase->center, run->center, cv::Scalar(255, 255, 255));

        DrawKeypointOnCVMat(mat, run->center, false, cv::Scalar(255, 255, 255));
        DrawKeypointOnCVMat(mat, chase->center, false, cv::Scalar(127, 127, 127));
    }
}

void EventCircleTracking::DrawCircleCluster(cv::Mat& mat,
                                            const std::vector<CircleClusterInfo::Ptr>& clusters,
                                            CircleClusterType type,
                                            const EventNormFlow::NormFlowPack::Ptr& nfPack,
                                            double scale) {
    cv::Scalar color;
    switch (type) {
        case CircleClusterType::CHASE:
            color = cv::Scalar(127, 127, 127);
            break;
        case CircleClusterType::RUN:
            color = cv::Scalar(255, 255, 255);
            break;
        case CircleClusterType::OTHER:
            color = cv::Scalar(0, 0, 0);
            break;
    }

    for (const auto& c : clusters) {
        if (type == CircleClusterType::OTHER) {
            DrawCluster(mat, c->nfCluster, nfPack, ns_viewer::Colour(1.0f, 1.0f, 1.0f));
        } else {
            DrawCluster(mat, c->nfCluster, nfPack);
        }

        DrawLineOnCVMat(mat, c->center, c->center + scale * c->dir, color);

        const double f = c->polarity ? 1 : -1;
        Eigen::Vector2d p = c->center + f * scale * Eigen::Vector2d(-c->dir(1), c->dir(0));
        DrawLineOnCVMat(mat, c->center, p, color);

        DrawKeypointOnCVMat(mat, c->center, false, color);
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
                                      const EventNormFlow::NormFlowPack::Ptr& nfPack,
                                      std::optional<ns_viewer::Colour> color) {
    if (color == std::nullopt) {
        color = ns_viewer::Entity::GetUniqueColour();
    }
    cv::Vec3b c(color->b * 255, color->g * 255, color->r * 255);
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
