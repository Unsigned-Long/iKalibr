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
#include "viewer/viewer.h"
#include "factor/time_varying_circle_fitting.hpp"
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <tiny-viewer/entity/line.h>

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

void EventCircleTracking::Process(const EventNormFlow::NormFlowPack::Ptr& nfPack,
                                  const Viewer::Ptr& viewer) const {
    auto evsInEachCircleClusterPair = this->ExtractPotentialCircleClusters(
        nfPack, this->CLUSTER_AREA_THD, this->DIR_DIFF_DEG_THD);

    if (viewer != nullptr) {
        // cluster results
        for (const auto& [evs1, evs2] : evsInEachCircleClusterPair) {
            viewer->AddEventData(evs1, nfPack->timestamp, Viewer::VIEW_MAP, {0.01, 100}, {}, 2.0f);
            viewer->AddEventData(evs2, nfPack->timestamp, Viewer::VIEW_MAP, {0.01, 100}, {}, 2.0f);
        }
        viewer->AddEventData(nfPack->ActiveEvents(-1.0), nfPack->timestamp, Viewer::VIEW_MAP,
                             {0.01, 100}, ns_viewer::Colour::Black(), 1.0f);
    }

    std::vector<TimeVaryingCircle::Ptr> tvCircles;
    tvCircles.reserve(evsInEachCircleClusterPair.size());
    for (const auto& [evs1, evs2] : evsInEachCircleClusterPair) {
        if (auto c = FitTimeVaryingCircle(evs1, evs2); c != nullptr) {
            tvCircles.push_back(c);
        }
    }

    // filter overlapped circles
    std::vector<TimeVaryingCircle::Circle> circles;
    circles.reserve(tvCircles.size());
    for (const auto& c : tvCircles) {
        circles.push_back(c->CircleAt(nfPack->timestamp));
    }

    std::list<std::size_t> indices;
    for (int i = 0; i < static_cast<int>(tvCircles.size()); ++i) {
        const auto& curCircle = circles.at(i);
        for (int j = i + 1; j < static_cast<int>(tvCircles.size()); ++j) {
            const auto& tarCircle = circles.at(j);
            const double distCenter = (curCircle.first - tarCircle.first).norm();
            if (distCenter < curCircle.second + tarCircle.second) {
                // overlapped
                indices.push_back(i);
                break;
            }
        }
    }
    std::size_t backValidIdx = circles.size() - 1;
    for (std::size_t index : indices) {
        circles.at(index) = circles.at(backValidIdx);
        tvCircles.at(index) = tvCircles.at(backValidIdx);
        --backValidIdx;
    }
    circles.resize(backValidIdx + 1);
    tvCircles.resize(backValidIdx + 1);

    if (viewer != nullptr) {
        // time-varying ellipse fitting
        for (const auto& ellipse : tvCircles) {
            viewer->AddSpatioTemporalTrace(ellipse->PosVecAt(1E-3), nfPack->timestamp,
                                           Viewer::VIEW_MAP, 2.0f, ns_viewer::Colour::Green(),
                                           {0.01, 100});
        }
    }

    cv::Mat cMat = nfPack->tsImg.clone();
    for (const auto& c : circles) {
        cv::circle(cMat, cv::Point2d(c.first(0), c.first(1)), c.second, cv::Scalar(0, 255, 0), 1);
        DrawKeypointOnCVMat(cMat, c.first, false, cv::Scalar(0, 255, 0));
    }
    cv::imshow("Extracted Circles", cMat);
}

std::vector<std::pair<EventArray::Ptr, EventArray::Ptr>>
EventCircleTracking::ExtractPotentialCircleClusters(const EventNormFlow::NormFlowPack::Ptr& nfPack,
                                                    double CLUSTER_AREA_THD,
                                                    double DIR_DIFF_DEG_THD) {
    auto [pNormFlowCluster, nNormFlowCluster] = ClusterNormFlowEvents(nfPack, CLUSTER_AREA_THD);
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
    const double DIR_DIFF_COS_THD = std::cos(DIR_DIFF_DEG_THD /*degree*/ * DEG2RAD);
    auto pairs = SearchMatchesInRunChasePair(clusters, DIR_DIFF_COS_THD);
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
    auto newPairs = ReSearchMatchesCirclesOtherPair(clusters, alreadyMatched, DIR_DIFF_COS_THD);
    RemovingAmbiguousMatches(newPairs);
    pairs.insert(newPairs.begin(), newPairs.end());

    /**
     * Finally, we attempt to build new potential matching relationships within clusters that have
     * neither been registered nor matched.
     */
    for (const auto& [c1, c2] : newPairs) {
        alreadyMatched.insert(c1), alreadyMatched.insert(c2);
    }
    auto newPairs2 = ReSearchMatchesOtherOtherPair(clusters, alreadyMatched, DIR_DIFF_COS_THD);
    RemovingAmbiguousMatches(newPairs2);
    pairs.insert(newPairs2.begin(), newPairs2.end());

    auto mat2 = nfPack->tsImg.clone();
    DrawCircleClusterPair(mat2, pairs, nfPack);
    cv::Mat mat3;
    cv::hconcat(mat0, mat1, mat3);
    cv::hconcat(mat3, mat2, mat3);
    cv::imshow(
        "ExtractPotentialCircleClusters: Initial Cluster | Identification | Matching Results",
        mat3);

    /**
     * Next, for all the potential circular clusters that have been matched, we retrieve their
     * corresponding original events.
     */
    return RawEventsOfCircleClusterPairs(pairs, nfPack);
}

EventCircleTracking::TimeVaryingCircle::Ptr EventCircleTracking::FitTimeVaryingCircle(
    const EventArray::Ptr& ary1, const EventArray::Ptr& ary2) {
    double st = std::numeric_limits<double>::max();
    double et = std::numeric_limits<double>::min();

    auto ComputeCenter = [&st, &et](const EventArray::Ptr& ary) {
        Eigen::Vector2d c(0.0, 0.0);
        for (const auto& event : ary->GetEvents()) {
            c += event->GetPos().cast<double>();
            if (st > event->GetTimestamp()) {
                st = event->GetTimestamp();
            }
            if (et < event->GetTimestamp()) {
                et = event->GetTimestamp();
            }
        }
        return c / static_cast<double>(ary->GetEvents().size());
    };
    Eigen::Vector2d c1 = ComputeCenter(ary1), c2 = ComputeCenter(ary2);

    const Eigen::Vector2d c = 0.5 * (c1 + c2);
    const double r = (0.5 * (c1 - c2)).norm();

    auto circle = TimeVaryingCircle::Create(st, et, {0.0, c(0)}, {0.0, c(1)}, {0.0, 0.0, r * r});

    ceres::Problem problem;

    auto AddResidualsToProblem = [&circle, &problem](const EventArray::Ptr& ary) {
        for (const auto& event : ary->GetEvents()) {
            auto cf = TimeVaryingCircleFittingFactor::Create(event, 1.0);
            cf->AddParameterBlock(2);
            cf->AddParameterBlock(2);
            cf->AddParameterBlock(3);
            cf->SetNumResiduals(1);

            std::vector<double*> params;
            params.push_back(circle->cx.data());
            params.push_back(circle->cy.data());
            params.push_back(circle->r2.data());

            problem.AddResidualBlock(cf, nullptr, params);
        }
    };

    AddResidualsToProblem(ary1);
    AddResidualsToProblem(ary2);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    // options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    if (circle->RadiusAt(circle->et) < 1) {
        return nullptr;
    } else {
        return circle;
    }
}

std::vector<std::pair<EventArray::Ptr, EventArray::Ptr>>
EventCircleTracking::RawEventsOfCircleClusterPairs(
    const std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr>& pairs,
    const EventNormFlow::NormFlowPack::Ptr& nfPack) {
    cv::Mat occupyMat(nfPack->Rows(), nfPack->Cols(), CV_8UC1, cv::Scalar(0));
    auto RawEventsOfEachCircleClusterPairs = [&occupyMat,
                                              &nfPack](const CircleClusterInfo::Ptr& ccs) {
        std::list<Event::Ptr> clusters;
        for (const auto& nf : ccs->nfCluster) {
            for (const auto& [ex, ey, et] : nfPack->nfs.at(nf)) {
                if (auto& o = occupyMat.at<uchar>(ey, ex); o == 0) {
                    clusters.push_back(Event::Create(et, {ex, ey}, ccs->polarity));
                    o = 255;
                }
            }
        }
        clusters.sort([](const Event::Ptr& e1, const Event::Ptr& e2) {
            return e1->GetTimestamp() < e2->GetTimestamp();
        });
        return clusters;
    };
    auto RemoveOldEvents = [](std::list<Event::Ptr>& clusters, double timestamp) {
        auto it = clusters.begin();
        while (it != clusters.end() && (*it)->GetTimestamp() < timestamp) {
            ++it;
        }
        clusters.erase(clusters.begin(), it);
    };
    std::vector<std::pair<EventArray::Ptr, EventArray::Ptr>> eventsOfCluster;
    eventsOfCluster.reserve(pairs.size());
    for (const auto& [cCluster, rCluster] : pairs) {
        auto cEventAry = RawEventsOfEachCircleClusterPairs(cCluster);
        auto rEventAry = RawEventsOfEachCircleClusterPairs(rCluster);

        double st = std::max(cEventAry.front()->GetTimestamp(), rEventAry.front()->GetTimestamp());

        RemoveOldEvents(cEventAry, st);
        RemoveOldEvents(rEventAry, st);
        // todo: fit circle? rather than ellipse???

        if (cEventAry.empty() || rEventAry.empty()) {
            continue;
        }

        auto cAry = EventArray::Create(cEventAry.back()->GetTimestamp(),
                                       {cEventAry.cbegin(), cEventAry.cend()});

        auto rAry = EventArray::Create(rEventAry.back()->GetTimestamp(),
                                       {rEventAry.cbegin(), rEventAry.cend()});

        eventsOfCluster.push_back({cAry, rAry});
    }
    return eventsOfCluster;
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
                    // if 'bestCluster' already exists
                    auto iter = pairs.find(bestCluster);
                    if (iter == pairs.end()) {
                        pairs[bestCluster] = cluster;
                    } else {
                        const double oldDist = (iter->first->center - iter->second->center).norm();
                        if (bestDist < oldDist) {
                            iter->second = cluster;
                        }
                    }
                }
            }
        }
    }
    return pairs;
}

std::map<EventCircleTracking::CircleClusterInfo::Ptr, EventCircleTracking::CircleClusterInfo::Ptr>
EventCircleTracking::ReSearchMatchesOtherOtherPair(
    const std::map<CircleClusterType, std::vector<CircleClusterInfo::Ptr>>& clusters,
    const std::set<CircleClusterInfo::Ptr>& alreadyMatched,
    double DIR_DIFF_COS_THD) {
    if (clusters.find(CircleClusterType::OTHER) == clusters.end()) {
        return {};
    }
    // chase, run
    std::map<CircleClusterInfo::Ptr, CircleClusterInfo::Ptr> pairs;
    for (const auto& c1 : clusters.at(CircleClusterType::OTHER)) {
        if (alreadyMatched.find(c1) != alreadyMatched.end()) {
            // already been matched
            continue;
        }
        CircleClusterInfo::Ptr bestCluster = nullptr;
        double bestDist = std::numeric_limits<double>::max();
        CircleClusterType bestClusterType = CircleClusterType::OTHER;

        for (const auto& c2 : clusters.at(CircleClusterType::OTHER)) {
            // same
            if (c1 == c2) {
                continue;
            }
            if (alreadyMatched.find(c2) != alreadyMatched.end()) {
                // already been matched
                continue;
            }
            double dist = -1.0;
            const double d1 = TryMatchRunChaseCircleClusterPair(c1, c2, DIR_DIFF_COS_THD);
            const double d2 = TryMatchRunChaseCircleClusterPair(c2, c1, DIR_DIFF_COS_THD);
            if (d1 < 0.0 && d2 < 0.0) {
                // do nothing
            } else if (d1 < 0.0 && d2 > 0.0) {
                // c2: run, c1: chase
                dist = d2;
                bestClusterType = CircleClusterType::RUN;
            } else if (d2 < 0.0 && d1 > 0.0) {
                // c1: run, c2: chase
                dist = d1;
                bestClusterType = CircleClusterType::CHASE;
            } else {
                if (d1 < d2) {
                    // c1: run, c2: chase
                    dist = d1;
                    bestClusterType = CircleClusterType::CHASE;
                } else {
                    // c2: run, c1: chase
                    dist = d2;
                    bestClusterType = CircleClusterType::RUN;
                }
            }
            if (dist < 0.0) {
                continue;
            }
            if (dist < bestDist && !ClusterExistsInCurCircle(clusters, c1, c2)) {
                bestCluster = c2;
                bestDist = dist;
            }
        }

        if (bestCluster != nullptr) {
            if (bestClusterType == CircleClusterType::CHASE) {
                // c1: run, c2: chase
                // if 'bestCluster' already exists
                auto iter = pairs.find(bestCluster);
                if (iter == pairs.end()) {
                    pairs[bestCluster] = c1;
                } else {
                    const double oldDist = (iter->first->center - iter->second->center).norm();
                    if (bestDist < oldDist) {
                        iter->second = c1;
                    }
                }
            } else if (bestClusterType == CircleClusterType::RUN) {
                // c2: run, c1: chase
                pairs[c1] = bestCluster;
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

        auto iter = reversedPairs.find(type2);

        if (iter == reversedPairs.end()) {
            reversedPairs[type2] = type1;
        } else {
            CircleClusterInfo::Ptr prevType1 = iter->second;
            if ((type1->center - type2->center).squaredNorm() <
                (prevType1->center - type2->center).squaredNorm()) {
                iter->second = type1;
            }
        }
    }

    pairs.clear();
    for (auto& entry : reversedPairs) {
        pairs[entry.second] = entry.first;
    }

    // Next, we remove the singular relationships (non-inverse relationships) from the current
    // 'pairs' and 'reversedPairs'.
    for (auto it = pairs.begin(); it != pairs.end();) {
        auto it2 = pairs.find(it->second);
        if (it2 == pairs.cend()) {
            // it's reversible
            ++it;
            continue;
        }
        if (it2->second == it->first) {
            it = pairs.erase(it);
            continue;
        }

        const double d1 = (it->first->center - it->second->center).squaredNorm();
        const double d2 = (it2->first->center - it2->second->center).squaredNorm();
        if (d1 > d2) {
            it = pairs.erase(it);
            continue;
        } else {
            it = pairs.erase(it2);
            continue;
        }
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
                                                  double areaThd) {
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
