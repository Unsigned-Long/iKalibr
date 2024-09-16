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

#include "core/pts_association.h"
#include "factor/data_correspondence.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

// -------------------
// PointToSurfelCondition
// -------------------

PointToSurfelCondition::PointToSurfelCondition(double pointToSurfelMax,
                                               uint8_t queryDepthMin,
                                               uint8_t queryDepthMax,
                                               size_t surfelPointMin,
                                               double planarityMin)
    : pointToSurfelMax(pointToSurfelMax),
      queryDepthMin(queryDepthMin),
      queryDepthMax(queryDepthMax),
      surfelPointMin(surfelPointMin),
      planarityMin(planarityMin) {}

PointToSurfelCondition &PointToSurfelCondition::WithPointToSurfelMax(double val) {
    this->pointToSurfelMax = val;
    return *this;
}

PointToSurfelCondition &PointToSurfelCondition::WithQueryDepthMin(uint8_t val) {
    this->queryDepthMin = val;
    return *this;
}

PointToSurfelCondition &PointToSurfelCondition::WithQueryDepthMax(uint8_t val) {
    this->queryDepthMax = val;
    return *this;
}

PointToSurfelCondition &PointToSurfelCondition::WithSurfelPointMin(std::size_t val) {
    this->surfelPointMin = val;
    return *this;
}

PointToSurfelCondition &PointToSurfelCondition::WithPlanarityMin(double val) {
    this->planarityMin = val;
    return *this;
}

// -----------------------
// PointToSurfelAssociator
// -----------------------

PointToSurfelAssociator::PointToSurfelAssociator(const IKalibrPointCloud::Ptr &mapInW,
                                                 double resolution,
                                                 std::uint8_t depth) {
    _smp = ufo::map::SurfelMap(resolution, depth);
    InsertCloudToSurfelMap<IKalibrPoint>(_smp, *mapInW);
}

PointToSurfelAssociator::Ptr PointToSurfelAssociator::Create(const IKalibrPointCloud::Ptr &mapInW,
                                                             double resolution,
                                                             std::uint8_t depth) {
    return std::make_shared<PointToSurfelAssociator>(mapInW, resolution, depth);
}

double PointToSurfelAssociator::SurfelScore(const ufo::map::SurfelMap &m, const ufo::map::Node &n) {
    const auto &s = m.getSurfel(n);
    double score = s.getPlanarity();
    return score;
}

double PointToSurfelAssociator::PointToSurfel(const ufo::map::SurfelMap::Surfel &s,
                                              const ufo::map::Point3 &p) {
    return std::abs(SurfelCoeffs(s).dot(Eigen::Vector4d(p.x, p.y, p.z, 1.0)));
}

Eigen::Vector4d PointToSurfelAssociator::SurfelCoeffs(const ufo::map::SurfelMap::Surfel &s) {
    auto norm = s.getNormal();
    auto d = -norm.dot(s.getMean());
    return {norm.x, norm.y, norm.z, d};
}

const ufo::map::SurfelMap &PointToSurfelAssociator::GetSurfelMap() const { return _smp; }

std::vector<PointToSurfelCorr::Ptr> PointToSurfelAssociator::Association(
    const IKalibrPointCloud::Ptr &mapCloud,
    const IKalibrPointCloud::Ptr &rawCloud,
    const PointToSurfelCondition &condition) {
    if (mapCloud == nullptr || rawCloud == nullptr) {
        return {};
    }

    namespace ufopred = ufo::map::predicate;

    // get the width and height of this scan
    const int pts = static_cast<int>(rawCloud->size());

    std::vector<double> winScores(pts, -1.0);
    std::vector<ufo::map::Node> winNodes(pts, ufo::map::Node());

#pragma omp parallel for num_threads(omp_get_max_threads()) default(none) \
    shared(pts, mapCloud, condition, winNodes, winScores)
    for (int i = 0; i < pts; ++i) {
        const auto &mp = mapCloud->at(i);

        // nan point
        if (IS_POS_NAN(mp)) {
            continue;
        }

        // predicate
        auto pred = ufopred::HasSurfel()
                    // depth constraint
                    && ufopred::DepthMin(condition.queryDepthMin) &&
                    ufopred::DepthMax(condition.queryDepthMax)
                    // point num constraint
                    && ufopred::NumSurfelPointsMin(condition.surfelPointMin)
                    // planarity constraint
                    && ufopred::SurfelPlanarityMin(condition.planarityMin)
                    // geometry constraint
                    && ufopred::Contains(ufo::geometry::Point(mp.x, mp.y, mp.z));

        double winScore = -1.0;
        ufo::map::Node winNode;

        for (const auto &node : _smp.query(pred)) {
            double s = SurfelScore(_smp, node);
            if (winScore < 0.0 || s > winScore) {
                // this surfel is a good surfel, check point to surfel distance
                if (PointToSurfel(_smp.getSurfel(node), ufo::map::Point3(mp.x, mp.y, mp.z)) <
                    condition.pointToSurfelMax) {
                    winScore = s, winNode = node;
                }
            }
        }

        winScores.at(i) = winScore, winNodes.at(i) = winNode;
    }

    std::vector<PointToSurfelCorr::Ptr> corrs;
    corrs.reserve(pts);
    for (int i = 0; i < pts; ++i) {
        double winScore = winScores.at(i);
        // valid
        if (winScore > 0.0) {
            const auto &rp = rawCloud->at(i);
            const auto &mp = mapCloud->at(i);

            auto corr =
                PointToSurfelCorr::Create(rp.timestamp, Eigen::Vector3d(rp.x, rp.y, rp.z), winScore,
                                          SurfelCoeffs(_smp.getSurfel(winNodes.at(i))));

            corr->pInMap = Eigen::Vector3d(mp.x, mp.y, mp.z);
            corr->node = winNodes.at(i);

            corrs.push_back(corr);
        }
    }

    return corrs;
}
}  // namespace ns_ikalibr
