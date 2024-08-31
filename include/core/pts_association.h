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

#ifndef IKALIBR_PTS_ASSOCIATION_H
#define IKALIBR_PTS_ASSOCIATION_H

#include "config/configor.h"
#include "util/cloud_define.hpp"
#include "ufo/map/point_cloud.h"
#include "ufo/map/surfel_map.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

struct PointToSurfelCorr;
using PointToSurfelCorrPtr = std::shared_ptr<PointToSurfelCorr>;

struct PointToSurfelCondition {
    double pointToSurfelMax;
    std::uint8_t queryDepthMin;
    std::uint8_t queryDepthMax;
    std::size_t surfelPointMin;
    double planarityMin;

    explicit PointToSurfelCondition(
        double pointToSurfelMax = Configor::Prior::LiDARDataAssociate::PointToSurfelMax,
        std::uint8_t queryDepthMin = Configor::Prior::LiDARDataAssociate::QueryDepthMin,
        std::uint8_t queryDepthMax = Configor::Prior::LiDARDataAssociate::QueryDepthMax,
        std::size_t surfelPointMin = Configor::Prior::LiDARDataAssociate::SurfelPointMin,
        double planarityMin = Configor::Prior::LiDARDataAssociate::PlanarityMin);

    PointToSurfelCondition &WithPointToSurfelMax(double val);

    PointToSurfelCondition &WithQueryDepthMin(uint8_t val);

    PointToSurfelCondition &WithQueryDepthMax(uint8_t val);

    PointToSurfelCondition &WithSurfelPointMin(std::size_t val);

    PointToSurfelCondition &WithPlanarityMin(double val);
};

class PointToSurfelAssociator {
public:
    using Ptr = std::shared_ptr<PointToSurfelAssociator>;

protected:
    ufo::map::SurfelMap _smp;

public:
    explicit PointToSurfelAssociator(const IKalibrPointCloud::Ptr &mapInW,
                                     double resolution,
                                     std::uint8_t depth);

    static Ptr Create(const IKalibrPointCloud::Ptr &mapInW, double resolution, std::uint8_t depth);

    std::vector<PointToSurfelCorrPtr> Association(const IKalibrPointCloud::Ptr &mapCloud,
                                                  const IKalibrPointCloud::Ptr &rawCloud,
                                                  const PointToSurfelCondition &condition);

    static double SurfelScore(const ufo::map::SurfelMap &m, const ufo::map::Node &n);

    [[nodiscard]] const ufo::map::SurfelMap &GetSurfelMap() const;

protected:
    static double PointToSurfel(const ufo::map::SurfelMap::Surfel &s, const ufo::map::Point3 &p);

    static Eigen::Vector4d SurfelCoeffs(const ufo::map::SurfelMap::Surfel &s);

    template <typename PointType>
    void InsertCloudToSurfelMap(ufo::map::SurfelMap &map, pcl::PointCloud<PointType> &pclCloud) {
        int cloudSize = pclCloud.size();

        ufo::map::PointCloud ufoCloud;
        ufoCloud.resize(cloudSize);

#pragma omp parallel for num_threads(omp_get_max_threads()) default(none) \
    shared(ufoCloud, cloudSize, pclCloud)
        for (int i = 0; i < cloudSize; i++) {
            ufoCloud[i].x = pclCloud.points[i].x;
            ufoCloud[i].y = pclCloud.points[i].y;
            ufoCloud[i].z = pclCloud.points[i].z;
        }

        map.insertSurfelPoint(std::begin(ufoCloud), std::end(ufoCloud));
    }
};
}  // namespace ns_ikalibr
#endif  // IKALIBR_PTS_ASSOCIATION_H
