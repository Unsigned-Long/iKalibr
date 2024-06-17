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

#ifndef IKALIBR_UFOMAP_LEARNER_HPP
#define IKALIBR_UFOMAP_LEARNER_HPP

#include "core/pts_association.h"
#include "viewer/viewer.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
struct UFOMapLearner {
public:
    static void Learn() {
        ns_ikalibr::Configor::LoadConfigure(
            "/home/csl/ros_ws/iKalibr/src/ikalibr/config/config.yaml");
        auto viewer = Viewer::Create(nullptr, nullptr);

        IKalibrPointCloud::Ptr cloud(new IKalibrPointCloud);
        pcl::io::loadPCDFile(
            "/home/csl/ros_ws/iKalibr/thirdparty/ctraj/thirdparty/tiny-viewer/data/scan.pcd",
            *cloud);

        auto associator = PointToSurfelAssociator::Create(
            cloud, Configor::Prior::LiDARDataAssociate::MapResolution,
            Configor::Prior::LiDARDataAssociate::MapDepthLevels);
        auto condition = PointToSurfelCondition();
        viewer->AddSurfelMap(associator->GetSurfelMap(), condition, Viewer::VIEW_ASSOCIATION);
        viewer->AddCloud(cloud, Viewer::VIEW_ASSOCIATION,
                         ns_viewer::Colour::Black().WithAlpha(0.2f), DefaultPointSize);
        std::cin.get();
    }
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_UFOMAP_LEARNER_HPP
