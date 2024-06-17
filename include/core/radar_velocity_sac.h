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

#ifndef IKALIBR_RADAR_VELOCITY_SAC_H
#define IKALIBR_RADAR_VELOCITY_SAC_H

#include "sensor/radar.h"
#include "opengv/sac/SampleConsensusProblem.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
class RadarVelocitySacProblem : public opengv::sac::SampleConsensusProblem<Eigen::Vector3d> {
public:
    /** The model we are trying to fit (transformation) */
    typedef Eigen::Vector3d model_t;

public:
    /**
     * \brief Constructor.
     * \param[in] adapter Visitor holding bearing vectors, world points, etc.
     * \param[in] algorithm The algorithm we want to use.
     * \param[in] randomSeed Whether to seed the random number generator with
     *            the current time.
     */
    explicit RadarVelocitySacProblem(const RadarTargetArray::Ptr &data, bool randomSeed = true)
        : opengv::sac::SampleConsensusProblem<model_t>(randomSeed),
          _data(data) {
        setUniformIndices(static_cast<int>(_data->GetTargets().size()));
    };

    /**
     * Destructor.
     */
    ~RadarVelocitySacProblem() override = default;

    /**
     * \brief See parent-class.
     */
    bool computeModelCoefficients(const std::vector<int> &indices,
                                  model_t &outModel) const override;

    /**
     * \brief See parent-class.
     */
    void getSelectedDistancesToModel(const model_t &model,
                                     const std::vector<int> &indices,
                                     std::vector<double> &scores) const override;

    /**
     * \brief See parent-class.
     */
    void optimizeModelCoefficients(const std::vector<int> &inliers,
                                   const model_t &model,
                                   model_t &optimized_model) override;

    /**
     * \brief See parent-class.
     */
    [[nodiscard]] int getSampleSize() const override;

protected:
    /** The adapter holding all input data */
    const RadarTargetArray::Ptr &_data;
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_RADAR_VELOCITY_SAC_H
