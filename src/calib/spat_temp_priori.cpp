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

#include "calib/spat_temp_priori.h"
#include "config/configor.h"
#include "util/status.hpp"
#include "util/utils_tpl.hpp"
#include "sensor/camera_data_loader.h"
#include "calib/estimator.h"
#include "calib/calib_param_manager.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
SpatialTemporalPriori::Ptr SpatialTemporalPriori::Create() {
    return std::make_shared<SpatialTemporalPriori>();
}

const std::map<SpatialTemporalPriori::FromTo, Sophus::SO3d>& SpatialTemporalPriori::GetExtriSO3()
    const {
    return SO3_Sen1ToSen2;
}

const std::map<SpatialTemporalPriori::FromTo, Eigen::Vector3d>& SpatialTemporalPriori::GetExtriPOS()
    const {
    return POS_Sen1InSen2;
}

const std::map<SpatialTemporalPriori::FromTo, double>& SpatialTemporalPriori::GetTimeOffset()
    const {
    return TO_Sen1ToSen2;
}

const std::map<std::string, double>& SpatialTemporalPriori::GetReadout() const {
    return RS_READOUT;
}

void SpatialTemporalPriori::CheckValidityWithConfigor() const {
    // check map if its ambiguous
    if (auto [res, p] = IsMapAmbiguous(this->SO3_Sen1ToSen2); res) {
        throw Status(Status::ERROR,
                     "extrinsic rotation priori of the sensor pair '{}' and '{}' is ambiguous!!! "
                     "Check the spatiotemporal priori config file!!!",
                     p.first, p.second);
    }
    if (auto [res, p] = IsMapAmbiguous(this->POS_Sen1InSen2); res) {
        throw Status(Status::ERROR,
                     "extrinsic transaction priori of the sensor pair '{}' and '{}' is "
                     "ambiguous!!! Check the spatiotemporal priori config file!!!",
                     p.first, p.second);
    }
    if (auto [res, p] = IsMapAmbiguous(this->TO_Sen1ToSen2); res) {
        throw Status(Status::ERROR,
                     "time offset priori of the sensor pair '{}' and '{}' is ambiguous!!! Check "
                     "the spatiotemporal priori config file!!!",
                     p.first, p.second);
    }

    std::set<std::string> topics;
    // add topics to vector
    for (const auto& [topic, _] : Configor::DataStream::IMUTopics) {
        topics.insert(topic);
    }
    for (const auto& [topic, _] : Configor::DataStream::RadarTopics) {
        topics.insert(topic);
    }
    for (const auto& [topic, _] : Configor::DataStream::LiDARTopics) {
        topics.insert(topic);
    }
    for (const auto& [topic, _] : Configor::DataStream::CameraTopics) {
        topics.insert(topic);
    }
    auto CheckTopic = [&topics](const std::pair<std::string, std::string>& sensorPair,
                                const std::string& prioriDesc) {
        const auto& [sen1, sen2] = sensorPair;
        if (sen1 == sen2) {
            throw Status(Status::WARNING,
                         "invalid (same topic names) prior {}: from sensor '{}' to sensor '{}'!!!  "
                         "Check the spatiotemporal priori config file!!!",
                         prioriDesc, sen1, sen2);
        }
        // this topic does not exist
        if (topics.count(sen1) == 0) {
            throw Status(Status::WARNING,
                         "invalid prior {}: sensor '{}' does not exist in Configor!!!  Check the "
                         "spatiotemporal priori config file!!!",
                         prioriDesc, sen1);
        } else if (topics.count(sen2) == 0) {
            throw Status(Status::WARNING,
                         "invalid prior {}: sensor '{}' does not exist in Configor!!!  Check the "
                         "spatiotemporal priori config file!!!",
                         prioriDesc, sen2);
        }
    };
    for (const auto& [sensorPair, _] : SO3_Sen1ToSen2) {
        CheckTopic(sensorPair, "extrinsic rotation");
    }
    for (const auto& [sensorPair, _] : POS_Sen1InSen2) {
        CheckTopic(sensorPair, "extrinsic translation");
    }
    for (const auto& [sensorPair, _] : TO_Sen1ToSen2) {
        CheckTopic(sensorPair, "time offset");
    }
    const double RT_PADDING = Configor::Prior::ReadoutTimePadding;
    for (const auto& [sensor, readout] : RS_READOUT) {
        auto iter = Configor::DataStream::CameraTopics.find(sensor);
        // this topic does not exist
        if (iter == Configor::DataStream::CameraTopics.cend()) {
            throw Status(Status::WARNING,
                         "invalid prior readout time: camera '{}' does not exist in Configor!!!  "
                         "Check the spatiotemporal priori config file!!!",
                         sensor);
        }
        auto model = EnumCast::stringToEnum<CameraModelType>(iter->second.Type);
        if (!IsOptionWith(CameraModelType::RS, model)) {
            // is not a rs camera
            throw Status(Status::ERROR,
                         "prior readout time is set for camera '{}', but it's not a RS camera in "
                         "Configor!!!  Check the spatiotemporal priori config file!!!",
                         sensor);
        }
        // range check
        if (readout > RT_PADDING) {
            throw Status(Status::ERROR,
                         "prior readout time ('{}') for rs camera '{}' is out of range ([0.000, "
                         "{:.3f}]), set a larger padding for readout time!",
                         readout, sensor, RT_PADDING);
        }
    }
}

void SpatialTemporalPriori::AddSpatTempPrioriConstraint(Estimator& estimator,
                                                        CalibParamManager& parMagr) const {
    // extrinsic rotations
    std::map<std::string, Sophus::SO3d*> SO3Address;
    std::map<std::string, Eigen::Vector3d*> POSAddress;
    // time offsets
    std::map<std::string, double*> TOAddress;
    // extrinsic rotations
    for (auto& [topic, item] : parMagr.EXTRI.SO3_BiToBr) {
        SO3Address.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.SO3_RjToBr) {
        SO3Address.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.SO3_LkToBr) {
        SO3Address.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.SO3_CmToBr) {
        SO3Address.insert({topic, &item});
    }
    // extrinsic translations
    for (auto& [topic, item] : parMagr.EXTRI.POS_BiInBr) {
        POSAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.POS_RjInBr) {
        POSAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.POS_LkInBr) {
        POSAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.POS_CmInBr) {
        POSAddress.insert({topic, &item});
    }
    // time offsets
    for (auto& [topic, item] : parMagr.TEMPORAL.TO_BiToBr) {
        TOAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.TEMPORAL.TO_RjToBr) {
        TOAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.TEMPORAL.TO_LkToBr) {
        TOAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.TEMPORAL.TO_CmToBr) {
        TOAddress.insert({topic, &item});
    }
    auto RefIMU = Configor::DataStream::ReferIMU;

    for (const auto& [sensorPair, Sen1ToSen2] : this->SO3_Sen1ToSen2) {
        const auto& [sen1, sen2] = sensorPair;
        Sophus::SO3d *rot1 = SO3Address.at(sen1), *rot2 = SO3Address.at(sen2);
        if (sen2 == RefIMU) {
            // if this priori is with respect to the reference IMU, set param constant directly
            *rot1 = Sen1ToSen2;
            if (estimator.HasParameterBlock(rot1->data())) {
                estimator.SetParameterBlockConstant(rot1->data());
            }
        } else if (estimator.HasParameterBlock(rot1->data()) ||
                   estimator.HasParameterBlock(rot2->data())) {
            // only one of the param block has been added to problem, we then add the constraint,
            // to make sure a unique least-squares solution
            estimator.AddPriorExtriSO3Constraint(Sen1ToSen2, rot1, rot2, PrioriWeight);
        }
    }
    for (const auto& [sensorPair, Sen1InSen2] : this->POS_Sen1InSen2) {
        const auto& [sen1, sen2] = sensorPair;
        Eigen::Vector3d *pos1 = POSAddress.at(sen1), *pos2 = POSAddress.at(sen2);
        Sophus::SO3d* rot2 = SO3Address.at(sen2);
        if (sen2 == RefIMU) {
            // if this priori is with respect to the reference IMU, set param constant directly
            *pos1 = Sen1InSen2;
            if (estimator.HasParameterBlock(pos1->data())) {
                estimator.SetParameterBlockConstant(pos1->data());
            }
        } else if (estimator.HasParameterBlock(pos1->data()) ||
                   estimator.HasParameterBlock(pos2->data())) {
            estimator.AddPriorExtriPOSConstraint(Sen1InSen2, pos1, rot2, pos2, PrioriWeight);
        }
    }
    for (const auto& [sensorPair, Sen1ToSen2] : this->TO_Sen1ToSen2) {
        const auto& [sen1, sen2] = sensorPair;
        double *to1 = TOAddress.at(sen1), *to2 = TOAddress.at(sen2);
        if (sen2 == RefIMU) {
            // if this priori is with respect to the reference IMU, set param constant directly
            *to1 = Sen1ToSen2;
            if (estimator.HasParameterBlock(to1)) {
                estimator.SetParameterBlockConstant(to1);
            }
        } else if (estimator.HasParameterBlock(to1) || estimator.HasParameterBlock(to2)) {
            estimator.AddPriorTimeOffsetConstraint(Sen1ToSen2, to1, to2, PrioriWeight);
        }
    }
    // readout times (we set them as constraints in optimization)
    for (const auto& [sensor, readout] : this->RS_READOUT) {
        double* data = &parMagr.TEMPORAL.RS_READOUT.at(sensor);
        *data = readout;
        if (estimator.HasParameterBlock(data)) {
            estimator.SetParameterBlockConstant(data);
        }
    }
    spdlog::info("add spatial and temp priori constraint finished");
}

}  // namespace ns_ikalibr