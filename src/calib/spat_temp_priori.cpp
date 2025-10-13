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

std::optional<Eigen::Vector3d> SpatialTemporalPriori::GetGravity() const {
    if (GRAVITY == Eigen::Vector3d::Zero())
        return {};
    return GRAVITY;
}

const std::map<std::string, double>& SpatialTemporalPriori::GetMinVisualScale() const {
    return MIN_VISUAL_SCALE;
}

bool SpatialTemporalPriori::HasSO3ToBr(const std::string& topic) const {
    return hasSO3ToBr.find(topic) != hasSO3ToBr.end();
}

std::optional<Sophus::SO3d> SpatialTemporalPriori::GetSO3ToBr(const std::string& topic) const {
    if (!HasSO3ToBr(topic))
        return {};

    const auto& refImu = Configor::DataStream::ReferIMU;

    for (const auto& [fromTo, so3] : this->SO3_Sen1ToSen2) {
        const auto& [from, to] = fromTo;
        if (to == refImu && from == topic) {
            return so3;
        } else if (from == refImu && to == topic) {
            return so3.inverse();
        }
    }
    return {};
}

bool SpatialTemporalPriori::HasPosInBr(const std::string& topic) const {
    return hasPosToBr.find(topic) != hasPosToBr.end();
}

std::optional<Eigen::Vector3d> SpatialTemporalPriori::GetPosInBr(const std::string& topic) const {
    if (!HasPosInBr(topic))
        return {};

    const auto& refImu = Configor::DataStream::ReferIMU;

    for (const auto& [fromTo, pos] : this->POS_Sen1InSen2) {
        const auto& [from, to] = fromTo;
        if (to == refImu && from == topic) {
            return pos;
        } else if (from == refImu && to == topic) {
            // we know the so3 exists because of the HasPosToBr() check above
            const auto so3 = GetSO3ToBr(to);
            return (*so3) * (-pos);
        }
    }
    return {};
}

bool SpatialTemporalPriori::HasTOToBr(const std::string& topic) const {
    return hasTOToBr.find(topic) != hasTOToBr.end();
}

std::optional<double> SpatialTemporalPriori::GetTOToBr(const std::string& topic) const {
    if (!HasTOToBr(topic))
        return {};

    const auto& refImu = Configor::DataStream::ReferIMU;

    for (const auto& [fromTo, offset] : this->TO_Sen1ToSen2) {
        const auto& [from, to] = fromTo;
        if (to == refImu && from == topic) {
            return offset;
        } else if (from == refImu && to == topic) {
            return -offset;
        }
    }
    return {};
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

    // topic, camera type string
    std::map<std::string, std::string> optCamModelType;
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
    for (const auto& [topic, config] : Configor::DataStream::CameraTopics) {
        topics.insert(topic);
        optCamModelType.insert({topic, config.Type});
    }
    for (const auto& [topic, config] : Configor::DataStream::RGBDTopics) {
        topics.insert(topic);
        optCamModelType.insert({topic, config.Type});
    }
    for (const auto& [topic, _] : Configor::DataStream::EventTopics) {
        topics.insert(topic);
        // for event, rs exposure mode dose not make sense
        // optCamModelType.insert({topic, config.Type});
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
        }
        if (topics.count(sen2) == 0) {
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
        auto iter = optCamModelType.find(sensor);
        // this topic does not exist
        if (iter == optCamModelType.cend()) {
            throw Status(Status::WARNING,
                         "invalid prior readout time: camera '{}' does not exist in Configor!!!  "
                         "Check the spatiotemporal priori config file!!!",
                         sensor);
        }
        auto model = EnumCast::stringToEnum<CameraModelType>(iter->second);
        if (!IsOptionWith(CameraModelType::RS, model)) {
            // is not a rs camera
            throw Status(Status::ERROR,
                         "prior readout time is set for camera '{}', but it's not a RS camera in "
                         "Configor!!! Check the spatiotemporal priori config file!!!",
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

    if (GRAVITY != Eigen::Vector3d::Zero()) {
        if (std::abs(GRAVITY.norm() - Configor::Prior::GravityNorm) > 1e-3) {
            throw Status(Status::ERROR, "the given prior gravity vector [{}, {}, {}] does not have "
                         "norm equal to Prior::GravityNorm ({})! The vector's norm is: {}.",
                         GRAVITY.x(), GRAVITY.y(), GRAVITY.z(), Configor::Prior::GravityNorm,
                         GRAVITY.norm());
        }
    }

    for (const auto& [cam, _] : MIN_VISUAL_SCALE) {
        if (optCamModelType.count(cam) == 0) {
            throw Status(Status::ERROR, "MIN_VISUAL_SCALE defined for topic '{}' which is not an "
                         "optical camera topic!", cam);
        }
    }

    for (const auto& [topic, weight] : INTRI_WEIGHTS) {
        if (optCamModelType.count(topic) == 0 &&
            Configor::DataStream::EventTopics.count(topic) == 0 &&
            Configor::DataStream::IMUTopics.count(topic) == 0) {
            throw Status(Status::ERROR, "INTRI_WEIGHTS defined for topic '{}' which is neither IMU "
                         "nor camera topic!", topic);
        }
        if (weight < 0) {
            throw Status(Status::ERROR, "INTRI_WEIGHTS defined for topic '{}' is negative!");
        }
    }

    const auto& refImu = Configor::DataStream::ReferIMU;

    for (const auto& [fromTo, _] : SO3_Sen1ToSen2) {
        const auto& [from, to] = fromTo;
        if (from == refImu)
            hasSO3ToBr.insert(to);
        else if (to == refImu)
            hasSO3ToBr.insert(from);
    }

    for (const auto& [fromTo, _] : POS_Sen1InSen2) {
        const auto& [from, to] = fromTo;
        if (from == refImu)
            hasPosToBr.insert(to);
        else if (to == refImu && hasSO3ToBr.find(from) != hasSO3ToBr.end())
            hasPosToBr.insert(from);
    }

    for (const auto& [fromTo, _] : TO_Sen1ToSen2) {
        const auto& [from, to] = fromTo;
        if (from == refImu)
            hasTOToBr.insert(to);
        else if (to == refImu)
            hasTOToBr.insert(from);
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
    for (auto& [topic, item] : parMagr.EXTRI.SO3_DnToBr) {
        SO3Address.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.SO3_EsToBr) {
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
    for (auto& [topic, item] : parMagr.EXTRI.POS_DnInBr) {
        POSAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.EXTRI.POS_EsInBr) {
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
    for (auto& [topic, item] : parMagr.TEMPORAL.TO_DnToBr) {
        TOAddress.insert({topic, &item});
    }
    for (auto& [topic, item] : parMagr.TEMPORAL.TO_EsToBr) {
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
            if (sen1 == RefIMU && estimator.HasParameterBlock(rot1->data())) {
                *rot1 = Sophus::SO3d();
                estimator.SetParameterBlockConstant(rot1->data());
            }
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
            if (sen1 == RefIMU && estimator.HasParameterBlock(pos1->data())) {
                *pos1 = Eigen::Vector3d::Zero();
                estimator.SetParameterBlockConstant(pos1->data());
            }
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
            if (sen1 == RefIMU && estimator.HasParameterBlock(to1)) {
                *to1 = 0.0;
                estimator.SetParameterBlockConstant(to1);
            }
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
    const auto gravityPrior = this->GetGravity();
    if (gravityPrior) {
        auto gravity = &parMagr.GRAVITY;
        *gravity = *gravityPrior;
        if (estimator.HasParameterBlock(gravity->data())) {
            estimator.SetParameterBlockConstant(gravity->data());
        }
    }
    for (const auto& [topic, weight] : this->INTRI_WEIGHTS) {
        if (parMagr.INTRI.IMU.count(topic) > 0) {
            const auto& intri = parMagr.INTRI.IMU.at(topic);
            const auto& prioriIntri = parMagr.INTRI.PrioriIMU.at(topic);
            if (std::isinf(weight)) {
                // set params constant
                if (estimator.HasParameterBlock(intri->ACCE.BIAS.data())) {
                    estimator.SetParameterBlockConstant(intri->ACCE.BIAS.data());
                }
                if (estimator.HasParameterBlock(intri->ACCE.MAP_COEFF.data())) {
                    estimator.SetParameterBlockConstant(intri->ACCE.MAP_COEFF.data());
                }
                if (estimator.HasParameterBlock(intri->GYRO.BIAS.data())) {
                    estimator.SetParameterBlockConstant(intri->GYRO.BIAS.data());
                }
                if (estimator.HasParameterBlock(intri->GYRO.MAP_COEFF.data())) {
                    estimator.SetParameterBlockConstant(intri->GYRO.MAP_COEFF.data());
                }
                if (estimator.HasParameterBlock(intri->SO3_AtoG.data())) {
                    estimator.SetParameterBlockConstant(intri->SO3_AtoG.data());
                }
            } else if (std::isfinite(weight)) {
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->ACCE.BIAS, intri->ACCE.BIAS, weight);
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->ACCE.MAP_COEFF, intri->ACCE.MAP_COEFF, weight);
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->GYRO.BIAS, intri->GYRO.BIAS, weight);
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->GYRO.MAP_COEFF, intri->GYRO.MAP_COEFF, weight);
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->SO3_AtoG, intri->SO3_AtoG, weight);
            }
        }
        if (parMagr.INTRI.Camera.count(topic) > 0 || parMagr.INTRI.RGBD.count(topic) > 0) {
            const auto isRGBD = parMagr.INTRI.RGBD.count(topic) > 0;
            const auto& intri = isRGBD ?
                parMagr.INTRI.RGBD.at(topic)->intri : parMagr.INTRI.Camera.at(topic);
            const auto& prioriIntri = isRGBD ?
                parMagr.INTRI.PrioriRGBD.at(topic)->intri : parMagr.INTRI.PrioriCamera.at(topic);
            if (std::isinf(weight)) {
                // set params constant
                if (estimator.HasParameterBlock(intri->FXAddress())) {
                    estimator.SetParameterBlockConstant(intri->FXAddress());
                }
                if (estimator.HasParameterBlock(intri->FYAddress())) {
                    estimator.SetParameterBlockConstant(intri->FYAddress());
                }
                if (estimator.HasParameterBlock(intri->CXAddress())) {
                    estimator.SetParameterBlockConstant(intri->CXAddress());
                }
                if (estimator.HasParameterBlock(intri->CYAddress())) {
                    estimator.SetParameterBlockConstant(intri->CYAddress());
                }
            } else if (std::isfinite(weight)) {
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->FXAddress(), intri->FXAddress(), weight);
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->FYAddress(), intri->FYAddress(), weight);
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->CXAddress(), intri->CXAddress(), weight);
                estimator.AddPriorEqualityConstraint(
                    prioriIntri->CYAddress(), intri->CYAddress(), weight);
            }
        }
    }
    spdlog::info("add spatial and temp priori constraint finished");
}

}  // namespace ns_ikalibr