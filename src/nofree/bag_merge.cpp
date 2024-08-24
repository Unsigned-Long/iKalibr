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

#include "nofree/bag_merge.h"
#include "rosbag/view.h"
#include "rosbag/bag.h"
#include "util/status.hpp"
#include "spdlog/spdlog.h"
#include "filesystem"
#include "cereal/types/map.hpp"
#include "cereal/types/string.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
std::string MergeConfigor::BagMergeInfo::GetDstTopic(const std::string &srcTopic) const {
    auto iter = topicsToMerge.find(srcTopic);
    return (iter == topicsToMerge.end() ? srcTopic : iter->second);
}

std::vector<std::string> MergeConfigor::BagMergeInfo::GetSrcTopicVec() const {
    std::vector<std::string> topics;
    std::transform(topicsToMerge.begin(), topicsToMerge.end(), std::back_inserter(topics),
                   [](const auto &p) { return p.first; });
    return topics;
}

MergeConfigor::Ptr MergeConfigor::Create() { return std::make_shared<MergeConfigor>(); }

MergeConfigor::Ptr MergeConfigor::LoadConfigure(const std::string &filename,
                                                CerealArchiveType::Enum archiveType) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return nullptr;
    }
    auto archive = GetInputArchiveVariant(file, archiveType);
    auto configor = MergeConfigor::Create();
    try {
        SerializeByInputArchiveVariant(archive, archiveType,
                                       cereal::make_nvp("MergeConfigor", *configor));
    } catch (const cereal::Exception &exception) {
        throw Status(Status::CRITICAL,
                     "The configuration file '{}' for 'MergeConfigor' is "
                     "outdated or broken, and can not be loaded in iKalibr using cereal!!! "
                     "To make it right, please refer to our latest configuration file "
                     "template released at "
                     "https://github.com/Unsigned-Long/iKalibr/blob/master/config/tool/"
                     "config-bag-merge.yaml, and then fix your custom configuration "
                     "file. Detailed cereal "
                     "exception information: \n'{}'",
                     filename, exception.what());
    }
    return configor;
}

bool MergeConfigor::SaveConfigure(const std::string &filename,
                                  CerealArchiveType::Enum archiveType) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    auto archive = GetOutputArchiveVariant(file, archiveType);
    SerializeByOutputArchiveVariant(archive, archiveType, cereal::make_nvp("MergeConfigor", *this));
    return true;
}

void MergeConfigor::PrintMainFields() {
    for (const auto &bag : _bags) {
        std::stringstream stream;
        for (const auto &[srcTopic, dstTopic] : bag.topicsToMerge) {
            stream << "{'" << srcTopic << "' -> '" << dstTopic << "'}";
        }
        spdlog::info("rosbag path: '{}', topics to merge in this bag: {}.", bag.bagPath,
                     stream.str());
    }
    spdlog::info("output rosbag path: '{}'.", _outputBagPath);
}

BagMerger::BagMerger(MergeConfigor::Ptr mergeConfigor)
    : _configor(std::move(mergeConfigor)) {}

BagMerger::Ptr BagMerger::Create(const MergeConfigor::Ptr &mergeConfigor) {
    return std::make_shared<BagMerger>(mergeConfigor);
}

std::pair<ros::Time, ros::Time> BagMerger::Process() {
    for (const auto &bagInfo : _configor->_bags) {
        if (!std::filesystem::exists(bagInfo.bagPath)) {
            throw Status(Status::CRITICAL, "the bag dose not exist: '{}'", bagInfo.bagPath);
        }
    }
    spdlog::info("start merge {} bag(s) to '{}'.", _configor->_bags.size(),
                 _configor->_outputBagPath);
    auto dstParentPath = std::filesystem::path(_configor->_outputBagPath).parent_path();
    if (!std::filesystem::exists(dstParentPath)) {
        if (!std::filesystem::create_directories(dstParentPath)) {
            throw Status(Status::CRITICAL,
                         "the parent path of the dst bag, i.e., '{}', dose not exist and can not "
                         "be created: '{}'",
                         dstParentPath.string());
        }
    }
    auto dstBag = std::make_unique<rosbag::Bag>();
    dstBag->open(_configor->_outputBagPath, rosbag::BagMode::Write | rosbag::BagMode::Read);
    for (const auto &bagInfo : _configor->_bags) {
        spdlog::info("process bag at '{}'...", bagInfo.bagPath);
        // open current rosbag
        auto srcBag = std::make_unique<rosbag::Bag>();
        srcBag->open(bagInfo.bagPath, rosbag::BagMode::Read);
        // query
        auto view = rosbag::View();
        if (auto topicVec = bagInfo.GetSrcTopicVec(); topicVec.empty()) {
            view.addQuery(*srcBag);
        } else {
            view.addQuery(*srcBag, rosbag::TopicQuery(bagInfo.GetSrcTopicVec()));
        }
        for (const auto &item : view) {
            dstBag->write(bagInfo.GetDstTopic(item.getTopic()), item.getTime(), item,
                          item.getConnectionHeader());
        }
        srcBag->close();
    }
    auto view = rosbag::View();
    view.addQuery(*dstBag);
    auto st = view.getBeginTime(), et = view.getEndTime();
    dstBag->close();
    spdlog::info("process finished...");
    return std::pair{st, et};
}
}  // namespace ns_ikalibr