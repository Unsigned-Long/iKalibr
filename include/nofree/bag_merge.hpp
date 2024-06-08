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

#ifndef IKALIBR_BAG_MERGE_HPP
#define IKALIBR_BAG_MERGE_HPP

#include "rosbag/view.h"
#include "rosbag/bag.h"
#include "ros/ros.h"
#include "util/utils.hpp"
#include "util/status.hpp"
#include "spdlog/spdlog.h"
#include "util/cereal_archive_helper.hpp"

_3_

namespace ns_ikalibr {

    struct MergeConfigor {
    public:
        using Ptr = std::shared_ptr<MergeConfigor>;

        struct BagMergeInfo {
            std::string bagPath;
            std::map<std::string, std::string> topicsToMerge;

            [[nodiscard]] std::string GetDstTopic(const std::string &srcTopic) const {
                auto iter = topicsToMerge.find(srcTopic);
                return (iter == topicsToMerge.end() ? srcTopic : iter->second);
            }

            [[nodiscard]] std::vector<std::string> GetSrcTopicVec() const {
                std::vector<std::string> topics;
                std::transform(
                        topicsToMerge.begin(), topicsToMerge.end(), std::back_inserter(topics),
                        [](const auto &p) { return p.first; }
                );
                return topics;
            }

        public:
            template<class Archive>
            void serialize(Archive &ar) {
                ar(cereal::make_nvp("BagPath", bagPath), cereal::make_nvp("TopicsToMerge", topicsToMerge));
            }
        };

    public:
        std::vector<BagMergeInfo> _bags;
        std::string _outputBagPath;

        MergeConfigor() = default;

        static Ptr Create() {
            return std::make_shared<MergeConfigor>();
        }

        // load configure information from file
        static Ptr LoadConfigure(const std::string &filename, CerealArchiveType::Enum archiveType) {
            std::ifstream file(filename);
            if (!file.is_open()) { return nullptr; }
            auto archive = GetInputArchiveVariant(file, archiveType);
            auto configor = MergeConfigor::Create();
            SerializeByInputArchiveVariant(archive, archiveType, cereal::make_nvp("MergeConfigor", *configor));
            return configor;
        }

        // save configure information to file
        bool SaveConfigure(const std::string &filename, CerealArchiveType::Enum archiveType) {
            std::ofstream file(filename);
            if (!file.is_open()) { return false; }
            auto archive = GetOutputArchiveVariant(file, archiveType);
            SerializeByOutputArchiveVariant(archive, archiveType, cereal::make_nvp("MergeConfigor", *this));
            return true;
        }

        void PrintMainFields() {
            for (const auto &bag: _bags) {
                std::stringstream stream;
                for (const auto &[srcTopic, dstTopic]: bag.topicsToMerge) {
                    stream << "{'" << srcTopic << "' -> '" << dstTopic << "'}";
                }
                spdlog::info("rosbag path: '{}', topics to merge in this bag: {}.", bag.bagPath, stream.str());
            }
            spdlog::info("output rosbag path: '{}'.", _outputBagPath);
        }

    public:
        template<class Archive>
        void serialize(Archive &ar) {
            ar(cereal::make_nvp("Bags", _bags), cereal::make_nvp("OutputBagPath", _outputBagPath));
        }
    };

    class BagMerger {
    public:
        using Ptr = std::shared_ptr<BagMerger>;

    private:
        MergeConfigor::Ptr _configor;
    public:
        explicit BagMerger(MergeConfigor::Ptr mergeConfigor) : _configor(std::move(mergeConfigor)) {}

        static Ptr Create(const MergeConfigor::Ptr &mergeConfigor) {
            return std::make_shared<BagMerger>(mergeConfigor);
        }

        void Process() {
            spdlog::info("start merge {} bag(s) to '{}'.", _configor->_bags.size(), _configor->_outputBagPath);
            auto dstBag = std::make_unique<rosbag::Bag>();
            dstBag->open(_configor->_outputBagPath, rosbag::BagMode::Write);
            for (const auto &bagInfo: _configor->_bags) {
                spdlog::info("process bag at '{}'...", bagInfo.bagPath);
                // open current rosbag
                auto srcBag = std::make_unique<rosbag::Bag>();
                srcBag->open(bagInfo.bagPath, rosbag::BagMode::Read);
                // query
                auto view = rosbag::View();
                view.addQuery(*srcBag, rosbag::TopicQuery(bagInfo.GetSrcTopicVec()));
                for (const auto &item: view) {
                    dstBag->write(
                            bagInfo.GetDstTopic(item.getTopic()), item.getTime(), item, item.getConnectionHeader()
                    );
                }
                srcBag->close();
            }
            dstBag->close();
            spdlog::info("process finished...");
        }
    };
}

#endif //IKALIBR_BAG_MERGE_HPP
