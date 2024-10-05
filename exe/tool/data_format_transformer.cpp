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

#include "ros/ros.h"
#include "spdlog/spdlog.h"
#include "config/configor.h"
#include "util/status.hpp"
#include "util/utils_tpl.hpp"
#include "calib/calib_param_manager.h"
#include "ctraj/core/pose.hpp"
#include "ctraj/core/spline_bundle.h"
#include "cereal/types/utility.hpp"
#include "cereal/types/list.hpp"
#include "spdlog/fmt/bundled/color.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_data_format_transformer");
    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // load parameters
        auto preDir =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_data_format_transformer/pre_dir");
        spdlog::info("pre directory: '{}'", preDir);

        auto srcFormatStr =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_data_format_transformer/src_format");
        spdlog::info("source data format: '{}'", srcFormatStr);
        ns_ikalibr::CerealArchiveType::Enum srcFormat;
        try {
            srcFormat = ns_ikalibr::EnumCast::stringToEnum<ns_ikalibr::CerealArchiveType::Enum>(
                srcFormatStr);
        } catch (...) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::CRITICAL,
                                     "data format '{}' is unsupported!!!", srcFormatStr);
        }
        const auto &srcExt = ns_ikalibr::Configor::Preference::FileExtension.at(srcFormat);

        auto dstFormatStr =
            ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_data_format_transformer/dst_format");
        spdlog::info("destination data format: '{}'", dstFormatStr);
        ns_ikalibr::CerealArchiveType::Enum dstFormat;
        try {
            dstFormat = ns_ikalibr::EnumCast::stringToEnum<ns_ikalibr::CerealArchiveType::Enum>(
                dstFormatStr);
        } catch (...) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::CRITICAL,
                                     "data format '{}' is unsupported!!!", dstFormatStr);
        }
        const auto &dstExt = ns_ikalibr::Configor::Preference::FileExtension.at(dstFormat);

        auto wsDirVecSrc = ns_ikalibr::GetParamFromROS<std::vector<std::string>>(
            "/ikalibr_data_format_transformer/ws_dir_vec");

        for (const auto &dir : wsDirVecSrc) {
            auto ws = preDir + '/' + dir;
            spdlog::info("perform data format transform for '{}', from '{}' to '{}'", ws, srcExt,
                         dstExt);

            // -------------------------
            // spatiotemporal parameters
            // -------------------------
            auto rName = ws + "/ikalibr_param" + srcExt;
            auto wName = ws + "/ikalibr_param" + dstExt;
            if (std::filesystem::exists(rName)) {
                spdlog::info("perform transformation:\n   '{}'\n-> '{}'", rName, wName);
                ns_ikalibr::CalibParamManager::Load(rName, srcFormat)->Save(wName, dstFormat);
            }

            // -------
            // splines
            // -------
            rName = ws + "/splines/knots" + srcExt;
            wName = ws + "/splines/knots" + dstExt;
            if (std::filesystem::exists(rName)) {
                spdlog::info("perform transformation:\n   '{}'\n-> '{}'", rName, wName);
                auto bundle =
                    ns_ctraj::SplineBundle<ns_ikalibr::Configor::Prior::SplineOrder>::Create({});
                {
                    // load
                    std::ifstream file(rName);
                    auto ar = ns_ikalibr::GetInputArchiveVariant(file, srcFormat);
                    SerializeByInputArchiveVariant(ar, srcFormat,
                                                   cereal::make_nvp("splines", *bundle));
                }
                {
                    // output
                    std::ofstream file(wName);
                    auto ar = ns_ikalibr::GetOutputArchiveVariant(file, dstFormat);
                    SerializeByOutputArchiveVariant(ar, dstFormat,
                                                    cereal::make_nvp("splines", *bundle));
                }
            }

            rName = ws + "/splines/samples" + srcExt;
            wName = ws + "/splines/samples" + dstExt;
            if (std::filesystem::exists(rName)) {
                spdlog::info("perform transformation:\n   '{}'\n-> '{}'", rName, wName);
                std::vector<ns_ctraj::Posed> poseSeq;
                {
                    // load
                    std::ifstream file(rName);
                    auto ar = ns_ikalibr::GetInputArchiveVariant(file, srcFormat);
                    SerializeByInputArchiveVariant(ar, srcFormat,
                                                   cereal::make_nvp("pose_seq", poseSeq));
                }
                {
                    // output
                    std::ofstream file(wName);
                    auto ar = ns_ikalibr::GetOutputArchiveVariant(file, dstFormat);
                    SerializeByOutputArchiveVariant(ar, dstFormat,
                                                    cereal::make_nvp("pose_seq", poseSeq));
                }
            }

            // -------
            // hessian
            // -------
            rName = ws + "/hessian/hessian" + srcExt;
            wName = ws + "/hessian/hessian" + dstExt;
            if (std::filesystem::exists(rName)) {
                spdlog::info("perform transformation:\n   '{}'\n-> '{}'", rName, wName);
                int row, col;
                std::vector<std::pair<std::string, int>> parOrderSize;
                Eigen::MatrixXd hessian;
                // load
                {
                    std::ifstream file(rName);
                    auto ar = ns_ikalibr::GetInputArchiveVariant(file, srcFormat);
                    SerializeByInputArchiveVariant(
                        ar, srcFormat, cereal::make_nvp("row", row), cereal::make_nvp("col", col),
                        cereal::make_nvp("par_order_size", parOrderSize));
                    hessian.resize(row, col);
                    SerializeByInputArchiveVariant(ar, srcFormat,
                                                   cereal::make_nvp("hessian", hessian));
                }
                {
                    // output
                    std::ofstream file(wName);
                    auto ar = ns_ikalibr::GetOutputArchiveVariant(file, dstFormat);
                    SerializeByOutputArchiveVariant(
                        ar, dstFormat, cereal::make_nvp("row", row), cereal::make_nvp("col", col),
                        cereal::make_nvp("hessian", hessian),
                        cereal::make_nvp("par_order_size", parOrderSize));
                }
            }

            // ----------------------------
            // parameters in each iteration
            // ----------------------------
            {
                auto paramsEpochDir = ws + "/iteration/epoch";
                if (std::filesystem::exists(paramsEpochDir)) {
                    auto files = ns_ikalibr::FilesInDir(paramsEpochDir);
                    files.erase(std::remove_if(files.begin(), files.end(),
                                               [&srcExt](const std::string &str) {
                                                   return std::filesystem::path(str).extension() !=
                                                          srcExt;
                                               }),
                                files.end());
                    for (const auto &filename : files) {
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        ns_ikalibr::CalibParamManager::Load(filename, srcFormat)
                            ->Save(wName, dstFormat);
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    }
                }
            }

            {
                auto paramsStageDir = ws + "/iteration/stage";
                if (std::filesystem::exists(paramsStageDir)) {
                    auto files = ns_ikalibr::FilesInDir(paramsStageDir);
                    files.erase(std::remove_if(files.begin(), files.end(),
                                               [&srcExt](const std::string &str) {
                                                   return std::filesystem::path(str).extension() !=
                                                          srcExt;
                                               }),
                                files.end());
                    for (const auto &filename : files) {
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        ns_ikalibr::CalibParamManager::Load(filename, srcFormat)
                            ->Save(wName, dstFormat);
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    }
                }
            }

            // ---------
            // residuals
            // ---------

            // inertial errors
            auto subWS = ws + "/residuals/inertial_error";
            if (std::filesystem::exists(subWS)) {
                auto files = ns_ikalibr::FilesInDirRecursive(subWS);
                files.erase(std::remove_if(files.begin(), files.end(),
                                           [&srcExt](const std::string &str) {
                                               return std::filesystem::path(str).extension() !=
                                                      srcExt;
                                           }),
                            files.end());
                for (const auto &filename : files) {
                    if (std::filesystem::path(filename).filename() == "inertial_mes" + srcExt) {
                        std::list<ns_ikalibr::IMUFrame> rawMes, estMes, diff;
                        // load
                        std::ifstream ifile(filename);
                        auto iar = GetInputArchiveVariant(ifile, srcFormat);
                        SerializeByInputArchiveVariant(iar, srcFormat,
                                                       cereal::make_nvp("raw_inertial", rawMes),
                                                       cereal::make_nvp("est_inertial", estMes),
                                                       cereal::make_nvp("inertial_diff", diff));
                        // save
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        std::ofstream ofile(wName);
                        auto oar = GetOutputArchiveVariant(ofile, dstFormat);
                        SerializeByOutputArchiveVariant(oar, dstFormat,
                                                        cereal::make_nvp("raw_inertial", rawMes),
                                                        cereal::make_nvp("est_inertial", estMes),
                                                        cereal::make_nvp("inertial_diff", diff));
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    } else if (std::filesystem::path(filename).filename() ==
                               "aligned_mes_to_ref" + srcExt) {
                        std::list<ns_ikalibr::IMUFrame> estMes;
                        // load
                        std::ifstream ifile(filename);
                        auto iar = GetInputArchiveVariant(ifile, srcFormat);
                        SerializeByInputArchiveVariant(
                            iar, srcFormat, cereal::make_nvp("aligned_inertial", estMes));
                        // save
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        std::ofstream ofile(wName);
                        auto oar = GetOutputArchiveVariant(ofile, dstFormat);
                        SerializeByOutputArchiveVariant(
                            oar, dstFormat, cereal::make_nvp("aligned_inertial", estMes));
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    }
                }
            }

            // reprojection error
            subWS = ws + "/residuals/reproj_error";
            if (std::filesystem::exists(subWS)) {
                auto files = ns_ikalibr::FilesInDirRecursive(subWS);
                files.erase(std::remove_if(files.begin(), files.end(),
                                           [&srcExt](const std::string &str) {
                                               return std::filesystem::path(str).extension() !=
                                                      srcExt;
                                           }),
                            files.end());
                for (const auto &filename : files) {
                    if (std::filesystem::path(filename).filename() == "residuals" + srcExt) {
                        std::list<Eigen::Vector2d> reprojErrors;
                        // load
                        std::ifstream ifile(filename);
                        auto iar = GetInputArchiveVariant(ifile, srcFormat);
                        SerializeByInputArchiveVariant(
                            iar, srcFormat, cereal::make_nvp("reproj_errors", reprojErrors));
                        // save
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        std::ofstream ofile(wName);
                        auto oar = GetOutputArchiveVariant(ofile, dstFormat);
                        SerializeByOutputArchiveVariant(
                            oar, dstFormat, cereal::make_nvp("reproj_errors", reprojErrors));
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    }
                }
            }

            // optical flow error
            subWS = ws + "/residuals/optical_flow_error";
            if (std::filesystem::exists(subWS)) {
                auto files = ns_ikalibr::FilesInDirRecursive(subWS);
                files.erase(std::remove_if(files.begin(), files.end(),
                                           [&srcExt](const std::string &str) {
                                               return std::filesystem::path(str).extension() !=
                                                      srcExt;
                                           }),
                            files.end());
                for (const auto &filename : files) {
                    if (std::filesystem::path(filename).filename() == "residuals" + srcExt) {
                        std::list<Eigen::Vector2d> velErrors;
                        // load
                        std::ifstream ifile(filename);
                        auto iar = GetInputArchiveVariant(ifile, srcFormat);
                        SerializeByInputArchiveVariant(iar, srcFormat,
                                                       cereal::make_nvp("of_errors", velErrors));
                        // save
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        std::ofstream ofile(wName);
                        auto oar = GetOutputArchiveVariant(ofile, dstFormat);
                        SerializeByOutputArchiveVariant(oar, dstFormat,
                                                        cereal::make_nvp("of_errors", velErrors));
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    }
                }
            }

            // radar doppler error
            subWS = ws + "/residuals/doppler_error";
            if (std::filesystem::exists(subWS)) {
                auto files = ns_ikalibr::FilesInDirRecursive(subWS);
                files.erase(std::remove_if(files.begin(), files.end(),
                                           [&srcExt](const std::string &str) {
                                               return std::filesystem::path(str).extension() !=
                                                      srcExt;
                                           }),
                            files.end());
                for (const auto &filename : files) {
                    if (std::filesystem::path(filename).filename() == "residuals" + srcExt) {
                        std::list<double> dopplerErrors;
                        // load
                        std::ifstream ifile(filename);
                        auto iar = GetInputArchiveVariant(ifile, srcFormat);
                        SerializeByInputArchiveVariant(
                            iar, srcFormat, cereal::make_nvp("doppler_errors", dopplerErrors));
                        // save
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        std::ofstream ofile(wName);
                        auto oar = GetOutputArchiveVariant(ofile, dstFormat);
                        SerializeByOutputArchiveVariant(
                            oar, dstFormat, cereal::make_nvp("doppler_errors", dopplerErrors));
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    }
                }
            }

            // lidar point-to-surfel error
            subWS = ws + "/residuals/lidar_pts_error";
            if (std::filesystem::exists(subWS)) {
                auto files = ns_ikalibr::FilesInDirRecursive(subWS);
                files.erase(std::remove_if(files.begin(), files.end(),
                                           [&srcExt](const std::string &str) {
                                               return std::filesystem::path(str).extension() !=
                                                      srcExt;
                                           }),
                            files.end());
                for (const auto &filename : files) {
                    if (std::filesystem::path(filename).filename() == "residuals" + srcExt) {
                        std::list<double> ptsErrors;
                        // load
                        std::ifstream ifile(filename);
                        auto iar = GetInputArchiveVariant(ifile, srcFormat);
                        SerializeByInputArchiveVariant(iar, srcFormat,
                                                       cereal::make_nvp("pts_errors", ptsErrors));
                        // save
                        wName = std::filesystem::path(filename).replace_extension(dstExt).string();
                        std::ofstream ofile(wName);
                        auto oar = GetOutputArchiveVariant(ofile, dstFormat);
                        SerializeByOutputArchiveVariant(oar, dstFormat,
                                                        cereal::make_nvp("pts_errors", ptsErrors));
                        spdlog::info("perform transformation:\n   '{}'\n-> '{}'", filename, wName);
                    }
                }
            }
        }

    } catch (const ns_ikalibr::IKalibrStatus &status) {
        // if error happened, print it
        static const auto FStyle = fmt::emphasis::italic | fmt::fg(fmt::color::green);
        static const auto WECStyle = fmt::emphasis::italic | fmt::fg(fmt::color::red);
        switch (status.flag) {
            case ns_ikalibr::Status::FINE:
                // this case usually won't happen
                spdlog::info(fmt::format(FStyle, "{}", status.what));
                break;
            case ns_ikalibr::Status::WARNING:
                spdlog::warn(fmt::format(WECStyle, "{}", status.what));
                break;
            case ns_ikalibr::Status::ERROR:
                spdlog::error(fmt::format(WECStyle, "{}", status.what));
                break;
            case ns_ikalibr::Status::CRITICAL:
                spdlog::critical(fmt::format(WECStyle, "{}", status.what));
                break;
        }
    } catch (const std::exception &e) {
        // an unknown exception not thrown by this program
        static const auto WECStyle = fmt::emphasis::italic | fmt::fg(fmt::color::red);
        spdlog::critical(fmt::format(WECStyle, "unknown error happened: '{}'", e.what()));
    }

    ros::shutdown();
    return 0;
}