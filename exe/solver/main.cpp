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
#include "calib/calib_solver.h"
#include "spdlog/fmt/bundled/color.h"
#include "calib/calib_solver_io.h"

_3_

int main(int argc, char **argv) {
    ros::init(argc, argv, "ikalibr_prog");
    try {
        ns_ikalibr::ConfigSpdlog();

        ns_ikalibr::PrintIKalibrLibInfo();

        // load settings
        auto configPath = ns_ikalibr::GetParamFromROS<std::string>("/ikalibr_prog/config_path");
        spdlog::info("loading configure from yaml file '{}'...", configPath);
        if (!std::filesystem::exists(configPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::CRITICAL, "configure file dose not exist: '{}'", configPath);
        }

        if (!ns_ikalibr::Configor::LoadConfigure(configPath)) {
            throw ns_ikalibr::Status(ns_ikalibr::Status::CRITICAL, "load configure file from '{}' failed!", configPath);
        } else {
            ns_ikalibr::Configor::PrintMainFields();
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }

        // create parameter manager based on loaded configure information
        auto paramMagr = ns_ikalibr::CalibParamManager::InitParamsFromConfigor();
        paramMagr->ShowParamStatus();

        // create data manager
        auto dataMagr = ns_ikalibr::CalibDataManager::Create();
        dataMagr->LoadCalibData();

        // pass parameter manager and data manager to solver for solving
        auto solver = ns_ikalibr::CalibSolver::Create(dataMagr, paramMagr);
        // the calibration results are stored in 'paramMagr'
        solver->Process();

        // solve finished, save calibration results (file type: JSON | YAML | XML | BINARY)
        auto filename = ns_ikalibr::Configor::DataStream::OutputPath + "/ikalibr_param" +
                        ns_ikalibr::Configor::GetFormatExtension();
        paramMagr->Save(filename, ns_ikalibr::Configor::Preference::DataIOFormat());

        auto solverIO = ns_ikalibr::CalibSolverIO::Create(solver);
        // update maps if exists, and save them if needed
        solverIO->UpdateMapsToSave();

        if (ns_ikalibr::Configor::Preference::OutputBSplines) { solverIO->SaveBSplines(); }

        if (ns_ikalibr::Configor::Preference::OutputHessianMat) { solverIO->SaveHessianMatrix(); }

        if (ns_ikalibr::Configor::Preference::OutputVisualLiDARCovisibility) { solverIO->VerifyVisualLiDARConsistency(); }

        if (ns_ikalibr::Configor::Preference::OutputVisualKinematics) { solverIO->SaveVisualKinematics(); }

        if (ns_ikalibr::Configor::Preference::OutputVisualKinematics) { solverIO->SaveVisualReprojectionError(); }

        if (ns_ikalibr::Configor::Preference::OutputColorizedMap) { solverIO->SaveVisualColorizedMap(); }

        if (ns_ikalibr::Configor::Preference::OutputAlignedInertialMes) { solverIO->SaveAlignedInertialMes(); }

        static const auto FStyle = fmt::emphasis::italic | fmt::fg(fmt::color::green);
        spdlog::info(fmt::format(FStyle, "solving and outputting finished!!! Everything is fine!!!"));

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