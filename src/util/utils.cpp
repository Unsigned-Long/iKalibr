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

#include "util/utils.h"
#include "filesystem"
#include "ctraj/core/pose.hpp"
#include "spdlog/spdlog.h"
#include "regex"
#include "random"
#include "opencv2/imgproc.hpp"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

void ConfigSpdlog() {
    // [log type]-[thread]-[time] message
    spdlog::set_pattern("%^[%L]%$-[%t]-[%H:%M:%S.%e] %v");

    // set log level
    spdlog::set_level(spdlog::level::debug);
}

void PrintIKalibrLibInfo() {
    std::cout << "+----------------------------------------------------------+\n"
                 "|     ██▓ ██ ▄█▀▄▄▄       ██▓     ██▓ ▄▄▄▄    ██▀███       |\n"
                 "|    ▓██▒ ██▄█▒▒████▄    ▓██▒    ▓██▒▓█████▄ ▓██ ▒ ██▒     |\n"
                 "|    ▒██▒▓███▄░▒██  ▀█▄  ▒██░    ▒██▒▒██▒ ▄██▓██ ░▄█ ▒     |\n"
                 "|    ░██░▓██ █▄░██▄▄▄▄██ ▒██░    ░██░▒██░█▀  ▒██▀▀█▄       |\n"
                 "|    ░██░▒██▒ █▄▓█   ▓██▒░██████▒░██░░▓█  ▀█▓░██▓ ▒██▒     |\n"
                 "|    ░▓  ▒ ▒▒ ▓▒▒▒   ▓▒█░░ ▒░▓  ░░▓  ░▒▓███▀▒░ ▒▓ ░▒▓░     |\n"
                 "|     ▒ ░░ ░▒ ▒░ ▒   ▒▒ ░░ ░ ▒  ░ ▒ ░▒░▒   ░   ░▒ ░ ▒░     |\n"
                 "|     ▒ ░░ ░░ ░  ░   ▒     ░ ░    ▒ ░ ░    ░   ░░   ░      |\n"
                 "|     ░  ░  ░        ░  ░    ░  ░ ░   ░         ░          |\n"
                 "|                                          ░               |\n"
                 "+----------+-----------------------------------------------+\n"
                 "|  iKalibr | https://github.com/Unsigned-Long/iKalibr.git  |\n"
                 "+----------+----------------+--------+---------------------+\n"
                 "|  Author  | Shuolong Chen  | E-Mail | shlchen@whu.edu.cn  |\n"
                 "+----------+----------------+--------+---------------------+"
              << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
}

std::string UpperString(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](char c) { return std::toupper(c); });
    return s;
}

std::string LowerString(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](char c) { return std::tolower(c); });
    return s;
}

bool IsNotWhiteSpace(int character) {
    return character != ' ' && character != '\n' && character != '\r' && character != '\t';
}

void StringLeftTrim(std::string *str) {
    str->erase(str->begin(), std::find_if(str->begin(), str->end(), IsNotWhiteSpace));
}

void StringRightTrim(std::string *str) {
    str->erase(std::find_if(str->rbegin(), str->rend(), IsNotWhiteSpace).base(), str->end());
}

void StringTrim(std::string *str) {
    StringLeftTrim(str);
    StringRightTrim(str);
}

std::string GetIndexedFilename(int idx, int num) {
    std::stringstream stream;
    std::string filename;
    stream << std::setfill('0') << std::setw(static_cast<int>(std::log10(num)) + 1) << idx;
    stream >> filename;
    return filename;
}

Eigen::MatrixXd TangentBasis(const Eigen::Vector3d &g0) {
    Eigen::Vector3d b, c;
    Eigen::Vector3d a = g0.normalized();
    Eigen::Vector3d tmp(0, 0, 1);
    if (a == tmp) tmp << 1, 0, 0;
    b = (tmp - a * (a.transpose() * tmp)).normalized();
    c = a.cross(b);
    Eigen::MatrixXd bc(3, 2);
    bc.block<3, 1>(0, 0) = b;
    bc.block<3, 1>(0, 1) = c;
    return bc;
}

std::vector<Sophus::SE3d> GenerateUniformPoseOnSphere(int n, double r) {
    const double phi = (std::sqrt(5.0) - 1.0) * 0.5;
    std::vector<Sophus::SE3d> poses;
    for (int i = 1; i < n + 1; ++i) {
        auto z = (2.0 * i - 1.0) / n - 1;
        auto x = std::sqrt(1.0 - z * z) * std::cos(2.0 * M_PI * i * phi);
        auto y = std::sqrt(1.0 - z * z) * std::sin(2.0 * M_PI * i * phi);
        Eigen::Vector3d POS_BiInRef(x * r, y * r, z * r);
        Eigen::Vector3d zAxis = -POS_BiInRef.normalized();
        Eigen::Vector3d xAxis = TangentBasis(zAxis).block<3, 1>(0, 0);
        Eigen::Vector3d yAxis = zAxis.cross(xAxis);
        Eigen::Matrix3d rotMat;
        rotMat.col(0) = xAxis;
        rotMat.col(1) = yAxis;
        rotMat.col(2) = zAxis;
        poses.emplace_back(Sophus::SO3d(rotMat), POS_BiInRef);
    }
    return poses;
}

std::vector<std::size_t> SamplingWoutReplace(std::default_random_engine &engine,
                                             std::size_t num,
                                             std::size_t start,
                                             std::size_t end,
                                             std::size_t step) {
    // create the pool for sampling
    std::vector<std::size_t> idxPool((end - start) / step + 1);
    for (int i = 0; i != static_cast<int>(idxPool.size()); ++i) {
        idxPool.at(i) = start + i * step;
    }
    std::vector<std::size_t> res(num);
    // the engine
    for (std::size_t i = 0; i != num; ++i) {
        // generate the random index
        std::uniform_int_distribution<std::size_t> ui(0, idxPool.size() - 1);
        std::size_t ridx = ui(engine);
        // record it
        res.at(i) = idxPool.at(ridx);
        // remove it
        idxPool.at(ridx) = idxPool.back();
        idxPool.pop_back();
    }
    return res;
}

std::vector<std::size_t> SamplingWithReplace(std::default_random_engine &engine,
                                             std::size_t num,
                                             std::size_t start,
                                             std::size_t end,
                                             std::size_t step) {
    // create the pool for sampling
    std::vector<std::size_t> idxPool((end - start) / step + 1);
    for (int i = 0; i != static_cast<int>(idxPool.size()); ++i) {
        idxPool.at(i) = start + i * step;
    }
    std::vector<std::size_t> res(num);
    // the engine
    std::uniform_int_distribution<std::size_t> ui(0, idxPool.size() - 1);
    for (std::size_t i = 0; i != num; ++i) {
        // generate the random index
        std::size_t ridx = ui(engine);
        // record it
        res.at(i) = idxPool.at(ridx);
    }
    return res;
}

Sophus::SO3d ObtainAlignedWtoRef(const Sophus::SO3d &SO3_B0ToRef,
                                 const Eigen::Vector3d &gravityInRef) {
    Eigen::Vector3d zNegAxis = -SO3_B0ToRef.matrix().col(2);
    Eigen::Vector3d rotDir = zNegAxis.cross(gravityInRef).normalized();
    double angRad = std::acos(zNegAxis.dot(gravityInRef) / gravityInRef.norm());
    Sophus::SO3d SO3_WtoRef =
        Sophus::SO3d(Eigen::AngleAxisd(angRad, rotDir).toRotationMatrix()) * SO3_B0ToRef;
    return SO3_WtoRef;
}

Eigen::Vector3d RotMatToYPR(const Eigen::Matrix3d &R) {
    Eigen::Vector3d n = R.col(0);
    Eigen::Vector3d o = R.col(1);
    Eigen::Vector3d a = R.col(2);

    Eigen::Vector3d ypr(3);
    double y = atan2(n(1), n(0));
    double p = atan2(-n(2), n(0) * cos(y) + n(1) * sin(y));
    double r = atan2(a(0) * sin(y) - a(1) * cos(y), -o(0) * sin(y) + o(1) * cos(y));

    ypr(0) = y;
    ypr(1) = p;
    ypr(2) = r;
    return ypr / M_PI * 180.0;
}

double NormalizeAngle(double ang_degree) {
    if (ang_degree > 180.0) ang_degree -= 360.0;
    if (ang_degree < -180.0) ang_degree += 360.0;
    return ang_degree;
}

std::vector<std::string> FilesInDir(const std::string &directory) {
    std::vector<std::string> files;
    for (const auto &elem : std::filesystem::directory_iterator(directory))
        if (elem.status().type() != std::filesystem::file_type::directory)
            files.emplace_back(std::filesystem::canonical(elem.path()).c_str());
    std::sort(files.begin(), files.end());
    return files;
}

std::vector<std::string> FilesInDirRecursive(const std::string &directory) {
    std::vector<std::string> files;
    for (const auto &elem : std::filesystem::recursive_directory_iterator(directory))
        if (elem.status().type() != std::filesystem::file_type::directory)
            files.emplace_back(std::filesystem::canonical(elem.path()).c_str());
    std::sort(files.begin(), files.end());
    return files;
}

std::vector<std::string> SplitString(const std::string &str, char splitor, bool ignoreEmpty) {
    std::vector<std::string> vec;
    auto iter = str.cbegin();
    while (true) {
        auto pos = std::find(iter, str.cend(), splitor);
        auto elem = std::string(iter, pos);
        if (!(elem.empty() && ignoreEmpty)) {
            vec.push_back(elem);
        }
        if (pos == str.cend()) {
            break;
        }
        iter = ++pos;
    }
    return vec;
}

bool _1_(const std::string &a) {
    std::ifstream b(a);
    std::regex c(
        "\x41\x75\x74\x68\x6f\x72\x3a\x20\x53\x68\x75\x6f\x6c\x6f\x6e\x67\x20\x43\x68\x65\x6e\x20"
        "\\"
        "\x28\x73\x68\x6c\x63\x68\x65\x6e\x40\x77\x68\x75"
        "\\"
        "\x2e\x65\x64\x75"
        "\\"
        "\x2e\x63\x6e"
        "\\"
        "\x29");
    std::regex d(
        "\x47\x69\x74\x48\x75\x62\x3a\x20\x68\x74\x74\x70\x73\x3a\x2f\x2f\x67\x69\x74\x68\x75\x62"
        "\x2e\x63\x6f\x6d\x2f\x55\x6e\x73\x69\x67\x6e\x65\x64\x2d\x4c\x6f\x6e\x67");
    std::regex e(
        "\x4f\x52\x43\x49\x44\x3a\x20\x30\x30\x30\x30\x2d\x30\x30\x30\x32\x2d\x35\x32\x38\x33\x2d"
        "\x39\x30\x35\x37");
    std::regex f(
        "\x68\x74\x74\x70\x73\x3a\x2f\x2f\x67\x69\x74\x68\x75\x62\x2e\x63\x6f\x6d\x2f\x55\x6e\x73"
        "\x69\x67\x6e\x65\x64\x2d\x4c\x6f\x6e\x67\x2f\x69\x4b\x61\x6c\x69\x62\x72\x2e\x67\x69\x74");
    std::string g;
    int h = (0x1a8 + 6249 - 0x1a11), i = (0x545 + 245 - 0x63a);
    while (std::getline(b, g)) {
        if (std::regex_search(g, c)) {
            ++h;
        }
        if (std::regex_search(g, d)) {
            ++h;
        }
        if (std::regex_search(g, e)) {
            ++h;
        }
        if (std::regex_search(g, f)) {
            ++h;
        }
        ++i;
        if (h == (0x1ab + 384 - 0x327) || i > (0x70b + 7509 - 0x2456)) {
            break;
        }
    }
    b.close();
    if (h != (0x31a + 6724 - 0x1d5a)) {
        throw std::runtime_error(
            "\x44\x6f\x20\x6e\x6f\x74\x20\x74\x61\x6d\x70\x65\x72\x20\x77\x69\x74\x68\x20\x74\x68"
            "\x65\x20\x43\x6f\x70\x79\x72\x69\x67\x68\x74\x20\x6f\x66\x20\x69\x4b\x61\x6c\x69\x62"
            "\x72\x20\x69\x6e\x20\x66\x69\x6c\x65\x20\x27" +
            a + "\x27\x21\x21\x21");
    }
    return true;
}

void DrawKeypointOnCVMat(cv::Mat &img,
                         const Eigen::Vector2d &feat,
                         bool withBox,
                         const cv::Scalar &color) {
    DrawKeypointOnCVMat(img, cv::Point2d(feat(0), feat(1)), withBox, color);
}

void DrawKeypointOnCVMat(cv::Mat &img,
                                     const cv::Point2d &feat,
                                     bool withBox,
                                     const cv::Scalar &color) {
    // square
    if (withBox) {
        cv::drawMarker(img, feat, color, cv::MarkerTypes::MARKER_SQUARE, 10, 1);
    }
    // key point
    cv::drawMarker(img, feat, color, cv::MarkerTypes::MARKER_SQUARE, 2, 2);
}

void DrawLineOnCVMat(cv::Mat &img,
                     const Eigen::Vector2d &p1,
                     const Eigen::Vector2d &p2,
                     const cv::Scalar &color) {
    DrawLineOnCVMat(img, cv::Point2d(p1(0), p1(1)), cv::Point2d(p2(0), p2(1)), color);
}

void DrawLineOnCVMat(cv::Mat &img,
                                 const cv::Point2d &p1,
                                 const cv::Point2d &p2,
                                 const cv::Scalar &color) {
    cv::line(img, p1, p2, color, 1, cv::LINE_AA);
}
void PutTextOnCVMat(cv::Mat &img,
                    const std::string &str,
                    const cv::Point2d &pt,
                    double xBias,
                    double yBias,
                    const cv::Scalar &color) {
    cv::putText(img, str, cv::Point2d(pt.x + xBias, pt.y + yBias),
                cv::HersheyFonts::FONT_HERSHEY_PLAIN, 1.0, color, 2);
}
void PutTextOnCVMat(cv::Mat &img,
                                const std::string &str,
                                const Eigen::Vector2d &pt,
                                double xBias,
                                double yBias,
                                const cv::Scalar &color) {
    PutTextOnCVMat(img, str, cv::Point2d(pt(0), pt(1)), xBias, yBias, color);
}
}  // namespace ns_ikalibr