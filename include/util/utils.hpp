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

#ifndef IKALIBR_UTILS_HPP
#define IKALIBR_UTILS_HPP

#include "ctraj/utils/sophus_utils.hpp"
#include "spdlog/fmt/fmt.h"
#include "filesystem"
#include "ctraj/core/pose.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/bundled/color.h"
#include "magic_enum.hpp"
#include "ros/ros.h"
#include "regex"

namespace ns_ikalibr {

    // config the 'spdlog' log pattern
    inline void ConfigSpdlog() {
        // [log type]-[thread]-[time] message
        spdlog::set_pattern("%^[%L]%$-[%t]-[%H:%M:%S.%e] %v");

        // set log level
        spdlog::set_level(spdlog::level::debug);
    }

    inline void PrintIKalibrLibInfo() {
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
                     "+----------+----------------+--------+---------------------+" << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    template<class Type>
    Type GetParamFromROS(const std::string &param) {
        Type par;
        if (!ros::param::get(param, par)) {
            throw std::runtime_error("the ros param couldn't obtained from '" + param + "'.");
        }
        return par;
    }

    inline std::string UpperString(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](char c) { return std::toupper(c); });
        return s;
    }

    inline std::string LowerString(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](char c) { return std::tolower(c); });
        return s;
    }

    using namespace magic_enum::bitwise_operators;

    template<class EnumType>
    inline bool IsOptionWith(EnumType desired, EnumType current) {
        return (desired == (desired & current));
    }

    inline bool IsNotWhiteSpace(const int character) {
        return character != ' ' && character != '\n' && character != '\r' && character != '\t';
    }

    inline void StringLeftTrim(std::string *str) {
        str->erase(str->begin(), std::find_if(str->begin(), str->end(), IsNotWhiteSpace));
    }

    inline void StringRightTrim(std::string *str) {
        str->erase(std::find_if(str->rbegin(), str->rend(), IsNotWhiteSpace).base(), str->end());
    }

    inline void StringTrim(std::string *str) {
        StringLeftTrim(str);
        StringRightTrim(str);
    }

    inline std::string GetIndexedFilename(int idx, int num) {
        std::stringstream stream;
        std::string filename;
        stream << std::setfill('0') << std::setw(static_cast<int>(std::log10(num)) + 1) << idx;
        stream >> filename;
        return filename;
    }

    inline Eigen::MatrixXd TangentBasis(const Eigen::Vector3d &g0) {
        Eigen::Vector3d b, c;
        Eigen::Vector3d a = g0.normalized();
        Eigen::Vector3d tmp(0, 0, 1);
        if (a == tmp)
            tmp << 1, 0, 0;
        b = (tmp - a * (a.transpose() * tmp)).normalized();
        c = a.cross(b);
        Eigen::MatrixXd bc(3, 2);
        bc.block<3, 1>(0, 0) = b;
        bc.block<3, 1>(0, 1) = c;
        return bc;
    }

    inline std::vector<Sophus::SE3d> GenerateUniformPoseOnSphere(int n, double r) {
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

    /**
     * @brief sampling the samples without replacement
     *
     * @param num the num of the samples to sampling
     * @param engine the random engine
     * @param start the start index
     * @param end the end index
     * @param step the step
     * @attention range: [start, end](step) i.e. for [1, 5](2) -> pool: {1, 3, 5}
     * @return std::vector<std::size_t>
     */
    inline std::vector<std::size_t> SamplingWoutReplace(std::default_random_engine &engine,
                                                        std::size_t num,
                                                        std::size_t start,
                                                        std::size_t end,
                                                        std::size_t step = 1) {
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

    /**
     * @brief sampling the samples without replacement
     *
     * @tparam ElemType the element type
     * @param engine the random engine
     * @param dataVec the data vector
     * @param num the num of the samples to sampling
     * @return std::vector<std::size_t>
     */
    template<typename ElemType>
    std::vector<std::size_t> SamplingWoutReplace(std::default_random_engine &engine,
                                                 const std::vector<ElemType> &dataVec,
                                                 std::size_t num) {
        return SamplingWoutReplace(engine, num, 0, dataVec.size() - 1, 1);
    }

    /**
     * @brief sampling the samples without replacement
     *
     * @tparam ElemType the element type
     * @param engine the random engine
     * @param dataVec the data vector
     * @param num the num of the samples to sampling
     * @return std::vector<ElemType>
     */
    template<typename ElemType>
    std::vector<ElemType> SamplingWoutReplace2(std::default_random_engine &engine,
                                               const std::vector<ElemType> &dataVec,
                                               std::size_t num) {
        std::vector<std::size_t> res = SamplingWoutReplace(engine, dataVec, num);
        std::vector<ElemType> samples(num);
        for (int i = 0; i != static_cast<int>(num); ++i) {
            samples.at(i) = dataVec.at(res.at(i));
        }
        return samples;
    }

    /**
     * @brief sampling the samples with replacement
     *
     * @param engine the random engine
     * @param num the num of the samples to sampling
     * @param start the start index
     * @param end the end index
     * @param step the step
     * @attention range: [start, end](step) i.e. for [1, 5](2) -> pool: {1, 3, 5}
     * @return std::vector<std::size_t>
     */
    inline std::vector<std::size_t> SamplingWithReplace(std::default_random_engine &engine,
                                                        std::size_t num,
                                                        std::size_t start,
                                                        std::size_t end,
                                                        std::size_t step = 1) {
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

    /**
     * @brief sampling the samples with replacement
     *
     * @tparam ElemType the element type
     * @param engine the random engine
     * @param dataVec the data vector
     * @param num the num of the samples to sampling
     * @return std::vector<std::size_t>
     */
    template<typename ElemType>
    std::vector<std::size_t> SamplingWithReplace(std::default_random_engine &engine,
                                                 const std::vector<ElemType> &dataVec,
                                                 std::size_t num) {
        return SamplingWithReplace(engine, num, 0, dataVec.size() - 1, 1);
    }

    /**
     * @brief sampling the samples with replacement
     *
     * @tparam ElemType the element type
     * @param engine the random engine
     * @param dataVec the data vector
     * @param num the num of the samples to sampling
     * @return std::vector<ElemType>
     */
    template<typename ElemType>
    std::vector<ElemType> SamplingWithReplace2(std::default_random_engine &engine,
                                               const std::vector<ElemType> &dataVec,
                                               std::size_t num) {
        std::vector<std::size_t> res = SamplingWithReplace(engine, dataVec, num);
        std::vector<ElemType> samples(num);
        for (int i = 0; i != static_cast<int>(num); ++i) {
            samples.at(i) = dataVec.at(res.at(i));
        }
        return samples;
    }

    // obtain the rotation from aligned {w} to {ref} based on 'SO3_B0ToRef' and 'gravityInRef'
    inline Sophus::SO3d ObtainAlignedWtoRef(const Sophus::SO3d &SO3_B0ToRef, const Eigen::Vector3d &gravityInRef) {
        Eigen::Vector3d zNegAxis = -SO3_B0ToRef.matrix().col(2);
        Eigen::Vector3d rotDir = zNegAxis.cross(gravityInRef).normalized();
        double angRad = std::acos(zNegAxis.dot(gravityInRef) / gravityInRef.norm());
        Sophus::SO3d SO3_WtoRef = Sophus::SO3d(Eigen::AngleAxisd(angRad, rotDir).toRotationMatrix()) * SO3_B0ToRef;
        return SO3_WtoRef;
    }

    inline Eigen::Vector3d RotMatToYPR(const Eigen::Matrix3d &R) {
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

    inline double NormalizeAngle(double ang_degree) {
        if (ang_degree > 180.0) ang_degree -= 360.0;
        if (ang_degree < -180.0) ang_degree += 360.0;
        return ang_degree;
    }

    template<typename Derived>
    inline Eigen::Matrix<typename Derived::Scalar, 3, 3> SkewSymmetric(
            const Eigen::MatrixBase<Derived> &v3d) {
        /**
         *  0 -z  y
         *  z  0 -x
         * -y  x  0
         */
        Eigen::Matrix<typename Derived::Scalar, 3, 3> m;
        m << typename Derived::Scalar(0), -v3d.z(), v3d.y(), v3d.z(),
                typename Derived::Scalar(0), -v3d.x(), -v3d.y(), v3d.x(),
                typename Derived::Scalar(0);
        return m;
    }

    template<typename Derived>
    inline Eigen::Matrix<typename Derived::Scalar, 4, 4> LeftQuatMatrix(
            const Eigen::QuaternionBase<Derived> &q) {
        /**
         * p * q = left_mat(p) * q
         *
         *   1 -qz  qy  qx
         *  qz   1 -qx  qy
         * -qy  qx   1  qz
         * -qx -qy -qz  qw
         */
        Eigen::Matrix<typename Derived::Scalar, 4, 4> m;
        Eigen::Matrix<typename Derived::Scalar, 3, 1> vq = q.vec();
        typename Derived::Scalar q4 = q.w();
        m.block(0, 0, 3, 3) << q4 * Eigen::Matrix3d::Identity() + SkewSymmetric(vq);
        m.block(3, 0, 1, 3) << -vq.transpose();
        m.block(0, 3, 3, 1) << vq;
        m(3, 3) = q4;
        return m;
    }

    template<typename Derived>
    inline Eigen::Matrix<typename Derived::Scalar, 4, 4> RightQuatMatrix(
            const Eigen::QuaternionBase<Derived> &p) {
        /**
         * p * q = right_mat(q) * p
         *
         *   1  qz -qy  qx
         * -qz   1  qx  qy
         *  qy -qx   1  qz
         * -qx -qy -qz  qw
         */
        Eigen::Matrix<typename Derived::Scalar, 4, 4> m;
        Eigen::Matrix<typename Derived::Scalar, 3, 1> vp = p.vec();
        typename Derived::Scalar p4 = p.w();
        m.block(0, 0, 3, 3) << p4 * Eigen::Matrix3d::Identity() - SkewSymmetric(vp);
        m.block(3, 0, 1, 3) << -vp.transpose();
        m.block(0, 3, 3, 1) << vp;
        m(3, 3) = p4;
        return m;
    }

    struct RetrieveKey {
        template<typename T>
        typename T::first_type operator()(T keyValuePair) const {
            return keyValuePair.first;
        }
    };

    template<typename KeyType, typename ValueType>
    std::vector<KeyType> ExtractKeysAsVec(const std::map<KeyType, ValueType> &inputMap) {
        std::vector<KeyType> keys;
        std::transform(inputMap.begin(), inputMap.end(), std::back_inserter(keys), RetrieveKey());
        return keys;
    }

    template<typename KeyType, typename ValueType>
    std::set<KeyType> ExtractKeysAsSet(const std::map<KeyType, ValueType> &inputMap) {
        std::set<KeyType> keys;
        std::transform(inputMap.begin(), inputMap.end(), std::inserter(keys, keys.end()), RetrieveKey());
        return keys;
    }

    template<typename ScaleType>
    inline std::string FormatValueVector(const std::vector<const char *> &descVec,
                                         const std::vector<ScaleType> &valVec,
                                         const char *scaleFormatStr = "{:+011.6f}") {
        std::string str;
        const int M = static_cast<int>(descVec.size());
        for (int i = 0; i < (M - 1); ++i) {
            str += '\'' + std::string(descVec.at(i)) + "': " +
                   fmt::format(scaleFormatStr, valVec.at(i)) + ", ";
        }
        str += '\'' + std::string(descVec.at(M - 1)) + "': " +
               fmt::format(scaleFormatStr, valVec.at(M - 1));
        return str;
    }

    /**
     * @brief a function to get all the filenames in the directory
     * @param directory the directory
     * @return the filenames in the directory
     */
    inline std::vector<std::string> FilesInDir(const std::string &directory) {
        std::vector<std::string> files;
        for (const auto &elem: std::filesystem::directory_iterator(directory))
            if (elem.status().type() != std::filesystem::file_type::directory)
                files.emplace_back(std::filesystem::canonical(elem.path()).c_str());
        std::sort(files.begin(), files.end());
        return files;
    }

    inline std::vector<std::string> FilesInDirRecursive(const std::string &directory) {
        std::vector<std::string> files;
        for (const auto &elem: std::filesystem::recursive_directory_iterator(directory))
            if (elem.status().type() != std::filesystem::file_type::directory)
                files.emplace_back(std::filesystem::canonical(elem.path()).c_str());
        std::sort(files.begin(), files.end());
        return files;
    }

    /**
     * @brief a function to split a string to some string elements according the splitor
     * @param str the string to be split
     * @param splitor the splitor char
     * @param ignoreEmpty whether ignoring the empty string element or not
     * @return the split string vector
     */
    inline std::vector<std::string> SplitString(const std::string &str, char splitor, bool ignoreEmpty = true) {
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

    template<typename Scale, int Rows, int Cols>
    inline Eigen::Matrix<Scale, Rows, Cols>
    TrapIntegrationOnce(const std::vector<std::pair<Scale, Eigen::Matrix<Scale, Rows, Cols>>> &data) {
        Eigen::Matrix<Scale, Rows, Cols> sum = Eigen::Matrix<Scale, Rows, Cols>::Zero();
        for (int i = 0; i < static_cast<int>(data.size()) - 1; ++i) {
            int j = i + 1;
            const auto &di = data.at(i);
            const auto &dj = data.at(j);
            sum += (di.second + dj.second) * (dj.first - di.first) * Scale(0.5);
        }
        return sum;
    }

    template<typename Scale, int Rows, int Cols>
    inline Eigen::Matrix<Scale, Rows, Cols>
    TrapIntegrationTwice(const std::vector<std::pair<Scale, Eigen::Matrix<Scale, Rows, Cols>>> &data) {
        std::vector<std::pair<Scale, Eigen::Matrix<Scale, Rows, Cols>>> dataOnce;
        Eigen::Matrix<Scale, Rows, Cols> sum = Eigen::Matrix<Scale, Rows, Cols>::Zero();
        for (int i = 0; i < static_cast<int>(data.size()) - 1; ++i) {
            int j = i + 1;
            const auto &di = data.at(i);
            const auto &dj = data.at(j);
            sum += (di.second + dj.second) * (dj.first - di.first) * Scale(0.5);
            dataOnce.push_back({(dj.first + di.first) * Scale(0.5), sum});
        }
        return TrapIntegrationOnce(dataOnce);
    }

    template<typename EigenVectorType>
    inline auto EigenVecXToVector(const EigenVectorType &eigenVec) {
        std::vector<typename EigenVectorType::Scalar> vec(eigenVec.rows());
        for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
            vec.at(i) = eigenVec(i);
        }
        return vec;
    }

    template<class ScaleType>
    inline Sophus::SO3<ScaleType> ComputeKarcherMean(const std::vector<Sophus::SO3<ScaleType>> &so3Vec,
                                                     double tolerance = 1E-15) {
        if (so3Vec.empty()) {
            return {};
        }
        Sophus::SO3<ScaleType> X = so3Vec.front();
        while (true) {
            Eigen::Vector3<ScaleType> A = Eigen::Vector3<ScaleType>::Zero();
            for (const auto &item: so3Vec) {
                A += (X.inverse() * item).log();
            }
            A /= static_cast<double>(so3Vec.size());
            if (A.norm() < tolerance) {
                break;
            } else {
                X = X * Sophus::SO3<ScaleType>::exp(A);
            }
        }
        return X;
    }

    template<class Scale, int Rows, int Cols>
    inline Eigen::Matrix<Scale, Rows, Cols>
    ComputeMatVecMean(const std::vector<Eigen::Matrix<Scale, Rows, Cols>> &vec) {
        Eigen::Matrix<Scale, Rows, Cols> X = Eigen::Matrix<Scale, Rows, Cols>::Zero();
        for (const auto &item: vec) {
            X += item;
        }
        X /= static_cast<double>(vec.size());
        return X;
    }

    template<class Type>
    inline Type ComputeNumericalMean(const std::vector<Type> &vec) {
        Type X = static_cast<Type>(double{0.0});
        for (const auto &item: vec) {
            X += item;
        }
        X /= static_cast<double >(vec.size());
        return X;
    }

    template<class KeyType, class ValueType>
    inline std::vector<ValueType> ValueVecFromMap(const std::map<KeyType, ValueType> &m) {
        std::vector<ValueType> v;
        std::transform(m.begin(), m.end(), std::back_inserter(v), [](const std::pair<KeyType, ValueType> &p) {
            return p.second;
        });
        return v;
    }

#define CONCAT(a, b) CONCAT_INNER(a, b)
#define CONCAT_INNER(a, b) a##b
#define UNIQUE_NAME(base) CONCAT(base, __COUNTER__)

    inline bool _1_(const std::string a) {
        std::ifstream b(a);
        std::regex c(
                "\x41\x75\x74\x68\x6f\x72\x3a\x20\x53\x68\x75\x6f\x6c\x6f\x6e\x67\x20\x43\x68\x65\x6e\x20" "\\" "\x28\x73\x68\x6c\x63\x68\x65\x6e\x40\x77\x68\x75" "\\" "\x2e\x65\x64\x75" "\\" "\x2e\x63\x6e" "\\" "\x29"
        );
        std::regex d(
                "\x47\x69\x74\x48\x75\x62\x3a\x20\x68\x74\x74\x70\x73\x3a\x2f\x2f\x67\x69\x74\x68\x75\x62\x2e\x63\x6f\x6d\x2f\x55\x6e\x73\x69\x67\x6e\x65\x64\x2d\x4c\x6f\x6e\x67"
        );
        std::regex e(
                "\x4f\x52\x43\x49\x44\x3a\x20\x30\x30\x30\x30\x2d\x30\x30\x30\x32\x2d\x35\x32\x38\x33\x2d\x39\x30\x35\x37"
        );
        std::regex f(
                "\x68\x74\x74\x70\x73\x3a\x2f\x2f\x67\x69\x74\x68\x75\x62\x2e\x63\x6f\x6d\x2f\x55\x6e\x73\x69\x67\x6e\x65\x64\x2d\x4c\x6f\x6e\x67\x2f\x69\x4b\x61\x6c\x69\x62\x72\x2e\x67\x69\x74"
        );
        std::string g;
        int h = (0x1a8 + 6249 - 0x1a11), i = (0x545 + 245 - 0x63a);
        while (std::getline(b, g)) {
            if (std::regex_search(g, c)) { ++h; }
            if (std::regex_search(g, d)) { ++h; }
            if (std::regex_search(g, e)) { ++h; }
            if (std::regex_search(g, f)) { ++h; }
            ++i;
            if (h == (0x1ab + 384 - 0x327) || i > (0x70b + 7509 - 0x2456)) { break; }
        }
        b.close();
        if (h != (0x31a + 6724 - 0x1d5a)) {
            throw std::runtime_error(
                    "\x44\x6f\x20\x6e\x6f\x74\x20\x74\x61\x6d\x70\x65\x72\x20\x77\x69\x74\x68\x20\x74\x68\x65\x20\x43\x6f\x70\x79\x72\x69\x67\x68\x74\x20\x6f\x66\x20\x69\x4b\x61\x6c\x69\x62\x72\x20\x69\x6e\x20\x66\x69\x6c\x65\x20\x27"
                    + a + "\x27\x21\x21\x21"
            );
        }
        return true;
    }
// this is a unique tag of iKalibr, do not touch this!!! otherwise, unexpected exceptions would happen!!!
#define _3_ namespace {bool UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);}
    _3_
}

#endif //IKALIBR_UTILS_HPP
