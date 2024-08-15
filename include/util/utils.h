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

#ifndef IKALIBR_UTILS_H
#define IKALIBR_UTILS_H
/**
 * @attention Description about "veta/type_def.hpp"
 'EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST' is a macro used to handle Eigen objects
within STL containers. In C++, we often use standard library containers (e.g., std::vector) to store
data. However, when we attempt to store Eigen objects (e.g., Eigen::Vector2d) in these containers,
we encounter some issues. When using Eigen objects in STL containers, special handling is required
to ensure proper memory alignment and initialization.
'EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST' exists to address such issues.
 The macro allows us to use Eigen objects in STL containers without manually dealing with memory
alignment and initialization. With this macro, we can use std::vector<Eigen::Vector2d> just like any
other vector, without worrying about the underlying implementation details. This specialization
needs to be defined before using code snippets like std::vector<Vector2d>. The benefit is that you
don’t need to declare std::vector with Eigen::aligned_allocator everywhere in your code. However,
the downside is that the specialization must be defined before using code snippets like
std::vector<Vector2d.
 */
#include "veta/type_def.hpp"
#include "magic_enum.hpp"
#include "ros/ros.h"
#include "ctraj/utils/sophus_utils.hpp"
#include "spdlog/fmt/fmt.h"
#include "opencv2/core.hpp"

namespace ns_ikalibr {

// config the 'spdlog' log pattern
void ConfigSpdlog();

void PrintIKalibrLibInfo();

// given n points and a x value, compute the y value using lagrange polynomial
template <class Type, int N>
double LagrangePolynomial(Type xQuery,
                          const std::array<Type, N> &xData,
                          const std::array<Type, N> &yData) {
    Type y = 0.0;
    for (int i = 0; i < N; ++i) {
        Type li = 1.0;
        for (int j = 0; j < N; ++j) {
            if (j == i) {
                continue;
            }
            li *= (xQuery - xData[j]) / (xData[i] - xData[j]);
        }
        y += yData[i] * li;
    }
    return y;
}

// given three points, compute the first order of the middle point using lagrange polynomial
template <class Type>
Type LagrangePolynomialTripleMidFOD(const std::array<Type, 3> &xData,
                                    const std::array<Type, 3> &yData) {
    Type x0 = xData[0];
    Type x1 = xData[1];
    Type x2 = xData[2];
    Type y0 = yData[0];
    Type y1 = yData[1];
    Type y2 = yData[2];
    Type v1 = y0 * (x1 - x2) / (x0 - x1) / (x0 - x2);
    Type v2 = y1 * (1.0 / (x1 - x0) + 1.0 / (x1 - x2));
    Type v3 = y2 * (x1 - x0) / (x2 - x0) / (x2 - x1);
    return v1 + v2 + v3;
}

// given three points, compute the first order of the middle point using lagrange polynomial
template <class Type>
Type LagrangePolynomialTripleMidFOD(const std::array<std::pair<Type, Type>, 3> &data) {
    std::array<Type, 3> xData;
    std::array<Type, 3> yData;
    for (int i = 0; i < 3; ++i) {
        xData[i] = data[i].first;
        yData[i] = data[i].second;
    }
    return LagrangePolynomialTripleMidFOD(xData, yData);
}

template <class Type>
Type GetParamFromROS(const std::string &param) {
    Type par;
    if (!ros::param::get(param, par)) {
        throw std::runtime_error("the ros param couldn't obtained from '" + param + "'.");
    }
    return par;
}

std::string UpperString(std::string s);

std::string LowerString(std::string s);

void DrawKeypointOnCVMat(cv::Mat &img,
                         const Eigen::Vector2d &feat,
                         bool withBox = true,
                         const cv::Scalar &color = cv::Scalar(0, 255, 0));

void DrawKeypointOnCVMat(cv::Mat &img,
                         const cv::Point2d &feat,
                         bool withBox = true,
                         const cv::Scalar &color = cv::Scalar(0, 255, 0));

void DrawLineOnCVMat(cv::Mat &img,
                     const Eigen::Vector2d &p1,
                     const Eigen::Vector2d &p2,
                     const cv::Scalar &color = cv::Scalar(0, 255, 0));

void DrawLineOnCVMat(cv::Mat &img,
                     const cv::Point2d &p1,
                     const cv::Point2d &p2,
                     const cv::Scalar &color = cv::Scalar(0, 255, 0));

void PutTextOnCVMat(cv::Mat &img,
                    const std::string &str,
                    const cv::Point2d &pt,
                    double xBias = 10.0,
                    double yBias = 0.0,
                    const cv::Scalar &color = cv::Scalar(255, 0, 0));

void PutTextOnCVMat(cv::Mat &img,
                    const std::string &str,
                    const Eigen::Vector2d &pt,
                    double xBias = 10.0,
                    double yBias = 0.0,
                    const cv::Scalar &color = cv::Scalar(255, 0, 0));

using namespace magic_enum::bitwise_operators;

template <class EnumType>
bool IsOptionWith(EnumType desired, EnumType current) {
    return (desired == (desired & current));
}

bool IsNotWhiteSpace(int character);

void StringLeftTrim(std::string *str);

void StringRightTrim(std::string *str);

void StringTrim(std::string *str);

std::string GetIndexedFilename(int idx, int num);

Eigen::MatrixXd TangentBasis(const Eigen::Vector3d &g0);

std::vector<Sophus::SE3d> GenerateUniformPoseOnSphere(int n, double r);

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
std::vector<std::size_t> SamplingWoutReplace(std::default_random_engine &engine,
                                             std::size_t num,
                                             std::size_t start,
                                             std::size_t end,
                                             std::size_t step = 1);

/**
 * @brief sampling the samples without replacement
 *
 * @tparam ElemType the element type
 * @param engine the random engine
 * @param dataVec the data vector
 * @param num the num of the samples to sampling
 * @return std::vector<std::size_t>
 */
template <typename ElemType>
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
template <typename ElemType>
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
std::vector<std::size_t> SamplingWithReplace(std::default_random_engine &engine,
                                             std::size_t num,
                                             std::size_t start,
                                             std::size_t end,
                                             std::size_t step = 1);

/**
 * @brief sampling the samples with replacement
 *
 * @tparam ElemType the element type
 * @param engine the random engine
 * @param dataVec the data vector
 * @param num the num of the samples to sampling
 * @return std::vector<std::size_t>
 */
template <typename ElemType>
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
template <typename ElemType>
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
Sophus::SO3d ObtainAlignedWtoRef(const Sophus::SO3d &SO3_B0ToRef,
                                 const Eigen::Vector3d &gravityInRef);

Eigen::Vector3d RotMatToYPR(const Eigen::Matrix3d &R);

double NormalizeAngle(double ang_degree);

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 3> SkewSymmetric(const Eigen::MatrixBase<Derived> &v3d) {
    /**
     *  0 -z  y
     *  z  0 -x
     * -y  x  0
     */
    Eigen::Matrix<typename Derived::Scalar, 3, 3> m;
    m << typename Derived::Scalar(0), -v3d.z(), v3d.y(), v3d.z(), typename Derived::Scalar(0),
        -v3d.x(), -v3d.y(), v3d.x(), typename Derived::Scalar(0);
    return m;
}

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar, 4, 4> LeftQuatMatrix(
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

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar, 4, 4> RightQuatMatrix(
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
    template <typename T>
    typename T::first_type operator()(T keyValuePair) const {
        return keyValuePair.first;
    }
};
struct RetrieveVal {
    template <typename T>
    typename T::first_type operator()(T keyValuePair) const {
        return keyValuePair.second;
    }
};

template <typename KeyType, typename ValueType>
std::vector<KeyType> ExtractKeysAsVec(const std::map<KeyType, ValueType> &inputMap) {
    std::vector<KeyType> keys;
    std::transform(inputMap.begin(), inputMap.end(), std::back_inserter(keys), RetrieveKey());
    return keys;
}

template <typename KeyType, typename ValueType>
std::set<KeyType> ExtractKeysAsSet(const std::map<KeyType, ValueType> &inputMap) {
    std::set<KeyType> keys;
    std::transform(inputMap.begin(), inputMap.end(), std::inserter(keys, keys.end()),
                   RetrieveKey());
    return keys;
}

template <typename KeyType, typename ValueType>
std::vector<ValueType> ExtractValsAsVec(const std::map<KeyType, ValueType> &inputMap) {
    std::vector<ValueType> vals;
    std::transform(inputMap.begin(), inputMap.end(), std::back_inserter(vals), RetrieveVal());
    return vals;
}

template <typename KeyType, typename ValueType>
std::set<ValueType> ExtractValsAsSet(const std::map<KeyType, ValueType> &inputMap) {
    std::set<ValueType> vals;
    std::transform(inputMap.begin(), inputMap.end(), std::inserter(vals, vals.end()),
                   RetrieveVal());
    return vals;
}

template <typename ScaleType>
std::string FormatValueVector(const std::vector<const char *> &descVec,
                              const std::vector<ScaleType> &valVec,
                              const char *scaleFormatStr = "{:+011.6f}") {
    std::string str;
    const int M = static_cast<int>(descVec.size());
    for (int i = 0; i < (M - 1); ++i) {
        str += '\'' + std::string(descVec.at(i)) +
               "': " + fmt::format(scaleFormatStr, valVec.at(i)) + ", ";
    }
    str += '\'' + std::string(descVec.at(M - 1)) +
           "': " + fmt::format(scaleFormatStr, valVec.at(M - 1));
    return str;
}

/**
 * @brief a function to get all the filenames in the directory
 * @param directory the directory
 * @return the filenames in the directory
 */
std::vector<std::string> FilesInDir(const std::string &directory);

std::vector<std::string> FilesInDirRecursive(const std::string &directory);

/**
 * @brief a function to split a string to some string elements according the splitor
 * @param str the string to be split
 * @param splitor the splitor char
 * @param ignoreEmpty whether ignoring the empty string element or not
 * @return the split string vector
 */
std::vector<std::string> SplitString(const std::string &str, char splitor, bool ignoreEmpty = true);

template <typename Scale, int Rows, int Cols>
Eigen::Matrix<Scale, Rows, Cols> TrapIntegrationOnce(
    const std::vector<std::pair<Scale, Eigen::Matrix<Scale, Rows, Cols>>> &data) {
    Eigen::Matrix<Scale, Rows, Cols> sum = Eigen::Matrix<Scale, Rows, Cols>::Zero();
    for (int i = 0; i < static_cast<int>(data.size()) - 1; ++i) {
        int j = i + 1;
        const auto &di = data.at(i);
        const auto &dj = data.at(j);
        sum += (di.second + dj.second) * (dj.first - di.first) * Scale(0.5);
    }
    return sum;
}

template <typename Scale, int Rows, int Cols>
Eigen::Matrix<Scale, Rows, Cols> TrapIntegrationTwice(
    const std::vector<std::pair<Scale, Eigen::Matrix<Scale, Rows, Cols>>> &data) {
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

template <typename EigenVectorType>
auto EigenVecXToVector(const EigenVectorType &eigenVec) {
    std::vector<typename EigenVectorType::Scalar> vec(eigenVec.rows());
    for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
        vec.at(i) = eigenVec(i);
    }
    return vec;
}

template <class ScaleType>
Sophus::SO3<ScaleType> ComputeKarcherMean(const std::vector<Sophus::SO3<ScaleType>> &so3Vec,
                                          double tolerance = 1E-15) {
    if (so3Vec.empty()) {
        return {};
    }
    Sophus::SO3<ScaleType> X = so3Vec.front();
    while (true) {
        Eigen::Vector3<ScaleType> A = Eigen::Vector3<ScaleType>::Zero();
        for (const auto &item : so3Vec) {
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

template <class Scale, int Rows, int Cols>
Eigen::Matrix<Scale, Rows, Cols> ComputeMatVecMean(
    const std::vector<Eigen::Matrix<Scale, Rows, Cols>> &vec) {
    Eigen::Matrix<Scale, Rows, Cols> X = Eigen::Matrix<Scale, Rows, Cols>::Zero();
    for (const auto &item : vec) {
        X += item;
    }
    X /= static_cast<double>(vec.size());
    return X;
}

template <class Type>
Type ComputeNumericalMean(const std::vector<Type> &vec) {
    Type X = static_cast<Type>(double{0.0});
    for (const auto &item : vec) {
        X += item;
    }
    X /= static_cast<double>(vec.size());
    return X;
}

template <class KeyType, class ValueType>
std::vector<ValueType> ValueVecFromMap(const std::map<KeyType, ValueType> &m) {
    std::vector<ValueType> v;
    std::transform(m.begin(), m.end(), std::back_inserter(v),
                   [](const std::pair<KeyType, ValueType> &p) { return p.second; });
    return v;
}

#define IKALIBR_CONCAT(a, b) IKALIBR_CONCAT_INNER(a, b)
#define IKALIBR_CONCAT_INNER(a, b) a##b
#define IKALIBR_UNIQUE_NAME(base) IKALIBR_CONCAT(base, __COUNTER__)

bool _1_(const std::string &a);
// this is a unique tag of iKalibr, do not touch this!!! otherwise, unexpected exceptions would
// happen!!!
#define _3_                                                    \
    namespace {                                                \
    bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__); \
    }
_3_
}  // namespace ns_ikalibr

#endif  // IKALIBR_UTILS_H
