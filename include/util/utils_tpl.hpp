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

#ifndef IKALIBR_UTILS_TPL_HPP
#define IKALIBR_UTILS_TPL_HPP

#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {
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

template <class EnumType>
bool IsOptionWith(EnumType desired, EnumType current) {
    return (desired == (desired & current));
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
                              const char *scaleFormatStr) {
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
                                          double tolerance) {
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

}  // namespace ns_ikalibr

#endif  // IKALIBR_UTILS_TPL_HPP
