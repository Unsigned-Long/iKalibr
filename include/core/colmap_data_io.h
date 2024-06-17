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

#ifndef IKALIBR_COLMAP_DATA_IO_H
#define IKALIBR_COLMAP_DATA_IO_H

#include "util/utils.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace Eigen {
typedef Eigen::Matrix<float, 3, 4> Matrix3x4f;
typedef Eigen::Matrix<double, 3, 4> Matrix3x4d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<uint8_t, 3, 1> Vector3ub;
typedef Eigen::Matrix<uint8_t, 4, 1> Vector4ub;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
}  // namespace Eigen

namespace ns_ikalibr {

struct ColMapDataIO {
public:
    ////////////////////////////////////////////////////////////////////////////////
    // Index types, determines the maximum number of objects.
    ////////////////////////////////////////////////////////////////////////////////

    // Unique identifier for cameras.
    typedef uint32_t camera_t;

    // Unique identifier for images.
    typedef uint32_t image_t;

    // Each image pair gets a unique ID, see `Database::ImagePairToPairId`.
    typedef uint64_t image_pair_t;

    // Index per image, i.e. determines maximum number of 2D points per image.
    typedef uint32_t point2D_t;

    // Unique identifier per added 3D point. Since we add many 3D points,
    // delete them, and possibly re-add them again, the maximum number of allowed
    // unique indices should be large.
    typedef uint64_t point3D_t;

    // Values for invalid identifiers or indices.
    static constexpr camera_t kInvalidCameraId = std::numeric_limits<camera_t>::max();
    static constexpr image_t kInvalidImageId = std::numeric_limits<image_t>::max();
    static constexpr image_pair_t kInvalidImagePairId = std::numeric_limits<image_pair_t>::max();
    static constexpr point2D_t kInvalidPoint2DIdx = std::numeric_limits<point2D_t>::max();
    static constexpr point3D_t kInvalidPoint3DId = std::numeric_limits<point3D_t>::max();

public:
    struct Camera {
    public:
        // The unique identifier of the camera. If the identifier is not specified
        // it is set to `kInvalidCameraId`.
        camera_t camera_id_{};

        // The identifier of the camera model. If the camera model is not specified
        // the identifier is `kInvalidCameraModelId`.
        int model_id_{};

        // The dimensions of the image, 0 if not initialized.
        size_t width_{};
        size_t height_{};

        // The focal length, principal point, and extra parameters. If the camera
        // model is not specified, this vector is empty.
        std::vector<double> params_;

    public:
        Camera() = default;
    };

    struct Point2D {
    public:
        // The image coordinates in pixels, starting at upper left corner with 0.
        Eigen::Vector2d xy_;

        // The identifier of the 3D point. If the 2D point is not part of a 3D point
        // track the identifier is `kInvalidPoint3DId` and `HasPoint3D() = false`.
        point3D_t point3D_id_{};

    public:
        Point2D() = default;
    };

    struct Image {
    public:
        // Identifier of the image, if not specified `kInvalidImageId`.
        image_t image_id_{};

        // The name of the image, i.e. the relative path.
        std::string name_;

        // The identifier of the associated camera. Note that multiple images might
        // share the same camera. If not specified `kInvalidCameraId`.
        camera_t camera_id_{};

        // The pose of the image, defined as the transformation from world to image.
        // QW, QX, QY, QZ
        Eigen::Vector4d qvec_;
        Eigen::Vector3d tvec_;

        // All image points, including points that are not part of a 3D point track.
        std::vector<Point2D> points2D_;

    public:
        Image() = default;

        void SetPoints2D(const std::vector<Eigen::Vector2d> &points);

        void SetPoint3DForPoint2D(point2D_t point2D_idx, point3D_t point3D_id);

        [[nodiscard]] Eigen::Quaterniond QuatWorldToImg() const;
    };

    // Track class stores all observations of a 3D point.
    struct TrackElement {
    public:
        // The image in which the track element is observed.
        image_t image_id;
        // The point in the image that the track element is observed.
        point2D_t point2D_idx;

    public:
        TrackElement() = default;
    };

    struct Point3D {
    public:
        // The 3D position of the point.
        Eigen::Vector3d xyz_;

        // The color of the point in the range [0, 255].
        Eigen::Vector3ub color_;

        // The mean reprojection error in pixels.
        double error_{};

        // The track of the point as a list of image observations.
        // class Track track_;
        std::vector<TrackElement> track_;

    public:
        Point3D() = default;
    };

public:
    static std::map<camera_t, Camera> ReadCamerasText(const std::string &path);

    static std::map<image_t, Image> ReadImagesText(const std::string &path);

    static std::map<point3D_t, Point3D> ReadPoints3DText(const std::string &path);
};
}  // namespace ns_ikalibr

#endif  // IKALIBR_COLMAP_DATA_IO_H
