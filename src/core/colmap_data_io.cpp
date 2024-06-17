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

#include "core/colmap_data_io.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

// ------------
// ColMapDataIO
// ------------

std::map<ColMapDataIO::camera_t, ColMapDataIO::Camera> ColMapDataIO::ReadCamerasText(
    const std::string &path) {
    std::map<camera_t, Camera> cameras_;

    std::ifstream file(path);

    std::string line;
    std::string item;

    while (std::getline(file, line)) {
        StringTrim(&line);

        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::stringstream line_stream(line);

        Camera camera;

        // ID
        std::getline(line_stream, item, ' ');
        camera.camera_id_ = std::stoul(item);
        // camera.SetCameraId(std::stoul(item));

        // MODEL
        std::getline(line_stream, item, ' ');
        // camera.SetModelIdFromName(item);

        // WIDTH
        std::getline(line_stream, item, ' ');
        camera.width_ = std::stoll(item);
        // camera.SetWidth(std::stoll(item));

        // HEIGHT
        std::getline(line_stream, item, ' ');
        camera.height_ = std::stoll(item);
        // camera.SetHeight(std::stoll(item));

        // PARAMS
        // camera.Params().clear();
        camera.params_.clear();
        while (!line_stream.eof()) {
            std::getline(line_stream, item, ' ');
            // camera.Params().push_back(std::stold(item));
            camera.params_.push_back(std::stod(item));
        }

        cameras_.emplace(camera.camera_id_, camera);
    }
    return cameras_;
}

std::map<ColMapDataIO::image_t, ColMapDataIO::Image> ColMapDataIO::ReadImagesText(
    const std::string &path) {
    std::map<image_t, Image> images_;

    std::ifstream file(path);

    std::string line;
    std::string item;

    while (std::getline(file, line)) {
        StringTrim(&line);

        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::stringstream line_stream1(line);

        // ID
        std::getline(line_stream1, item, ' ');
        const image_t image_id = std::stoul(item);

        Image image;
        // image.SetImageId(image_id);
        image.image_id_ = image_id;

        // image.SetRegistered(true);
        // reg_image_ids_.push_back(image_id);

        // QVEC (qw, qx, qy, qz)
        std::getline(line_stream1, item, ' ');
        // image.Qvec(0) = std::stold(item);
        image.qvec_(0) = std::stod(item);

        std::getline(line_stream1, item, ' ');
        // image.Qvec(1) = std::stold(item);
        image.qvec_(1) = std::stod(item);

        std::getline(line_stream1, item, ' ');
        // image.Qvec(2) = std::stold(item);
        image.qvec_(2) = std::stod(item);

        std::getline(line_stream1, item, ' ');
        // image.Qvec(3) = std::stold(item);
        image.qvec_(3) = std::stod(item);

        // image.NormalizeQvec();
        image.qvec_.normalize();

        // TVEC
        std::getline(line_stream1, item, ' ');
        // image.Tvec(0) = std::stold(item);
        image.tvec_(0) = std::stod(item);

        std::getline(line_stream1, item, ' ');
        // image.Tvec(1) = std::stold(item);
        image.tvec_(1) = std::stod(item);

        std::getline(line_stream1, item, ' ');
        // image.Tvec(2) = std::stold(item);
        image.tvec_(2) = std::stod(item);

        // CAMERA_ID
        std::getline(line_stream1, item, ' ');
        // image.SetCameraId(std::stoul(item));
        image.camera_id_ = std::stoul(item);

        // NAME
        std::getline(line_stream1, item, ' ');
        // image.SetName(item);
        image.name_ = item;

        // POINTS2D
        if (!std::getline(file, line)) {
            break;
        }

        StringTrim(&line);
        std::stringstream line_stream2(line);

        std::vector<Eigen::Vector2d> points2D;
        std::vector<point3D_t> point3D_ids;

        if (!line.empty()) {
            while (!line_stream2.eof()) {
                Eigen::Vector2d point;

                std::getline(line_stream2, item, ' ');
                point.x() = std::stod(item);

                std::getline(line_stream2, item, ' ');
                point.y() = std::stod(item);

                points2D.push_back(point);

                std::getline(line_stream2, item, ' ');
                if (item == "-1") {
                    point3D_ids.push_back(kInvalidPoint3DId);
                } else {
                    point3D_ids.push_back(std::stoll(item));
                }
            }
        }

        // image.SetUp(Camera(image.CameraId()));
        image.SetPoints2D(points2D);

        for (point2D_t point2D_idx = 0;
             point2D_idx < static_cast<point2D_t>(image.points2D_.size()); ++point2D_idx) {
            if (point3D_ids[point2D_idx] != kInvalidPoint3DId) {
                image.SetPoint3DForPoint2D(point2D_idx, point3D_ids[point2D_idx]);
            }
        }

        images_.emplace(image.image_id_, image);
    }
    return images_;
}

std::map<ColMapDataIO::point3D_t, ColMapDataIO::Point3D> ColMapDataIO::ReadPoints3DText(
    const std::string &path) {
    std::map<point3D_t, Point3D> points3D_;

    std::ifstream file(path);

    std::string line;
    std::string item;

    while (std::getline(file, line)) {
        StringTrim(&line);

        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::stringstream line_stream(line);

        // ID
        std::getline(line_stream, item, ' ');
        const point3D_t point3D_id = std::stoll(item);

        // Make sure, that we can add new 3D points after reading 3D points
        // without overwriting existing 3D points.
        // num_added_points3D_ = std::max(num_added_points3D_, point3D_id);

        Point3D point3D;

        // XYZ
        std::getline(line_stream, item, ' ');
        // point3D.XYZ(0) = std::stold(item);
        point3D.xyz_(0) = std::stod(item);

        std::getline(line_stream, item, ' ');
        // point3D.XYZ(1) = std::stold(item);
        point3D.xyz_(1) = std::stod(item);

        std::getline(line_stream, item, ' ');
        // point3D.XYZ(2) = std::stold(item);
        point3D.xyz_(2) = std::stod(item);

        // Color
        std::getline(line_stream, item, ' ');
        // point3D.Color(0) = static_cast<uint8_t>(std::stoi(item));
        point3D.color_(0) = static_cast<uint8_t>(std::stoi(item));

        std::getline(line_stream, item, ' ');
        // point3D.Color(1) = static_cast<uint8_t>(std::stoi(item));
        point3D.color_(1) = static_cast<uint8_t>(std::stoi(item));

        std::getline(line_stream, item, ' ');
        // point3D.Color(2) = static_cast<uint8_t>(std::stoi(item));
        point3D.color_(2) = static_cast<uint8_t>(std::stoi(item));

        // ERROR
        std::getline(line_stream, item, ' ');
        // point3D.SetError(std::stold(item));
        point3D.error_ = std::stod(item);

        // TRACK
        while (!line_stream.eof()) {
            TrackElement track_el{};

            std::getline(line_stream, item, ' ');
            StringTrim(&item);
            if (item.empty()) {
                break;
            }

            track_el.image_id = std::stoul(item);

            std::getline(line_stream, item, ' ');
            track_el.point2D_idx = std::stoul(item);

            // point3D.Track().AddElement(track_el);
            point3D.track_.push_back(track_el);
        }

        // point3D.Track().Compress();

        points3D_.emplace(point3D_id, point3D);
    }
    return points3D_;
}

// -------------------------------
// camera, image, point3d, point3d
// -------------------------------
void ColMapDataIO::Image::SetPoints2D(const std::vector<Eigen::Vector2d> &points) {
    points2D_.resize(points.size());
    for (point2D_t point2D_idx = 0; point2D_idx < points.size(); ++point2D_idx) {
        points2D_[point2D_idx].xy_ = points[point2D_idx];
    }
}

void ColMapDataIO::Image::SetPoint3DForPoint2D(const ColMapDataIO::point2D_t point2D_idx,
                                               const ColMapDataIO::point3D_t point3D_id) {
    Point2D &point2D = points2D_.at(point2D_idx);
    point2D.point3D_id_ = point3D_id;
}

Eigen::Quaterniond ColMapDataIO::Image::QuatWorldToImg() const {
    return {qvec_(0), qvec_(1), qvec_(2), qvec_(3)};
}
}  // namespace ns_ikalibr
