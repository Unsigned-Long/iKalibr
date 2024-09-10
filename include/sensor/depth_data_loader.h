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

#ifndef IKALIBR_DEPTH_DATA_LOADER_H
#define IKALIBR_DEPTH_DATA_LOADER_H

#include "sensor/camera_data_loader.h"
#include "sensor/rgbd.h"

namespace {
bool IKALIBR_UNIQUE_NAME(_2_) = ns_ikalibr::_1_(__FILE__);
}

namespace ns_ikalibr {

class DepthDataLoader {
public:
    using Ptr = std::shared_ptr<DepthDataLoader>;

protected:
    CameraModelType _model;
    bool _isInverse;

public:
    explicit DepthDataLoader(CameraModelType model, bool isInverse);

    virtual DepthFrame ::Ptr UnpackFrame(const rosbag::MessageInstance &msgInstance) = 0;

    static DepthDataLoader::Ptr GetLoader(const std::string &modelStr, bool isInverse);

    [[nodiscard]] CameraModelType GetCameraModel() const;

    virtual ~DepthDataLoader() = default;

protected:
    template <class MsgType>
    void CheckMessage(typename MsgType::ConstPtr msg) {
        if (msg == nullptr) {
            throw std::runtime_error(
                "Wrong sensor model: '" + std::string(EnumCast::enumToString(GetCameraModel())) +
                "' for rgbd cameras! It's incompatible with the type of ros message to load in!");
        }
    }

    static void InverseMat(cv::Mat &floatMat);
};

class DepthSensorImageLoader : public DepthDataLoader {
public:
    using Ptr = std::shared_ptr<DepthSensorImageLoader>;

public:
    explicit DepthSensorImageLoader(CameraModelType model, bool isInverse);

    static DepthSensorImageLoader::Ptr Create(CameraModelType model, bool isInverse);

    DepthFrame::Ptr UnpackFrame(const rosbag::MessageInstance &msgInstance) override;
};

class DepthSensorImageCompLoader : public DepthDataLoader {
public:
    using Ptr = std::shared_ptr<DepthSensorImageCompLoader>;

public:
    explicit DepthSensorImageCompLoader(CameraModelType model, bool isInverse);

    static DepthSensorImageCompLoader::Ptr Create(CameraModelType model, bool isInverse);

    DepthFrame::Ptr UnpackFrame(const rosbag::MessageInstance &msgInstance) override;
};

}  // namespace ns_ikalibr

#endif  // IKALIBR_DEPTH_DATA_LOADER_H
