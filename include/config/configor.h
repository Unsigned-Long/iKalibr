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

#ifndef IKALIBR_CONFIGOR_H
#define IKALIBR_CONFIGOR_H

/**
 * @attention Description about "veta/type_def.hpp"
 'EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST' is a macro used to handle Eigen objects within STL containers.
 In C++, we often use standard library containers (e.g., std::vector) to store data.
However, when we attempt to store Eigen objects (e.g., Eigen::Vector2d) in these containers, we encounter some issues.
 When using Eigen objects in STL containers, special handling is required to ensure proper memory alignment and initialization.
'EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST' exists to address such issues.
 The macro allows us to use Eigen objects in STL containers without manually dealing with memory alignment and initialization.
With this macro, we can use std::vector<Eigen::Vector2d> just like any other vector, without worrying about the underlying implementation details.
 This specialization needs to be defined before using code snippets like std::vector<Vector2d>.
The benefit is that you don’t need to declare std::vector with Eigen::aligned_allocator everywhere in your code. However, the downside is that the specialization must be defined before using code snippets like std::vector<Vector2d.
 */
#include "veta/type_def.hpp"
#include "memory"
#include "cereal/archives/json.hpp"
#include "cereal/cereal.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/set.hpp"
#include "util/utils.hpp"
#include "util/enum_cast.hpp"
#include "sensor/imu_intrinsic.hpp"
#include "thread"
#include "ros/package.h"

_3_

namespace ns_ikalibr {

    struct Configor {
    public:
        using Ptr = std::shared_ptr<Configor>;

    public:

        static struct DataStream {
            struct IMUConfig {
            public:
                std::string Type;
                std::string Intrinsics;
                double AcceWeight;
                double GyroWeight;

                IMUConfig() : Type(), Intrinsics(), AcceWeight(), GyroWeight() {};

            public:
                template<class Archive>
                void serialize(Archive &ar) {
                    ar(CEREAL_NVP(Type), CEREAL_NVP(Intrinsics), CEREAL_NVP(AcceWeight), CEREAL_NVP(GyroWeight));
                }
            };

            struct RadarConfig {
            public:
                std::string Type;
                double Weight;

                RadarConfig() : Type(), Weight() {};

            public:
                template<class Archive>
                void serialize(Archive &ar) {
                    ar(CEREAL_NVP(Type), CEREAL_NVP(Weight));
                }
            };

            struct LiDARConfig {
            public:
                std::string Type;
                double Weight;

                LiDARConfig() : Type(), Weight() {};

            public:
                template<class Archive>
                void serialize(Archive &ar) {
                    ar(CEREAL_NVP(Type), CEREAL_NVP(Weight));
                }
            };

            struct CameraConfig {
            public:
                std::string Type;
                std::string Intrinsics;
                double Weight;
                int TrackLengthMin;

                CameraConfig() : Type(), Intrinsics(), Weight(), TrackLengthMin() {};

            public:
                template<class Archive>
                void serialize(Archive &ar) {
                    ar(CEREAL_NVP(Type), CEREAL_NVP(Intrinsics), CEREAL_NVP(Weight), CEREAL_NVP(TrackLengthMin));
                }
            };

            static std::optional<std::string> CreateImageStoreFolder(const std::string &camTopic);

            static std::optional<std::string> CreateSfMWorkspace(const std::string &camTopic);

            static std::string GetImageStoreInfoFile(const std::string &camTopic);

            static std::map<std::string, IMUConfig> IMUTopics;
            static std::map<std::string, RadarConfig> RadarTopics;
            static std::map<std::string, LiDARConfig> LiDARTopics;
            static std::map<std::string, CameraConfig> CameraTopics;

            static std::string ReferIMU;

            static std::string BagPath;
            static double BeginTime;
            static double Duration;

            static std::string OutputPath;
            const static std::string PkgPath;
            const static std::string DebugPath;

        public:
            template<class Archive>
            void serialize(Archive &ar) {
                ar(
                        CEREAL_NVP(IMUTopics), CEREAL_NVP(RadarTopics), CEREAL_NVP(LiDARTopics),
                        CEREAL_NVP(CameraTopics), CEREAL_NVP(ReferIMU), CEREAL_NVP(BagPath),
                        CEREAL_NVP(BeginTime), CEREAL_NVP(Duration), CEREAL_NVP(OutputPath)
                );
            }
        } dataStream;

        static struct Prior {
            static double GravityNorm;
            static constexpr int SplineOrder = 4;
            static bool OptTemporalParams;
            static double TimeOffsetPadding;
            static double ReadoutTimePadding;

            static struct KnotTimeDist {
                static double SO3Spline;
                static double ScaleSpline;

            public:
                template<class Archive>
                void serialize(Archive &ar) {
                    ar(CEREAL_NVP(SO3Spline), CEREAL_NVP(ScaleSpline));
                }
            } knotTimeDist;

            static struct NDTLiDAROdometer {
                static double Resolution;
                static double KeyFrameDownSample;

            public:
                template<class Archive>
                void serialize(Archive &ar) {
                    ar(CEREAL_NVP(Resolution), CEREAL_NVP(KeyFrameDownSample));
                }
            } ndtLiDAROdometer;

            static struct LiDARDataAssociate {
                static double MapDownSample;
                static double PointToSurfelMax;
                static double PlanarityMin;

                static constexpr double QueryDepthMin = 1;
                static constexpr double QueryDepthMax = 2;
                static constexpr double SurfelPointMin = 100;
                //   0,   1,   2,   3,   4, ...
                // 0.1, 0.2, 0.4, 0.8, 1.6, ...
                static constexpr double MapResolution = 0.1;
                static constexpr std::uint8_t MapDepthLevels = 16;

                static constexpr double PointToSurfelCountInScan = 200;

            public:
                template<class Archive>
                void serialize(Archive &ar) {
                    ar(CEREAL_NVP(MapDownSample), CEREAL_NVP(PointToSurfelMax), CEREAL_NVP(PlanarityMin));
                }
            } lidarDataAssociate;

            static double CauchyLossForRadarFactor;
            static double CauchyLossForLiDARFactor;
            static double CauchyLossForCameraFactor;
        public:
            template<class Archive>
            void serialize(Archive &ar) {
                ar(CEREAL_NVP(GravityNorm), CEREAL_NVP(OptTemporalParams), CEREAL_NVP(TimeOffsetPadding),
                   CEREAL_NVP(ReadoutTimePadding), cereal::make_nvp("KnotTimeDist", knotTimeDist),
                   cereal::make_nvp("NDTLiDAROdometer", ndtLiDAROdometer),
                   cereal::make_nvp("LiDARDataAssociate", lidarDataAssociate),
                   CEREAL_NVP(CauchyLossForRadarFactor), CEREAL_NVP(CauchyLossForLiDARFactor),
                   CEREAL_NVP(CauchyLossForCameraFactor));
            }
        } prior;

        static struct Preference {
            static bool UseCudaInSolving;
            static bool OutputParamInEachIter;
            static bool OutputBSplines;
            static bool OutputMaps;
            static bool OutputHessianMat;
            static bool OutputVisualLiDARCovisibility;
            static bool OutputVisualKinematics;
            static bool OutputColorizedMap;
            static bool OutputAlignedInertialMes;
            static std::string OutputDataFormat;
            const static std::map<CerealArchiveType::Enum, std::string> FileExtension;
            static int ThreadsToUse;

            const static std::string SO3_SPLINE, SCALE_SPLINE;

            // in visualizator
            static double SplineScaleInViewer;
            static double CoordSScaleInViewer;

            static int AvailableThreads();

            static CerealArchiveType::Enum DataIOFormat();

        public:
            template<class Archive>
            void serialize(Archive &ar) {
                ar(
                        CEREAL_NVP(UseCudaInSolving),
                        CEREAL_NVP(OutputParamInEachIter),
                        CEREAL_NVP(OutputBSplines),
                        CEREAL_NVP(OutputMaps),
                        CEREAL_NVP(OutputHessianMat),
                        CEREAL_NVP(OutputVisualLiDARCovisibility),
                        CEREAL_NVP(OutputVisualKinematics),
                        CEREAL_NVP(OutputColorizedMap),
                        CEREAL_NVP(OutputAlignedInertialMes),
                        CEREAL_NVP(OutputDataFormat),
                        CEREAL_NVP(ThreadsToUse),
                        CEREAL_NVP(SplineScaleInViewer),
                        CEREAL_NVP(CoordSScaleInViewer)
                );
            }
        } preference;

    public:
        Configor();

        static Ptr Create();

        // load configure information from file
        template<class CerealArchiveType=CerealArchiveType::YAML>
        static bool LoadConfigure(const std::string &filename) {
            std::ifstream file(filename);
            auto archive = GetInputArchive<CerealArchiveType>(file);
            auto configor = Configor::Create();
            (*archive)(cereal::make_nvp("Configor", *configor));
            configor->CheckConfigure();
            return true;
        }

        // save configure information to file
        template<class CerealArchiveType=CerealArchiveType::YAML>
        bool SaveConfigure(const std::string &filename) {
            std::ofstream file(filename);
            auto archive = GetOutputArchive<CerealArchiveType>(file);
            (*archive)(cereal::make_nvp("Configor", *this));
            return true;
        }

        // load configure information from file
        static bool LoadConfigure(const std::string &filename, CerealArchiveType::Enum archiveType);

        // save configure information to file
        bool SaveConfigure(const std::string &filename, CerealArchiveType::Enum archiveType);

        // print the main fields
        static void PrintMainFields();

        static std::string GetFormatExtension();

        [[nodiscard]] static bool IsLiDARIntegrated();

        [[nodiscard]] static bool IsCameraIntegrated();

        [[nodiscard]] static bool IsRadarIntegrated();

    protected:
        // check the input configure
        static void CheckConfigure();

    public:
        template<class Archive>
        void serialize(Archive &ar) {
            ar(cereal::make_nvp("DataStream", dataStream), cereal::make_nvp("Prior", prior),
               cereal::make_nvp("Preference", preference));
        }
    };
}

#endif //IKALIBR_CONFIGOR_H
