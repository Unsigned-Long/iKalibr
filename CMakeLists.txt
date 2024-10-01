# iKalibr: Unified Targetless Spatiotemporal Calibration Framework
# Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China
# https://github.com/Unsigned-Long/iKalibr.git
#
# Author: Shuolong Chen (shlchen@whu.edu.cn)
# GitHub: https://github.com/Unsigned-Long
#  ORCID: 0000-0002-5283-9057
#
# Purpose: See .h/.hpp file.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * The names of its contributors can not be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.0.2)
project(ikalibr)

set(USE_CMAKE_UNITY_BUILD ON CACHE STRINGS "whether use 'CMAKE_UNITY_BUILD' in building")

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16" AND USE_CMAKE_UNITY_BUILD)
    set(CMAKE_UNITY_BUILD ON)
    message(STATUS "use 'CMAKE_UNITY_BUILD' in building!")
else ()
    message(STATUS "do not use 'CMAKE_UNITY_BUILD' in building!")
endif ()

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++17)
set(CMAKE_BUILD_TYPE "RELEASE")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall   -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall   -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

if (NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif ()

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
        roscpp
        rospy
        std_msgs
        geometry_msgs
        sensor_msgs
        rosbag
        cv_bridge
        message_generation
        velodyne_msgs
        velodyne_pointcloud
)

## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)


## Uncomment this if the package has a setup.py. This macro ensures
## modules and global scripts declared therein get installed
## See http://ros.org/doc/api/catkin/html/user_guide/setup_dot_py.html
# catkin_python_setup()

################################################
## Declare ROS messages, services and actions ##
################################################

## To declare and build messages, services or actions from within this
## package, follow these steps:
## * Let MSG_DEP_SET be the set of packages whose message types you use in
##   your messages/services/actions (e.g. std_msgs, actionlib_msgs, ...).
## * In the file package.xml:
##   * add a build_depend tag for "message_generation"
##   * add a build_depend and a exec_depend tag for each package in MSG_DEP_SET
##   * If MSG_DEP_SET isn't empty the following dependency has been pulled in
##     but can be declared for certainty nonetheless:
##     * add a exec_depend tag for "message_runtime"
## * In this file (CMakeLists.txt):
##   * add "message_generation" and every package in MSG_DEP_SET to
##     find_package(catkin REQUIRED COMPONENTS ...)
##   * add "message_runtime" and every package in MSG_DEP_SET to
##     catkin_package(CATKIN_DEPENDS ...)
##   * uncomment the add_*_files sections below as needed
##     and list every .msg/.srv/.action file to be processed
##   * uncomment the generate_messages entry below
##   * add every package in MSG_DEP_SET to generate_messages(DEPENDENCIES ...)

## Generate messages in the 'msg' folder
add_message_files(
        FILES

        AinsteinRadarTarget.msg
        AinsteinRadarTargetArray.msg

        SbgImuStatus.msg
        SbgImuData.msg

        LivoxCustomPoint.msg
        LivoxCustomMsg.msg

        AWR1843RadarScan.msg
        AWR1843RadarScanCustom.msg
)

## Generate services in the 'srv' folder
# add_service_files(
#   FILES
#   Service1.srv
#   Service2.srv
# )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
generate_messages(
        DEPENDENCIES
        std_msgs
        geometry_msgs
)

################################################
## Declare ROS dynamic reconfigure parameters ##
################################################

## To declare and build dynamic reconfigure parameters within this
## package, follow these steps:
## * In the file package.xml:
##   * add a build_depend and a exec_depend tag for "dynamic_reconfigure"
## * In this file (CMakeLists.txt):
##   * add "dynamic_reconfigure" to
##     find_package(catkin REQUIRED COMPONENTS ...)
##   * uncomment the "generate_dynamic_reconfigure_options" section below
##     and list every .cfg file to be processed

## Generate dynamic reconfigure parameters in the 'cfg' folder
# generate_dynamic_reconfigure_options(
#   cfg/DynReconf1.cfg
#   cfg/DynReconf2.cfg
# )

###########
## Build ##
###########

set(tiny-viewer_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ctraj/thirdparty/tiny-viewer-install/lib/cmake/tiny-viewer)
find_package(tiny-viewer)
set(ctraj_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ctraj-install/lib/cmake/ctraj)
find_package(ctraj)
set(ufomap_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ufomap-install/lib/cmake/ufomap)
find_package(ufomap)
set(ufomap_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ufomap-install/include)
include_directories(${ufomap_INCLUDE_DIR})
set(veta_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/veta-install/lib/cmake/veta)
find_package(veta)
set(opengv_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/opengv-install/lib/cmake/opengv-1.0)
find_package(opengv)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

find_package(magic_enum)
# magic_enum works in range of [MAGIC_ENUM_RANGE_MIN, MAGIC_ENUM_RANGE_MAX]
add_definitions(-DMAGIC_ENUM_RANGE_MIN=0)
# 1024 = 2^10
add_definitions(-DMAGIC_ENUM_RANGE_MAX=1023)
find_package(yaml-cpp)
find_package(PCL REQUIRED)

#define PCL_NO_PRECOMPILE [it's very important for Custom pcl Point type]
add_definitions(-DPCL_NO_PRECOMPILE)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/calib CALIB_SRC_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/solver SOLVER_SRC_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/core CORE_SRC_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/sensor SENSOR_SRC_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/config CONFIG_SRC_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/viewer VIEWER_SRC_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/util UTIL_SRC_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src/factor FACTOR_SRC_FILES)


## Specify additional locations of header files
## Your package locations should be listed before other locations
###################
# libikalibr_util #
###################
add_library(
        ${PROJECT_NAME}_util SHARED
        # source files
        ${UTIL_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_util PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_util PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty
        ctraj
        veta
)

#####################
# libikalibr_sensor #
#####################
add_library(
        ${PROJECT_NAME}_sensor SHARED
        # source files
        ${SENSOR_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_sensor PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_sensor PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty
        ctraj
        veta
)

#####################
# libikalibr_config #
#####################
add_library(
        ${PROJECT_NAME}_config SHARED
        # source files
        ${CONFIG_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_config PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_config PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty
        ctraj
        veta
)

#####################
# libikalibr_factor #
#####################
add_library(
        ${PROJECT_NAME}_factor SHARED
        # source files
        ${FACTOR_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_factor PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_factor PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty
        ctraj
        veta
        UFO::Map
        opengv
)
######################
# precompile headers #
######################
#file(GLOB FACTOR_HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/factor/*.h"
#        "${CMAKE_CURRENT_SOURCE_DIR}/include/factor/*.hpp")
#target_precompile_headers(
#        ${PROJECT_NAME}_factor PUBLIC
#        ${FACTOR_HEADER_FILES}
#)

#####################
# libikalibr_viewer #
#####################
add_library(
        ${PROJECT_NAME}_viewer SHARED
        # source files
        ${VIEWER_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_viewer PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_viewer PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty
        ctraj
        veta
        UFO::Map
)
#target_precompile_headers(
#        ${PROJECT_NAME}_viewer REUSE_FROM
#        ${PROJECT_NAME}_factor
#)
###################
# libikalibr_core #
###################
add_library(
        ${PROJECT_NAME}_core SHARED
        # source files
        ${CORE_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_core PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_core PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty
        ctraj
        veta
        opengv
        UFO::Map
)
#target_precompile_headers(
#        ${PROJECT_NAME}_core REUSE_FROM
#        ${PROJECT_NAME}_factor
#)
####################
# libikalibr_calib #
####################
add_library(
        ${PROJECT_NAME}_calib SHARED
        # source files
        ${CALIB_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_calib PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_calib PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty
        ctraj
        veta
        opengv
)
#target_precompile_headers(
#        ${PROJECT_NAME}_calib REUSE_FROM
#        ${PROJECT_NAME}_factor
#)
#####################
# libikalibr_solver #
#####################
add_library(
        ${PROJECT_NAME}_solver SHARED
        # source files
        ${SOLVER_SRC_FILES}
)

target_include_directories(
        ${PROJECT_NAME}_solver PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
        ${PROJECT_NAME}_solver PUBLIC
        ${catkin_LIBRARIES}
        # thirdparty

        veta
        ctraj
        opengv
)
#target_precompile_headers(
#        ${PROJECT_NAME}_solver REUSE_FROM
#        ${PROJECT_NAME}_factor
#)
###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if your package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
        # INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
        #        LIBRARIES ${PROJECT_NAME}
        CATKIN_DEPENDS roscpp rospy std_msgs geometry_msgs message_runtime
        #  DEPENDS system_lib
)

## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/${PROJECT_NAME}/ikalibr.cpp
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide
# add_executable(${PROJECT_NAME}_node src/ikalibr_node.cpp)

add_executable(
        ${PROJECT_NAME}_prog
        exe/solver/main.cpp
)
add_executable(
        ${PROJECT_NAME}_learn
        exe/nofree/learn.cpp
        src/nofree/logo_svg.cpp
        src/nofree/data_collect_demo.cpp
)
add_executable(
        ${PROJECT_NAME}_lidar_map_viewer
        exe/tool/lidar_map_viewer.cpp
)
add_executable(
        ${PROJECT_NAME}_data_format_transformer
        exe/tool/data_format_transformer.cpp
)
add_executable(
        ${PROJECT_NAME}_bag_merge
        exe/tool/bag_merge.cpp
        src/nofree/bag_merge.cpp
)
add_executable(
        ${PROJECT_NAME}_bag_topic_downsample
        exe/tool/bag_topic_downsample.cpp
)
add_executable(
        ${PROJECT_NAME}_imgs_to_bag
        exe/tool/imgs_to_bag.cpp
)
add_executable(
        ${PROJECT_NAME}_imu_intri_calib
        exe/tool/imu_intri_calib.cpp
        src/nofree/imu_intri_calib.cpp
)
add_executable(
        ${PROJECT_NAME}_raw_inertial_to_bag
        exe/tool/raw_inertial_to_bag.cpp
)

## Rename C++ executable without prefix
## The above recommended prefix causes long target names, the following renames the
## target back to the shorter version for ease of user use
## e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
# target_link_libraries(${PROJECT_NAME}_node
#   ${catkin_LIBRARIES}
# )
###################
# libikalibr_prog #
###################
target_include_directories(
        ${PROJECT_NAME}_prog PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_prog PRIVATE

        # the dependent library is placed after the library that depends on it.
        ${PROJECT_NAME}_solver
        ${PROJECT_NAME}_calib
        ${PROJECT_NAME}_factor
        ${PROJECT_NAME}_core
        ${PROJECT_NAME}_viewer
        ${PROJECT_NAME}_sensor
        ${PROJECT_NAME}_config
        ${PROJECT_NAME}_util

        # thirdparty
        ${YAML_CPP_LIBRARIES}
)
####################
# libikalibr_learn #
####################
target_include_directories(
        ${PROJECT_NAME}_learn PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_learn PRIVATE

        # iKalibr libraries
        ${PROJECT_NAME}_util

        # thirdparty
        ${YAML_CPP_LIBRARIES}
)
###############################
# libikalibr_lidar_map_viewer #
###############################
target_include_directories(
        ${PROJECT_NAME}_lidar_map_viewer PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_lidar_map_viewer PRIVATE

        # thirdparty
        ${PROJECT_NAME}_util
)
######################################
# libikalibr_data_format_transformer #
######################################
target_include_directories(
        ${PROJECT_NAME}_data_format_transformer PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_data_format_transformer PRIVATE

        # the dependent library is placed after the library that depends on it.
        ${PROJECT_NAME}_calib
        ${PROJECT_NAME}_factor
        ${PROJECT_NAME}_core
        ${PROJECT_NAME}_viewer
        ${PROJECT_NAME}_sensor
        ${PROJECT_NAME}_config
        ${PROJECT_NAME}_util

        # thirdparty
        ${YAML_CPP_LIBRARIES}
)
########################
# libikalibr_bag_merge #
########################
target_include_directories(
        ${PROJECT_NAME}_bag_merge PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_bag_merge PRIVATE

        # thirdparty
        ${PROJECT_NAME}_util
        ${YAML_CPP_LIBRARIES}
)
###################################
# libikalibr_bag_topic_downsample #
###################################
target_include_directories(
        ${PROJECT_NAME}_bag_topic_downsample PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_bag_topic_downsample PRIVATE

        # thirdparty
        ${PROJECT_NAME}_util
        ${YAML_CPP_LIBRARIES}
)
##########################
# libikalibr_imgs_to_bag #
##########################
target_include_directories(
        ${PROJECT_NAME}_imgs_to_bag PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_imgs_to_bag PRIVATE

        # thirdparty
        ${PROJECT_NAME}_util
)
##############################
# libikalibr_imu_intri_calib #
##############################
target_include_directories(
        ${PROJECT_NAME}_imu_intri_calib PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_imu_intri_calib PRIVATE

        # the dependent library is placed after the library that depends on it.
        ${PROJECT_NAME}_calib
        ${PROJECT_NAME}_factor
        ${PROJECT_NAME}_core
        ${PROJECT_NAME}_viewer
        ${PROJECT_NAME}_sensor
        ${PROJECT_NAME}_config
        ${PROJECT_NAME}_util

        # thirdparty
        ${YAML_CPP_LIBRARIES}
)
##################################
# libikalibr_raw_inertial_to_bag #
##################################
target_include_directories(
        ${PROJECT_NAME}_raw_inertial_to_bag PUBLIC
        # include
        ${catkin_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
        ${PROJECT_NAME}_raw_inertial_to_bag PRIVATE

        # thirdparty
        ${PROJECT_NAME}_util
)

#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# catkin_install_python(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables for installation
## See http://docs.ros.org/melodic/api/catkin/html/howto/format1/building_executables.html
# install(TARGETS ${PROJECT_NAME}_node
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark libraries for installation
## See http://docs.ros.org/melodic/api/catkin/html/howto/format1/building_libraries.html
# install(TARGETS ${PROJECT_NAME}
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_ikalibr.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
