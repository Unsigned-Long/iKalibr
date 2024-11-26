<div style="text-align: center;">
    <img src="../img/logo.svg" style="width: 100vw; height: auto;">
</div>

---

<h3 align="center">Tutorial: Prepare iKalibr Environment for Calibration</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>

---

**Attention, attention!** If you are using `Ubuntu 20.04`, please directly jump to [**A Tested Install Pipeline On Ubuntu 20.04**](#A Tested Install Pipeline On Ubuntu 20.04), where you will find the tested environment dependencies and versions of third-party libraries. If you are using a different version of `Ubuntu`, you can also prepare `iKalibr` environment in a similar way, but you may need to further adapt the library versions.

<p align="left">
    <a><strong>Install Required Third Libraries »</strong></a>
</p> 

The following libraries need to be installed to support `iKalibr`. If you have already installed some of them, just skip corresponding installation. Some libraries may have some details that need to be paid attention to.

+ install `ROS1` (Ubuntu **20.04** is suggested, Ubuntu **18.04** (ros melodic) is also available), requirements: **ROS1** & **C++17** support.

+ install `Ceres`: see the `GitHub` Profile of **[Ceres](https://github.com/ceres-solver/ceres-solver.git)** library, clone it, compile it, and install it. Make sure that the version of `Ceres` contains the `Manifold` module and `Cuda` support. **`Ceres` version equals to 2.2.0.**

+ install `Sophus`: see the `GitHub` Profile of **[Sophus](https://github.com/strasdat/Sophus.git)** library, clone it, compile it, and install it. Set the cmake option `SOPHUS_USE_BASIC_LOGGING` as `ON` when compile `Sophus`!

+ install `magic-enum`: see the `GitHub` Profile of **[magic-enum](https://github.com/Neargye/magic_enum.git)** library, clone it, compile it, and install it.

+ install `Pangolin`: see the `GitHub` Profile of **[Pangolin](https://github.com/stevenlovegrove/Pangolin.git)** library, clone it, compile it, and install it.

+ install `spdlog`: see the `GitHub` Profile of **[spdlog](https://github.com/gabime/spdlog.git)** library, clone it, compile it, and install it.

+ install `cereal`, `yaml-cpp`, and `colmap`. If possible, installing `colmap` from the [source](https://github.com/colmap/colmap.git) is recommend.

+ install required ros packages.


**Key point** (you can't skip this part): 

+ both `spdlog` and `Sophus` would involve the `fmt` library, and would lead to conflict if the versions of involved `fmt` in `spdlog` and `Sophus` are different. To solve this, we recommend that you set the cmake option `SOPHUS_USE_BASIC_LOGGING` as `ON` when compile `Sophus`, this would avoid to involve `fmt` logger dependency in `Sophus`. If conflict about `fmt` still exists, it is recommended that you try to compile the `spdlog` library from [source code](https://github.com/gabime/spdlog.git) (directly compile and install), so that `spdlog` will use the internal `fmt` library.
+ Third-party dependency libraries that do not specify version numbers use the system default version. For example, the default `pcl` version of `Ubuntu 20.04` is `1.10` (`pcl` can be installed together with `ROS`).

<p align="left">
    <a><strong>Clone iKalibr and Its Modules »</strong></a>
</p> 

+ create a ros workspace if needed and clone `iKalibr` to `src` directory as `ikalibr`:

  ```sh
  mkdir -p ~/iKalibr/src
  cd ~/iKalibr/src
  
  git clone --recursive https://github.com/Unsigned-Long/iKalibr.git ikalibr
  ```

  **Key point** (you can't skip this part): the ros package name is `ikalibr` (the letter `k` is lowercase), rather than `iKalibr` (the letter `K` is uppercase).

  change directory to '`ikalibr`', and run '`build_thirdparty.sh`'.

  ```sh
  cd ikalibr
  chmod +x build_thirdparty.sh
  ./build_thirdparty.sh
  ```

  this would build sub module libraries, namely `tiny-viewer`, `ctraj`, `ufomap`, `veta`, and `opengv`.
  
  **Key point** (you can't skip this part): the sub module `ufomap` would be checkout to `origin/devel_surfel` branch in `build_thirdparty.sh`. If errors about `ufomap` happen when compiling `iKalibr`, make sure the branch of `ufomap` is correct.

<p align="left">
    <a><strong>Compile iKalibr »</strong></a>
</p> 

+ change directory to catkin workspace, and generate the ros self-defined messages in `ikalibr`:

  ```sh
  cd ../..
  catkin_make ikalibr_generate_messages
  ```

+ compile `ikalibr` package:

  ```sh
  catkin_make -j8 -DUSE_CMAKE_UNITY_BUILD=ON
  ```

Congratulations :clap:, if everything goes well and no error happened. At the end, you would obtain several binary programs, such as `ikalibr_prog`, `ikalibr_imu_intri_calib`, etc. Each program is exactly an executable ros node, and can be launched by `rosrun` or provided `roslaunch` (recommend).



# A Tested Install Pipeline On Ubuntu 20.04

Considering the common issues reported by most people regarding the configuration and compilation of the `iKalibr` environment, I will provide here a step-by-step guide for setting up the `iKalibr` environment from scratch on an `Ubuntu 20.04` system. All steps have been tested and successfully compiled `iKalibr` on Ubuntu 20.04.
Please note that `COLMAP` and `GLOMAP` are not required when compiling `iKalibr`, as they are only used to provide `SfM` results within `iKalibr` and are decoupled from it. If you need to calibrate the camera, then you should compile them (for mapping-based visual-inertial calibration in `iKalibr`); otherwise, there is no need to compile or install them.

**Attention:** This install pipeline has already been, and only been tested on `Ubuntu 20.04`.

```sh
# ros 1 (noetic), add install source and key first, see https://wiki.ros.org/noetic/Installation/Ubuntu
sudo apt install ros-noetic-desktop-full
echo "source /opt/ros/noetic/setup.bash" >> ~/.bashrc
source ~/.bashrc

# cmake (3.27)
wget https://github.com/Kitware/CMake/releases/download/v3.27.0/cmake-3.27.0.tar.gz && \
    tar -zxvf cmake-3.27.0.tar.gz && \
    cd cmake-3.27.0 && \
    ./bootstrap && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf cmake-3.27.0 cmake-3.27.0.tar.gz

# eigen (use default one in ubuntu 20.04)

# ceres (2.2.0)
git clone --branch 2.2.0 --single-branch --recurse-submodules https://github.com/ceres-solver/ceres-solver
# then cmake, make, and install

# Sophus (1.22.10)
git clone --branch 1.22.10 --single-branch  https://github.com/strasdat/Sophus.git
cmake .. -DSOPHUS_USE_BASIC_LOGGING=ON
# then make and install

# magic_enum (v0.9.6)
git clone --branch v0.9.6 --single-branch https://github.com/Neargye/magic_enum.git
# then cmake, make, and install

# pangolin (v0.8)
git clone --branch v0.8 --single-branch --recursive https://github.com/stevenlovegrove/Pangolin.git
# then cmake, make, and install

# spdlog (the newest, use internal fmt)
git clone https://github.com/gabime/spdlog.git
# then cmake, make, and install

# install cereal and yaml-cpp
sudo apt-get install libcereal-dev
sudo apt-get install libyaml-cpp-dev

# ros packages
sudo apt-get install ros-noetic-cv-bridge
sudo apt-get install ros-noetic-velodyne

# install colmap and glomap from source, if you want to calibrate your optical cameras using mapping-based visual-inertial calibration (for mapping-free visual-inertial calibration, colmap and glomap are not required)
https://github.com/colmap/glomap.git
https://github.com/colmap/colmap.git
# then cmake, make, and install

# clone iKalibr
mkdir -p ~/iKalibr/src
cd ~/iKalibr/src
git clone --recursive https://github.com/Unsigned-Long/iKalibr.git ikalibr

# build thirdparty
cd ikalibr
chmod +x build_thirdparty.sh
./build_thirdparty.sh

# build iKalibr
cd ../..
catkin_make ikalibr_generate_messages
catkin_make -j8 -DUSE_CMAKE_UNITY_BUILD=ON

# test
source ./devel/setup.bash
roslaunch ikalibr ikalibr-learn.launch
```
