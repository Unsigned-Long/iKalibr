<div style="text-align: center;">
    <img src="../img/logo.svg" style="width: 100vw; height: auto;">
</div>

---

<h3 align="center">Tutorial: Prepare iKalibr Environment for Calibration</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>

---

<p align="left">
    <a><strong>Install Required Third Libraries »</strong></a>
</p> 

The following libraries need to be installed to support `iKalibr`. If you have already installed some of them, just skip corresponding installation. Some libraries may have some details that need to be paid attention to.

+ install `ROS1` (Ubuntu **20.04** is suggested, Ubuntu **18.04** (ros melodic) is also available), requirements: **ROS1** & **C++17** support.

  ```bash
  sudo apt install ros-noetic-desktop-full
  echo "source /opt/ros/noetic/setup.bash" >> ~/.bashrc
  source ~/.bashrc
  ```

+ install `Ceres`: see the `GitHub` Profile of **[Ceres](https://github.com/ceres-solver/ceres-solver.git)** library, clone it, compile it, and install it. Make sure that the version of `Ceres` contains the `Manifold` module and `Cuda` support. **`Ceres` version equals to 2.2.0 or higher.**

+ install `Sophus`: see the `GitHub` Profile of **[Sophus](https://github.com/strasdat/Sophus.git)** library, clone it, compile it, and install it. Set the cmake option `SOPHUS_USE_BASIC_LOGGING` as `ON` when compile `Sophus`!

+ install `magic-enum`: see the `GitHub` Profile of **[magic-enum](https://github.com/Neargye/magic_enum.git)** library, clone it, compile it, and install it.

+ install `Pangolin`: see the `GitHub` Profile of **[Pangolin](https://github.com/stevenlovegrove/Pangolin.git)** library, clone it, compile it, and install it.

+ install `cereal`, `yaml-cpp`, `spdlog`, and `colmap`. If possible, installing `colmap` from the [source](https://github.com/colmap/colmap.git) is recommend.

  ```bash
  sudo apt-get install libcereal-dev
  sudo apt-get install libyaml-cpp-dev
  # installing spdlog from source is a better way, to avoid fmt conflict
  sudo apt-get install libspdlog-dev
  # installing colmap from source is a better way, If the sensor suite you want to calibrate does not include optical cameras, you do not need to install colmap or glomap.
  sudo apt-get install colmap
  # if you want to use glomap to perform SfM, than clone it at (https://github.com/colmap/glomap.git), then build and install it. Generally speaking, glomap is recommanded for SfM (faster than colmap).
  ```

+ install ros packages:

  ```sh
  sudo apt-get install ros-noetic-cv-bridge
  sudo apt-get install ros-noetic-velodyne
  ```

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

  **Key point ** (you can't skip this part): the ros package name is `ikalibr` (the letter `k` is lowercase), rather than `iKalibr` (the letter `K` is uppercase).

  change directory to '`ikalibr`', and run '`build_thirdparty.sh`'.

  ```sh
  cd ikalibr
  chmod +x build_thirdparty.sh
  ./build_thirdparty.sh
  ```

  this would build sub module libraries, namely `tiny-viewer`, `ctraj`, `ufomap`, `veta`, and `opengv`.
  
  **Key point ** (you can't skip this part): the sub module `ufomap` would be checkout to `origin/devel_surfel` branch in `build_thirdparty.sh`. If errors about `ufomap` happen when compiling `iKalibr`, make sure the branch of `ufomap` is correct.

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

