<div align=center><img src="../img/logo.svg" width =100%></div>

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

+ install `Sophus`: see the `GitHub` Profile of **[Sophus](https://github.com/strasdat/Sophus.git)** library, clone it, compile it, and install it.

+ install `magic-enum`: see the `GitHub` Profile of **[magic-enum](https://github.com/Neargye/magic_enum.git)** library, clone it, compile it, and install it.

+ install `Pangolin`: see the `GitHub` Profile of **[Pangolin](https://github.com/stevenlovegrove/Pangolin.git)** library, clone it, compile it, and install it.

+ install `cereal`, `yaml-cpp`, `spdlog`, and `colmap`. If possible, installing `colmap` from the [source](https://github.com/colmap/colmap.git) is recommend.

  ```bash
  sudo apt-get install libcereal-dev
  sudo apt-get install libyaml-cpp-dev
  sudo apt-get install libspdlog-dev
  # installing colmap from source is a better way
  sudo apt-get install colmap
  ```

**Key point**: both `spdlog` and `Sophus` would involve the `fmt` library, and would lead to conflict if the versions of involved `fmt` in `spdlog` and `Sophus` are different. To solve this, we recommend that you set the cmake option `SOPHUS_USE_BASIC_LOGGING` as `ON` when compile `Sophus`, this would avoid to involve `fmt` logger dependency in `Sophus`.

<p align="left">
    <a><strong>Clone iKalibr and Its Modules »</strong></a>
</p> 

+ create a ros workspace if needed and clone `iKalibr` to `src` directory as `ikalibr`:

  ```sh
  mkdir -p ~/iKalibr/src
  cd ~/iKalibr/src
  
  git clone --recursive https://github.com/Unsigned-Long/iKalibr.git ikalibr
  ```

  **Key point**: the ros package name is `ikalibr` (the letter `k` is lowercase), rather than `iKalibr` (the letter `K` is uppercase).

  change directory to '`ikalibr`', and run '`build_thirdparty.sh`'.

  ```sh
  cd ikalibr
  chmod +x build_thirdparty.sh
  ./build_thirdparty.sh
  ```

  this would build sub module libraries, namely `tiny-viewer`, `ctraj`, `ufomap`, `veta`, and `opengv`.

<p align="left">
    <a><strong>Compile iKalibr »</strong></a>
</p> 

+ generate the ros self-defined messages in `ikalibr`:

  ```sh
  catkin_make ikalibr_generate_messages
  ```

+ compile `ikalibr` package:

  ```sh
  catkin_make -j8
  ```

Congratulations :clap: , if everything goes well and no error happened. At the end, you would obtain several binary programs, such as `ikalibr_prog`, `ikalibr_imu_intri_calib`, etc. Each program is exactly an executable ros node, and can be launched by `rosrun` or provided `roslaunch` (recommend).

