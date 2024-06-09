<div align=center><img src="../img/logo.svg" width =100%></div>

---

<h3 align="center">Additional Helpful Tools in iKalibr</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>

---

<p align="left">
    <a><strong>IMU Intrinsic Calibration »</strong></a>
</p> 

If you want to perform IMU-only multi-IMU spatiotemporal calibration in `iKalibr`, multiple IMUs are required to be pre-calibrated (intrinsics). This can be conducted based on `ikalibr-imu-intri-calib` program. Specifically, you need to modify the [config-imu-intri-calib.yaml](../../config/tool/config-imu-intri-calib.yaml) file to adapt. Here are the specific steps to implement it:

+ collect inertial measurements from the IMU, which undergoes several **stationary** poses (generally six symmetric poses).

+ based on the collected rosbag, find the time segment corresponding to stationary moments. [PlotJuggler](https://github.com/facontidavide/PlotJuggler.git) is suggested  for static time piece checking, and you can install it by:

  ```sh
  sudo apt-get install ros-noetic-plotjuggler
  sudo apt-get install ros-noetic-plotjuggler-ros
  # launch plotjuggler and drag your rosbag into it
  rosrun plotjuggler plotjuggler
  ```

  Afterwards, the specific time period can be determined based on the distribution of inertial data.

+ write the corresponding specific time period information into the configuration file [config-imu-intri-calib.yaml](../../config/tool/config-imu-intri-calib.yaml). Then, launch [ikalibr-imu-intri-calib](../../launch/tool/ikalibr-imu-intri-calib.launch) by:

  ```sh
  roslaunch ikalibr ikalibr-imu-intri-calib.launch
  ```

  The final output result file can be directly used for the IMU internal parameter setting in `ikalibr_prog`. 

**Key point**: If you calibrate sensor suite that has additional sensors beside the IMU, such as radar-inertial suite, then **you do not need to calibrate the IMU intrinsic** parameters separately, because they will be optimized during the spatiotemporal calibration.

<p align="left">
    <a><strong>Data Format Transformer »</strong></a>
</p> 

`iKalibr ` supports four data formats (this benefits from the [cereal](https://github.com/USCiLab/cereal.git) serialization library), namely `yaml`, `xml`, `json`, and `binary`. If you specify one of the data formats during the solution, you can convert between these four data formats after the solving. After the operation, all result files in the workspace will be converted (generating files with the same name but different extensions). 

Specifically, you can achieve this by [ikalibr-data-format-transformer](../../launch/tool/ikalibr-data-format-transformer.launch). For specific configuration, see the instructions inside. Then, you can run:

```sh
roslaunch ikalibr ikalibr-data-format-transformer.launch
```

<p align="left">
    <a><strong>Rosbag Merge and Split »</strong></a>
</p> 

If you collect data from each sensor separately into a separate rosbag, you can merge them into one rosbag using [ikalibr-bag-merge](../../launch/tool/ikalibr-bag-merge.launch) before solving. Specifically, you should first configure the [config-bag-merge.yaml](../../config/tool/config-bag-merge.yaml) file, and then run:

```sh
roslaunch ikalibr ikalibr-bag-merge.launch
```

<p align="left">
    <a><strong>Images To Rosbag »</strong></a>
</p> 

If you store the image directly as a single file on disk when you collect sensor data, you may need to make them as a rosbag. In this case, you can use [ikalibr-imgs-to-bag](../../launch/tool/ikalibr-imgs-to-bag.launch). After you have configured the contents in the launch file, run:

```sh
roslaunch ikalibr ikalibr-imgs-to-bag.launch
```


<p align="left">
    <a><strong>Raw Inertial Measurements To Rosbag »</strong></a>
</p> 

If you store the inertial measurements directly as a raw text file on disk when you collect sensor data, you may need to make them as a rosbag. In this case, you can use [ikalibr-raw-inertial-to-bag](../../launch/tool/ikalibr-raw-inertial-to-bag.launch). After you have configured the contents in the launch file, run:

```sh
roslaunch ikalibr ikalibr-raw-inertial-to-bag.launch
```


<p align="left">
    <a><strong>LiDAR Map Viewer »</strong></a>
</p> 

It's just a little toy, don't pay any attention to it.
