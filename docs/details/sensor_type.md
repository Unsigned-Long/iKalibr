

<div align=center><img src="../img/logo.svg" width =100%></div>

---

<h3 align="center">How To Select Sensor Types In Configure File</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>


---

When you want to calibrate a sensor suite, you need to write their information in the configure file, such as their topics, their types, and their weights when performing solving. **Selecting a right type for each sensor is important**, as different type points to different data loader to ensure decode the binary data in the rosbag correctly. 

<p align="left">
    <a><strong>Currently Supported Sensor Type »</strong></a>
</p> 

Currently, the following types of sensor are supported in `iKalibr`:

+ **IMU** type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/imu_data_loader.h#L46-L50).
  + `SENSOR_IMU`, `SBG_IMU`, `SENSOR_IMU_G`
+ **Radar** Type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/radar_data_loader.h#L46-L53). ROS message definition can be found [here](https://github.com/Unsigned-Long/iKalibr/tree/master/msg).
  + `AINSTEIN_RADAR`, `AWR1843BOOST_RAW`, `AWR1843BOOST_CUSTOM`, `POINTCLOUD2_POSV`, `POINTCLOUD2_POSIV`, `POINTCLOUD2_XRIO`
+ **LiDAR** type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/lidar_data_loader.h#L56-L71). The type of LiDAR is more prone to error than others, most of the time error happens when decoding time stamps of points if the wrong LiDAR type is passed in `iKalibr`.
  + *Velodyne LiDARs*: `VLP_16_PACKET`, `VLP_16_POINTS`, `VLP_32E_POINTS`

  + *Ouster LiDARs*: `OUSTER_16_POINTS`, `OUSTER_32_POINTS`, `OUSTER_64_POINTS`, `OUSTER_128_POINTS`

  + *Hesai Pandar XT LiDARs*: `PANDAR_XT_16`, `PANDAR_XT_32`

  + *Livox LiDARs*: `LIVOX_MID_360`, `LIVOX_AVIA`
+ **Camera** type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/camera_data_loader.h#L45-L65). If you have an RS camera but use the GS camera option, it will not cause a fatal error, but it will reduce the calibration accuracy.
  + `SENSOR_IMAGE_GS`, `SENSOR_IMAGE_RS_FIRST`, `SENSOR_IMAGE_RS_MID`, `SENSOR_IMAGE_RS_LAST`
  + for RS cameras, make sure understand how the RS images are timestamped: 
    + `SENSOR_IMAGE_RS_FIRST`: the timestamp of an image is denoted as the time of the first row;
    + `SENSOR_IMAGE_RS_MID`: the timestamp of an image is denoted as the time of the middle row;
    + `SENSOR_IMAGE_RS_LAST`: the timestamp of an image is denoted as the time of the last row;
+ **RGBD Camera** type: same as **Camera** type above. If your optical camera has depth information, you can calibrate it as an RGBD camera. Compared with calibrating ordinary optical cameras, calibrate a RGBD camera is faster as no SfM is required.

<p align="left">
    <a><strong>Can Not Find Suitable Sensor Type »</strong></a>
</p> 

If the sensor type you are interested in is not in this list, or the message format of your data collection is different from `iKalibr`, you can:

+ If it is an **official** data format and we do not implement it, please [issue](https://github.com/Unsigned-Long/iKalibr/issues) us.
+ If it is a **custom** data format, please implement the data interface yourself. I believe you can solve it yourself soon.
