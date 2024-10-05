<div style="text-align: center;">
    <img src="../img/logo.svg" style="width: 100vw; height: auto;">
</div>

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

+ **IMU** type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/sensor_model.h#L82-L105).
  + `SENSOR_IMU`, `SENSOR_IMU_G`, `SENSOR_IMU_G_NEG`, `SENSOR_IMU_DEG`, `SENSOR_IMU_DEG_G`, `SENSOR_IMU_DEG_G_NEG`.
+ **Radar** Type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/sensor_model.h#L126-L141). ROS message definition can be found [here](https://github.com/Unsigned-Long/iKalibr/tree/master/msg).
  + `AINSTEIN_RADAR`, `AWR1843BOOST_RAW`, `AWR1843BOOST_CUSTOM`, `POINTCLOUD2_POSV`, `POINTCLOUD2_POSIV`, `POINTCLOUD2_XRIO`.
+ **LiDAR** type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/sensor_model.h#L109-L122). The type of LiDAR is more prone to error than others, most of the time error happens when decoding time stamps of points if the wrong LiDAR type is passed in `iKalibr`.
  + *Velodyne LiDARs*: `VLP_16_PACKET`, `VLP_POINTS`.

  + *Ouster LiDARs*: `OUSTER_POINTS`.

  + *Hesai Pandar XT LiDARs*: `PANDAR_XT_POINTS`.

  + *Livox LiDARs*: `LIVOX_CUSTOM`.
+ **Camera** type: corresponding type definition can be found [here](https://github.com/Unsigned-Long/iKalibr/blob/master/include/sensor/sensor_model.h#L46-L78). If you have an RS camera but use the GS camera option, it will not cause a fatal error, but it will reduce the calibration accuracy.
  + `SENSOR_IMAGE_GS`, `SENSOR_IMAGE_COMP_GS`
  +  `SENSOR_IMAGE_RS_FIRST`, `SENSOR_IMAGE_COMP_RS_FIRST`
  + `SENSOR_IMAGE_RS_MID`, `SENSOR_IMAGE_COMP_RS_MID` 
  + `SENSOR_IMAGE_RS_LAST`, `SENSOR_IMAGE_COMP_RS_LAST`
  + for RS cameras, make sure understand how the RS images are timestamped: 
    + `SENSOR_IMAGE_RS_FIRST`: the timestamp of an image is denoted as the time of the first row;
    + `SENSOR_IMAGE_RS_MID`: the timestamp of an image is denoted as the time of the middle row;
    + `SENSOR_IMAGE_RS_LAST`: the timestamp of an image is denoted as the time of the last row;
+ **RGBD Camera** type: same as **Camera** type above. If your optical camera has depth information, you can calibrate it as an RGBD camera.

<p align="left">
    <a><strong>Can Not Find Suitable Sensor Type »</strong></a>
</p> 

If the sensor type you are interested in is not in this list, or the message format of your data collection is different from `iKalibr`, you can:

+ If it is an **official** data format and we do not implement it, please [issue](https://github.com/Unsigned-Long/iKalibr/issues) us.
+ If it is a **custom** data format, please implement the data interface yourself. I believe you can solve it yourself soon.
