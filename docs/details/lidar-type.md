<div align=center><img src="../img/logo.svg" width =100%></div>

---

<h3 align="center">How To Select Sensor Types In Configure File</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>

---


When you want to calibrate a sensor suite that integrates LiDARs, you need to write their information in the configure file, such as their topics, their types, and their weights when performing solving.

```yaml
LiDARTopics:
  # the 'key' is exactly the message topic
  - key: "/lidar0/scan"
    value:
      Type: "VLP_32E_POINTS"
      Weight: 100.0
  - key: "/lidar1/scan"
    value:
      Type: "LIVOX_MID_360"
      Weight: 80.0
  ...
```

Selecting a right type for each LiDAR is important, as different type points to different data loader to ensure decode the binary data in the rosbag correctly. Currently, the following types of LiDARs are supported in `iKalibr` (see corresponding [codes](https://github.com/Unsigned-Long/iKalibr/blob/master/src/sensor/lidar_data_loader.cpp#L30-L52)):

+ Velodyne LiDARs: [VLP_16_PACKET](https://github.com/Unsigned-Long/iKalibr/blob/master/src/sensor/lidar_data_loader.cpp#L82-L190), [VLP_16_POINTS](https://github.com/Unsigned-Long/iKalibr/blob/master/src/sensor/lidar_data_loader.cpp#L410-L510), [VLP_32E_POINTS](https://github.com/Unsigned-Long/iKalibr/blob/master/src/sensor/lidar_data_loader.cpp#L410-L510)
+ Ouster LiDARs: [OUSTER_16_POINTS](https://github.com/Unsigned-Long/iKalibr/blob/master/include/util/cloud_define.hpp#L61-L69), [OUSTER_32_POINTS](https://github.com/Unsigned-Long/iKalibr/blob/master/include/util/cloud_define.hpp#L61-L69), [OUSTER_64_POINTS](https://github.com/Unsigned-Long/iKalibr/blob/master/include/util/cloud_define.hpp#L61-L69), [OUSTER_128_POINTS](https://github.com/Unsigned-Long/iKalibr/blob/master/include/util/cloud_define.hpp#L61-L69)
+ Hesai Pandar XT LiDARs: [PANDAR_XT_16](https://github.com/Unsigned-Long/iKalibr/blob/master/include/util/cloud_define.hpp#L79-L86), [PANDAR_XT_32](https://github.com/Unsigned-Long/iKalibr/blob/master/include/util/cloud_define.hpp#L79-L86)
+ Livox LiDARs: [LIVOX_MID_360](https://github.com/Unsigned-Long/iKalibr/blob/master/msg/LivoxCustomPoint.msg#L1-L10), [LIVOX_AVIA](https://github.com/Unsigned-Long/iKalibr/blob/master/msg/LivoxCustomPoint.msg#L1-L10)

If your LiDAR type is not included in the currently supported types, please issue us.
