<div align=center><img src="../img/logo.svg" width =100%></div>

---

<h3 align="center">Tutorial: General Calibration Procedure in iKalibr</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>

---

<p align="left">
    <a><strong>Collect Sufficiently Excited Sensor Data »</strong></a>
</p> 

Sufficiently excited motion is required for accurate spatiotemporal calibration in `iKalibr` (to ensure the observability of spatiotemporal parameters). This is a less noticeable but more important point, and we have to be ready at the starting line.

+ todo...

<p align="left">
    <a><strong>Write An Adaptable Configure File »</strong></a>
</p> 
A template of configure file has been provided [here](config_template_note.md). The detailed notes are also provided. We highly recommend you to read it before performing configuring and further solving, if it's your first time.

<p align="left">
    <a><strong>Perform Calibration Using iKalibr »</strong></a>
</p> 

Finally, we can launch the `ikalibr` ros node to perform spatiotemporal determination for a given sensor suite. Just run the following command:

+ ***way** **one***: use `rosrun`, the `roscore` should be ready. The `path_of_your_config_file` is the filename pointing to your own `yaml`-format configure file.

  ```sh
  rosrun ikalibr ikalibr_prog _config_path:="path_of_your_config_file"
  ```

+ **way** **two**: use `roslaunch`. Before launch, you have to pass the filename of your own `yaml`-format configure file to `ikalibr-prog.launch` in directory `ikalibr/launch/solver`.

  ```sh
  roslaunch ikalibr ikalibr-prog.launch
  ```

  Here is the detailed launch file:

  ```xml
  <?xml version="1.0" encoding="UTF-8" ?>
  <launch>
      <node pkg="ikalibr" type="ikalibr_prog" name="ikalibr_prog" output="screen">
          <!-- change the value of this field to the path of your self-defined config file -->
          <param name="config_path" value="path_of_your_config_file" type="string"/>
      </node>
  </launch>
  ```

  

