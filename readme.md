<div align=center><img src="docs/img/logo.svg" width =100%></div>

---

<h3 align="center">iKalibr: Unified Targetless Spatiotemporal Calibration Framework</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>
<p align="center"><i>"The naming of <strong>iKalibr</strong> is inspired by <a href="https://github.com/ethz-asl/kalibr.git">Kalibr</a>, but not an extension of it. Wonder <a href="https://github.com/Unsigned-Long/iKalibr/blob/master/docs/details/why_ikalibr.md">why</a>?"</i></p>
<p align="center">
    :tada: <a href="https://github.com/Unsigned-Long/iKalibr/blob/master/docs/details/news.md"><i><strong>News »</strong> sopport spatial & temporal priori constraints, if needed</i></a>
</p>

---

`iKalibr` is a spatiotemporal calibration framework focusing on resilient integrated inertial systems (sensor suite integrates at least one IMU), the features of `iKalibr` are listed as follows:

+ ***Targetless***: requires no additional artificial targets or facilities. This is perhaps the biggest difference between `iKalibr` and [Kalibr](https://github.com/ethz-asl/kalibr.git) (Kalibr is a chessboard-based visual-inertial calibrator).
+ ***Spatiotemporal***: determines both spatial (extrinsic rotations and translations) and temporal (time offsets, readout time of RS cameras) parameters.
+ ***Resilient and compact***: supports a wide range of sensor suites for one-shot resilient calibration. The <u>IMUs</u>, <u>radars</u>, <u>LiDARs</u>, and <u>optical cameras</u> (both [GS](https://www.arducam.com/global-shutter-camera/) and [RS](https://en.wikipedia.org/wiki/Rolling_shutter) modes) are supported currently. See the following illustration of the full pipeline. "Compact" means that no additional sensors are required to calibrate a given sensor suite.
+ ***Easy to extend***: `iKalibr` is a general spatiotemporal calibration pipeline, and easy to extend to support other sensor types, such as the popular [event](https://en.wikipedia.org/wiki/Event_camera) cameras.

For more details, please refer to our article :point_down::

+ **S. Chen**, X. Li*, S. Li, Y. Zhou, and X. Yang. iKalibr: Unified Targetless Spatiotemporal Calibration for Resilient Integrated Inertial Systems[J]. arXiv:2407.11420 [cs.RO]. [[paper-arXiv](https://arxiv.org/abs/2407.11420)] [[video](https://www.bilibili.com/video/BV1Dm421G7FV/?vd_source=87245258ec5029cca67d77cef1a6201d)]

---

<p align="middle">
    <a href="https://www.bilibili.com/video/BV1Dm421G7FV/?vd_source=87245258ec5029cca67d77cef1a6201d"><strong>« Demo Video of iKalibr (Click To Jump) »</strong></a>
</p> 
<p align="center"><i>"Targetless, Spatial & Temporal, Resilient, Easy To Use, You Only Calibrate Once"</i></p>

---

<div align=center><img src="docs/img/overview.drawio.jpg" width =100%></div>

---

<p align="left">
    <a href="https://github.com/Unsigned-Long/iKalibr/blob/master/docs/details/build_ikalibr.md"><strong>Tutorial: Prepare iKalibr Environment for Calibration »</strong></a>
</p> 


+ install available operation systems and required third libraries.
+ clone `iKalibr` and its submodules on github, compile submodules.
+ compile `iKalibr` (this would require some time).

<p align="left">
    <a href="https://github.com/Unsigned-Long/iKalibr/blob/master/docs/details/use_ikalibr.md"><strong>Tutorial: General Calibration Procedure in iKalibr »</strong></a>
</p> 


+ collect sensor data dynamically, sufficiently excited motion is needed.
+ write adaptable configure file (a template file has been provided).
+ run `iKalibr` based on the configure file.

<p align="left">
    <a href="https://github.com/Unsigned-Long/iKalibr/blob/master/docs/details/tools.md"><strong>Additional Helpful Tools in iKalibr »</strong></a>
</p> 


+ rosbag assembly before solving, such as raw inertial csv file or images to rosbag, merge or split rosbag.
+ data format transformation after solving, some visualization scripts.

<p align="left">
    <a href="https://github.com/Unsigned-Long/iKalibr/blob/master/docs/details/dataset.md"><strong>Dataset Utilized in iKalibr for Evaluation »</strong></a>
</p> 



+ our dataset: two cameras, a Livox Avia lidar (with built-in IMU), a velodyne VLP-32C lidars, two mmWave 3D radars, a MTI IMU.
+ [LI-Calib (OA-Calib)](https://github.com/APRIL-ZJU/lidar_IMU_calib.git) dataset: a velodyne VLP-16 lidars, three IMUs.
+ [River](https://github.com/Unsigned-Long/River.git) dataset: two mmWave 3D radars, a MTI IMU.
+ [TUM GS-RS](https://cvg.cit.tum.de/data/datasets/rolling-shutter-dataset) dataset: a GS camera, a RS camera, a Bosch IMU.

---

<h5 align="center">Copyright</h5>

---

<p align="center"><i>iKalibr: Unified Targetless Spatiotemporal Calibration Framework<br />
    Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China<br />
    https://github.com/Unsigned-Long/iKalibr.git<br /></i></p>
<p align="center"><i>Author: <strong>Shuolong Chen</strong> (shlchen@whu.edu.cn)<br />
    GitHub: https://github.com/Unsigned-Long<br />
    ORCID: 0000-0002-5283-9057<br /></i></p>

*Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:*

* *Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.*
* *Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.*
* *Neither the name of the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.*

*THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*
