<div align=center><img src="../img/logo.svg" width =100%></div>

---

<h3 align="center">Lastest News About iKalibr</h3>
<p align="center">
    <a href="https://github.com/Unsigned-Long"><strong>Author » Shuolong Chen</strong></a>
</p>



---

<p align="left">
    <a><strong>Version 1.2.0 » Compatible With GLOMAP For SfM Reconstruction</strong></a>
</p>

`GLOMAP` is a general purpose global structure-from-motion pipeline for image-based reconstruction. `GLOMAP` requires a `COLMAP` database as input and outputs a `COLMAP` sparse reconstruction. As compared to `COLMAP`, this project provides a much more efficient and scalable reconstruction process, *typically 1-2 orders of magnitude faster, with on-par or superior reconstruction quality*.

Below is a comparison of SfM reconstruction between colmap mapper and glomap mapper, focusing on reconstruction effect and reconstruction speed.

|     300 IMAGES |              COLMAP MAPPER              |        GLOMAP MAPPER (RECOMMAND)        |
| -------------: | :-------------------------------------: | :-------------------------------------: |
| **START TIME** |             17:04:09.440652             |             16:56:57.489804             |
|   **END TIME** |             17:21:00.542973             |             17:00:54.970989             |
|   **RUN TIME** |            16.852 [minutes]             |           **3.570 [minutes]**           |
|  **FINAL MAP** | ![COLMAP MAPPER](../img/colmap_sfm.png) | ![GLOMAP MAPPER](../img/glomap_sfm.png) |



**Attention:** If cameras are integrated in sensor suite to be calibrated by `iKalibr`, structure from motion (SfM) is required for each camera.



<p align="left">
    <a><strong>Version 1.1.0 » Support Spatial & Temporal Priori Constraints</strong></a>
</p>

For a multi-sensor kit to be calibrated, if *the user already knows any spatiotemporal parameters between some sensors*, the user can pass these parameters into `iKalibr` by configuring this file, and this would treat the spatiotemporal calibration problem in `iKalibr` as an equality constraint optimization problem with prior information, to ensure that the solution meets the prior spatiotemporal parameters.

A typical example of using this file is to calibrate a multi-sensor kit containing `Livox Avia`. Since the manufacturer has provided the extrinsics of the IMU relative to the LiDAR in Avia, these extrinsics can be passed into ikalibr as a priori. 

*"The origin O' of IMU coordinate is defined in the point cloud coordinates as (-41.65, -23.26, 28.40) (Unit: mm).", from **Livox AVIA User Manual***.

To pass the spatiotemporal priori, please edit the [config file](../../config/spat-temp-priori.yaml) and give its path to the `SpatTempPrioriPath` field in the main config file of `iKalibr`.

<div align=center><img src="../../docs/img/priori_constraints.jpg" width =50%></div>

**Attention:** 

+ Please note that this prior knowledge is not necessary in `iKalibr`. If *you have them and are very sure that they are correct*, then provide it to `iKalibr` through this file. If you don't have them, don't need to provide  the corresponding configure file.
+ The new feature in `iKalibr` (support for prior constraints) is theoretically a nonlinear least squares problem with equality constraints. Technically, it can be implemented through [Augmented Lagrangian](https://en.wikipedia.org/wiki/Augmented_Lagrangian_method) or [Sequential Quadratic Programming (SQP)](https://en.wikipedia.org/wiki/Sequential_quadratic_programming). Unfortunately, `Ceres` does not currently support this type of constrained optimization problem. Therefore, in terms of implementation, we directly treat this prior constraint as a residual with a large weight, which is also the suggestion given by `Ceres` developers.

