#  iKalibr: Unified Targetless Spatiotemporal Calibration Framework
#  Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China
#  https://github.com/Unsigned-Long/iKalibr.git
#  Author: Shuolong Chen (shlchen@whu.edu.cn)
#  GitHub: https://github.com/Unsigned-Long
#   ORCID: 0000-0002-5283-9057
#  Purpose: See .h/.hpp file.
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * The names of its contributors can not be
#    used to endorse or promote products derived from this software without
#    specific prior written permission.
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
import math

file_path = '/home/csl/ros_ws/iKalibr/src/ikalibr/debug/event_local_planes78.yaml'

import yaml
from scipy.spatial import KDTree

def parse_planes_from_yaml(file_path):
    with open(file_path, "r") as file:
        data = yaml.safe_load(file)

    planes_data = []
    event_local_planes = data.get("event_local_planes", [])

    for plane in event_local_planes:
        first = plane["first"]
        A, B, C = first["r0c0"], first["r1c0"], first["r2c0"]

        second = plane["second"]
        points = [
            (point["tuple_element0"], point["tuple_element1"], point["tuple_element2"])
            for point in second
        ]

        planes_data.append({
            "parameters": {"A": A, "B": B, "C": C},
            "points": points
        })

    return planes_data


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from plt_utils import drawer


def visualize_plane_and_samples(ax, plane):
    A, B, C = plane["parameters"].values()
    points = np.array(plane["points"])

    x_points = points[:, 0]
    y_points = points[:, 1]
    t_points = points[:, 2]

    x_range = np.linspace(min(x_points), max(x_points), 50)
    y_range = np.linspace(min(y_points), max(y_points), 50)
    x_grid, y_grid = np.meshgrid(x_range, y_range)
    t_grid = -(A * x_grid + B * y_grid + C)

    ax.plot_surface(x_grid, y_grid, t_grid, alpha=0.5, rstride=100, cstride=100, edgecolor='none')
    ax.scatter(x_points, y_points, t_points, color='red', label="Sample Points")

    distances = []
    foot_points = []
    dist_sum = 0.0
    for x, y, t in zip(x_points, y_points, t_points):
        # f = A * x + B * y + t + C
        # denominator = A ** 2 + B ** 2 + 1
        # distance = abs(f) / np.sqrt(denominator)
        t_pred = -(A * x + B * y + C)
        dist = math.fabs(t_pred - t)
        distances.append(dist)
        dist_sum += dist

        # x_foot = x - (A * f) / denominator
        # y_foot = y - (B * f) / denominator
        # t_foot = t - f / denominator
        x_foot = x
        y_foot = y
        t_foot = t_pred
        foot_points.append((x_foot, y_foot, t_foot))
    print("   points: {}".format(len(points)))
    print("     dist: {}".format(dist_sum))
    print(" avg dist: {}".format(dist_sum / len(points)))

    for (x, y, t), (x_foot, y_foot, t_foot) in zip(points, foot_points):
        ax.plot([x, x_foot], [y, y_foot], [t, t_foot], color='blue', linestyle='-')

    for i, (x, y, t) in enumerate(points):
        ax.text(x, y, t, f"{distances[i] * 1E3:.2f}", color="black", fontsize=8)

    center_x = np.mean(x_points)
    center_y = np.mean(y_points)
    center_t = -(A * center_x + B * center_y + C)
    ax.scatter([center_x], [center_y], [center_t], color='k', label="Center Points")

    n = A ** 2 + B ** 2
    norm_flow_vec = np.array([-A / n, -B / n, 0.0])
    # print(np.linalg.norm(norm_flow_vec))
    # normal_vector = np.array([A, B, 1])
    norm_flow_vec /= np.linalg.norm(norm_flow_vec)
    normal_start = np.array([center_x, center_y, center_t])
    normal_end = normal_start + norm_flow_vec
    x_coords = [normal_start[0] - norm_flow_vec[0], normal_end[0]]
    y_coords = [normal_start[1] - norm_flow_vec[1], normal_end[1]]
    z_coords = [normal_start[2] - norm_flow_vec[2], normal_end[2]]
    ax.plot(x_coords, y_coords, z_coords, color='green', linewidth=2, label="Norm Flow Vector")
    # ax.quiver(
    #     normal_start[0], normal_start[1], normal_start[2],
    #     norm_flow_vec[0], norm_flow_vec[1], norm_flow_vec[2],
    #     color='green', length=2.0, normalize=True, label="Norm Flow Vector"
    # )


if __name__ == "__main__":
    planes = parse_planes_from_yaml(file_path)

    for i, plane in enumerate(planes):

        A, B, C = plane["parameters"].values()
        n = A ** 2 + B ** 2
        norm_flow_vec = np.array([-A / n, -B / n, 0.0])
        norm = np.linalg.norm(norm_flow_vec)
        if norm < 2000:
            continue
        print(f"Plane {i + 1}:")

        print("  nf norm: {}".format(norm))
        # continue
        # print(f"  Parameters: {plane['parameters']}")
        # print(f"  Points: {plane['points']}")
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        visualize_plane_and_samples(ax, plane)

        ax.set_title("Visualization of Fitted Plane and Sample Points", fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel("x (pixel)")
        ax.set_ylabel("y (pixel)")
        ax.set_zlabel("t (sec)")
        ax.legend()
        plt.show()
        print()
