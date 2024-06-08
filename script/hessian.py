#  iKalibr: Unified Targetless Spatiotemporal Calibration Framework
#  Copyright 2024, the School of Geodesy and Geomatics (SGG), Wuhan University, China
#  https://github.com/Unsigned-Long/iKalibr.git
#
#  Author: Shuolong Chen (shlchen@whu.edu.cn)
#  GitHub: https://github.com/Unsigned-Long
#   ORCID: 0000-0002-5283-9057
#
#  Purpose: See .h/.hpp file.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * The names of its contributors can not be
#    used to endorse or promote products derived from this software without
#    specific prior written permission.
#
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
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from helper import get_array_fields
from plt_utils import drawer

# ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

filename = '/home/csl/ros_ws/iKalibr/src/ikalibr/data/li-calib/Court-01/hessian/hessian.json'
# size = -1
size = 30

colors = ['#ee1d23', '#3b4ba8', '#231f20', '#97D077', '#B5739D']
z_color = 'w'
p_color = colors[3]
n_color = colors[4]

line_color = 'k'
line_alpha = 0.4
line_width = 1.5


def recovery_mat(data, row, col):
    mat = [[] for elem in range(row)]
    idx = 0
    for v in data.values():
        mat[idx % col].append(v)
        idx += 1
    return mat


if __name__ == '__main__':
    data = get_array_fields(filename, [])
    rows = data['row']
    cols = data['col']
    HMat = data['hessian']
    order_size_raw = data['par_order_size']
    order_size = []
    for elem in order_size_raw:
        order_size.append([elem['first'], elem['second']])

    hessian = recovery_mat(HMat, rows, cols)
    for i in range(len(hessian)):
        for j in range(len(hessian[i])):
            val = hessian[i][j]
            if val > 0.0:
                hessian[i][j] = 1
            elif val < 0.0:
                hessian[i][j] = -1

    drawer.set_fig_size(10.0, 10.0)
    fig, ax = plt.subplots(1, 1)

    ax.xaxis.set_ticks_position('top')
    ax.invert_yaxis()
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    custom_cmap = ListedColormap([n_color, z_color, p_color])
    if size <= 0:
        size = rows
    ax.pcolormesh(np.array(hessian)[:size, :size], cmap=custom_cmap, shading='flat')

    for i in range(size):
        ax.vlines(i, 0, size, colors='w', lw=line_width)
        ax.hlines(i, 0, size, colors='w', lw=line_width)

    pos = 0
    for elem in order_size:
        ax.vlines(pos, 0, size, colors=line_color, alpha=line_alpha, lw=line_width)
        ax.hlines(pos, 0, size, colors=line_color, alpha=line_alpha, lw=line_width)

        pos += elem[1]

        if pos > size:
            pos -= elem[1]
            break
        print(elem, pos)

    ax.vlines(pos, 0, size, colors=line_color, alpha=line_alpha, lw=line_width)
    ax.hlines(pos, 0, size, colors=line_color, alpha=line_alpha, lw=line_width)

    ax.set_box_aspect(1)
    ax.set_title('Hessian matrix')

    base, ext = os.path.splitext(filename)
    save_filename = base + '.png'
    drawer.show_figure(save_filename)
