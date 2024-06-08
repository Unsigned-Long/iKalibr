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

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.patches import Rectangle
import math
import matplotlib.path as mpath

config = {
    # "text.usetex": True,
    "font.family": 'Times New Roman',  # sans-serif/serif/cursive/fantasy/monospace
    "font.size": 18.0,  # medium/large/small
    'font.style': 'normal',  # normal/italic/oblique
    'font.weight': 'bold',  # bold
    "mathtext.fontset": 'cm',
    "font.serif": ['cmb10'],
    "axes.unicode_minus": False,
}

plt.rcParams.update(config)


class drawer:

    @staticmethod
    def show_figure(save_path=None):
        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path)
        plt.show()

    @staticmethod
    def set_fig_size(width, height):
        plt.rcParams['figure.figsize'] = (width, height)

    @staticmethod
    def math_symbols(s: str):
        return r'$' + s + '$'

    @staticmethod
    def set_xticks(ax, xmin, xmax, tick_count):
        ax.set_xticks(np.arange(xmin, xmax + 1E-10, (xmax - xmin) / tick_count))
        ax.set_xlim(xmin, xmax)

    @staticmethod
    def set_yticks(ax, ymin, ymax, tick_count):
        ax.set_yticks(np.arange(ymin, ymax + 1E-10, (ymax - ymin) / tick_count))
        ax.set_ylim(ymin, ymax)

    @staticmethod
    def set_zticks(ax, zmin, zmax, tick_count):
        ax.set_zticks(np.arange(zmin, zmax + 1E-10, (zmax - zmin) / tick_count))
        ax.set_zlim(zmin, zmax)

    @staticmethod
    def add_grids(ax, alpha=0.5, axis='both'):
        ax.grid(ls='-.', alpha=alpha, axis=axis)

    @staticmethod
    def add_grids_3d(ax, alpha=0.5):
        ax.xaxis._axinfo["grid"]['linestyle'] = "-."
        ax.yaxis._axinfo["grid"]['linestyle'] = "-."
        ax.zaxis._axinfo["grid"]['linestyle'] = "-."

        ax.xaxis._axinfo["grid"]['alpha'] = alpha
        ax.yaxis._axinfo["grid"]['alpha'] = alpha
        ax.zaxis._axinfo["grid"]['alpha'] = alpha

    @staticmethod
    def add_coordinate(ax, x, y, z, arrow_length=1, line_width=3, mark_size=5, marker='8'):
        # draw x, y, z axis
        ax.plot([x, x + arrow_length], [y, y], [z, z], linestyle='-', ms=mark_size, marker=marker, c='red',
                lw=line_width)

        ax.plot([x, x], [y, y + arrow_length], [z, z], linestyle='-', ms=mark_size, marker=marker, c='green',
                lw=line_width)

        ax.plot([x, x], [y, y], [z, z + arrow_length], linestyle='-', ms=mark_size, marker=marker, c='blue',
                lw=line_width)

        ax.plot([x], [y], [z], linestyle='-', ms=mark_size, marker=marker, c='black', lw=line_width)

    @staticmethod
    def set_label_decimal(ax, format_str, axis='y'):
        if axis == 'y':
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter(format_str))
        elif axis == 'x':
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter(format_str))
        elif axis == 'z':
            ax.zaxis.set_major_formatter(mtick.FormatStrFormatter(format_str))

    @staticmethod
    def set_legend_top(ax, cols):
        ax.legend(bbox_to_anchor=(0.0, 1.05, 1.0, 0.0), loc=3, ncol=cols, mode="expand", borderaxespad=0.0)

    @staticmethod
    def get_cmap(color_name, len, start=0.0, end=1.0):
        """
        Perceptually Uniform Sequential:
        ['viridis', 'plasma', 'inferno', 'magma', 'cividis']

        Sequential:
        ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
         'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

        Sequential(2):
        ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
         'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper']

        Diverging:
        ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

        Cyclic:
        ['twilight', 'twilight_shifted', 'hsv']

        Qualitative:
        ['Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c']

        Miscellaneous：
        ['flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
         'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral', 'gist_ncar']
        """

        return plt.get_cmap(color_name)(np.linspace(start, end, len))

    @staticmethod
    def set_sci_label(ax, axis='y'):
        ax.ticklabel_format(style='sci', scilimits=(0, 0), axis=axis)

    @staticmethod
    def use_style(idx=0):
        styles = ['default', 'classic'] + sorted(
            style for style in plt.style.available if style != 'classic' and not style.startswith('_')
        )
        plt.style.use(styles[idx])

    @staticmethod
    def landmark_marker():
        s1 = [[-1, 0], [0, 0], [1, 0]]
        s2 = [[0, -1], [0, 0], [0, 1]]
        s3 = [[-math.sqrt(2) / 2, -math.sqrt(2) / 2], [0, 0], [math.sqrt(2) / 2, math.sqrt(2) / 2]]
        s4 = [[math.sqrt(2) / 2, -math.sqrt(2) / 2], [0, 0], [-math.sqrt(2) / 2, math.sqrt(2) / 2]]
        c = [1, 2, 2]
        marker = mpath.Path(vertices=np.concatenate([s1, s2, s3, s4]), codes=np.concatenate([c, c, c, c]))
        return marker

    @staticmethod
    def set_axis_label_color(ax, color, axis='y'):
        if axis == 'y':
            ax.yaxis.get_label().set_color(color)
            # ax.yaxis.get_ticks().set_color(color)
            ax.tick_params(axis, colors=color)
        elif axis == 'x':
            ax.xaxis.get_label().set_color(color)

    @staticmethod
    def erase_ax(ax, style_idx=1):
        ax.set_yticks([])
        x_min, x_max = ax.get_xlim()
        hatches = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']
        ax.add_patch(Rectangle((x_min, 0), x_max - x_min, 1, fill=False, hatch=hatches[1]))


if __name__ == '__main__':
    drawer.set_fig_size(12.0, 4.0)

    host = host_subplot(121)

    twin = host.twinx()

    # fig, axs = plt.subplots(nrows=1, ncols=2, sharex='all', sharey='all')
    t = np.arange(0.0, 1.0 + 0.01, 0.01)
    s = np.sin(4 * np.pi * t)
    tan = np.tan(4 * np.pi * t)

    colors = drawer.get_cmap('plasma', 2, 0.3, 0.8)

    sin_plot, = host.plot(t, s, label=drawer.math_symbols('y=sin(4\pi t)'), color=colors[0])
    tan_plot, = twin.plot(t, tan, label=drawer.math_symbols('y=tan(4\pi t)'), c=colors[1])

    host.set_xlabel('t' + drawer.math_symbols('(s)'))
    host.set_ylabel('y' + drawer.math_symbols('(m)'))
    twin.set_ylabel('y' + drawer.math_symbols('(m)'))

    drawer.set_xticks(host, 0.0, 1.0, 8)
    drawer.set_yticks(host, -2.0, 2.0, 5)
    drawer.set_yticks(twin, -20.0, 20.0, 5)

    drawer.add_grids(host)

    drawer.set_axis_label_color(host, sin_plot.get_color())
    drawer.set_axis_label_color(twin, tan_plot.get_color())

    # drawer.set_sci_label(axs[0])
    drawer.set_label_decimal(host, '%.1f', axis='x')

    drawer.set_legend_top(host, 2)

    drawer.erase_ax(host_subplot(122))

    drawer.show_figure()
