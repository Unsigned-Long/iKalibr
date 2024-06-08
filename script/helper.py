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

import json
import os


def get_array_fields(filename, fields):
    array_fields = {}
    file = open(filename, "r")
    lines = file.readlines()
    content = ''
    for line in lines:
        content += line

    array_buffer = json.loads(content)

    for field in fields:
        array_buffer = array_buffer[field]

    # data = {}
    #
    # for i in range(len(array_buffer)):
    #     topic = array_buffer[i]['key']
    #     value = array_buffer[i]['value']
    #     data[topic] = value

    return array_buffer


def get_files_in_dir(dir, format):
    file_paths = []
    for folder, subs, files in os.walk(dir):
        for filename in files:
            file_paths.append(os.path.abspath(os.path.join(folder, filename)))
    files = [file for file in file_paths if file.endswith(format)]
    return files


def sort_param_files(files):
    sorted_files = sorted(files, key=lambda s: int(s[s.rfind('_') + 1:s.rfind('.')]))
    return sorted_files
