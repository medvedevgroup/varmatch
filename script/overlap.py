#!/usr/bin/env python

# Copyright 2015, Chen Sun
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
    Authors:
    Chen Sun(chensun@cse.psu.edu)
    Paul Medvedev(pashadag@cse.psu.edu)
"""

import sys
import textwrap as _textwrap
import argparse
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

versionError = "You are using an old version of python, please upgrade to python 2.7+\n"
if sys.hexversion < 0x02070000:
    print (versionError)
    exit()

RUN = True

author_email = 'chensun@cse.psu.edu'

class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        paragraphs = text.split('\n')
        #return paragraphs
        multiline_text = []
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.wrap(paragraph, width)
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

    def _fill_text(self, text, width, indent):
        return ''.join(indent + line for line in text.splitlines(True))

citation = 'Please cite our paper.'

parser = argparse.ArgumentParser(prog="overlap", epilog = citation, formatter_class=lambda prog: SmartFormatter(prog,max_help_position=8))
parser.add_argument('match_files', nargs='+', metavar='File List', help='.match file list, should take the same VCF as baseline')
args = parser.parse_args()

baselinevar_num = {}
# key is baseline variant, identified by chr_pos_ref_alt_phasing
variant_num_list = []


def read_file(filename):
    global baselinevar_num
    global variant_num
    variant_num = 0
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip()
            columns = line.split('\t')
            chrname = columns[0]
            baseline_columns = columns[4].split(';')
            if baseline_columns[0] == '.':
                baseline_key = chrname + ',' + columns[1]# + ',' + columns[2] + ',' + columns[3]
                #print baseline_key
                if baseline_key in baselinevar_num:
                    baselinevar_num[baseline_key] += 1
                else:
                    baselinevar_num[baseline_key] = 1
                variant_num += 1
                continue

            for baseline in baseline_columns:
                baseline_key = chrname + ',' + baseline.split(',')[0]
                #print baseline_key
                if baseline_key in baselinevar_num:
                    baselinevar_num[baseline_key] += 1
                else:
                    baselinevar_num[baseline_key] = 1
                variant_num += 1
    variant_num_list.append(variant_num)

for filename in args.match_files:
    read_file(filename)

overlap_threshold = len(args.match_files)

print overlap_threshold
overlap_num = 0
for baseline_key in baselinevar_num:
    if baselinevar_num[baseline_key] >= overlap_threshold:
        overlap_num += 1

print overlap_num, variant_num_list

# Subset sizes
s = (
    variant_num_list[0]-overlap_num,  # Ab
    variant_num_list[1]-overlap_num,  # aB
    overlap_num,  # AB
)

v = venn2(subsets=s, set_labels=('bwa-fb', 'pt'))

# Subset labels
v.get_label_by_id('10').set_text(format(s[0], ',d'))
v.get_label_by_id('01').set_text(format(s[1], ',d'))
v.get_label_by_id('11').set_text(format(s[2], ',d'))

# Subset colors
#v.get_patch_by_id('10').set_color('red')
#v.get_patch_by_id('01').set_color('yellow')
#v.get_patch_by_id('11').set_color('blue')

# Subset alphas
#v.get_patch_by_id('10').set_alpha(0.4)
#v.get_patch_by_id('01').set_alpha(1.0)
#v.get_patch_by_id('11').set_alpha(0.7)

for text in v.set_labels:
    text.set_fontsize(32)
for text in v.subset_labels:
    text.set_fontsize(32)

# Border styles
c = venn2_circles(subsets=s, linestyle='solid', linewidth='0')
#c[0].set_ls('dashed')  # Line style
#c[0].set_lw(2.0)       # Line width

#plt.show()
plt.tight_layout()
plt.savefig('./vm_venn.png')
