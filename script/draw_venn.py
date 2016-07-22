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

from sys import argv

total_a = int(argv[1])
total_b = int(argv[2])
overlap_num = int(argv[3])

s = (
    #total_a-overlap_num,  # Ab
    #total_b-overlap_num,  # aB
    total_a,
	total_b,
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
    text.set_fontsize(32.0)
for text in v.subset_labels:
    text.set_fontsize(32.0)

# Border styles
c = venn2_circles(subsets=s, linestyle='solid', linewidth=0.0)
#c[0].set_ls('dashed')  # Line style
#c[0].set_lw(2.0)       # Line width

#plt.show()
plt.tight_layout()
plt.savefig('./de_venn.eps', format='eps', dpi=350)
