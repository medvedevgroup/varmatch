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
    Paul Medvedev(pashadag@cse.psu.edu)
    Chen Sun(chensun@cse.psu.edu)
"""

import sys
import textwrap as _textwrap
import argparse
import os
from intervaltree import Interval, IntervalTree

versionError = "You are using an old version of python, please upgrade to python 2.7+\n"
if sys.hexversion < 0x02070000:
    print (versionError)
    exit()

RUN = True

author_email = 'chensun@cse.psu.edu'

chr_2_index = {}
index_2_chr = {}

for i in range(22):
    chr_2_index[str(i+1)] = i
    index_2_chr[i] = str(i+1)

chr_2_index['X'] = 22
index_2_chr[22] = 'X'

chr_2_index['Y'] = 23
index_2_chr[23] = 'Y'


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

parser = argparse.ArgumentParser(prog="filter_hc", epilog = citation, formatter_class=lambda prog: SmartFormatter(prog,max_help_position=8))
parser.add_argument('vcf_files', nargs='+', metavar='File List', help='VCF file list')
parser.add_argument('-b', '--bed_file', help='bed hc file', metavar='File')
parser.add_argument('-o', '--output', help='output directory', metavar='DIRECTORY',default='./nhc')
args = parser.parse_args()

if not os.path.exists(args.output):
    os.mkdir(args.output)

interval_tree_list = []

for i in range(24):
    interval_tree_list.append(IntervalTree())

with open (args.bed_file) as bed:
    for line in bed.readlines():
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        chr_name = columns[0]
        if chr_name not in chr_2_index:
            print 'BED: ' + line
            continue
        chr_index = chr_2_index[chr_name]
        start_p = int(columns[1])
        end_p = int(columns[2]) # 0 based, exclude end position
        interval_tree_list[chr_index][start_p: end_p] = (start_p, end_p)

for vcf_filename in args.vcf_files:
    match_basename = os.path.basename(vcf_filename)
    nhc_filename = args.output + '/' + match_basename + '.nhc.vcf'

    output_list = []
    with open(vcf_filename) as vcf_file:
        for line in vcf_file.readlines():
            if line.startswith('#'):
                output_list.append(line)
                continue
            columns = line.split('\t')
            chr_name = columns[0]
            if chr_name not in chr_2_index:
                output_list.append(line)
                continue
            chr_index = chr_2_index[chr_name]
            var_pos = int(columns[1]) - 1 # 1 based system to 0 based system
            query_result = interval_tree_list[chr_index][var_pos]
            if len(query_result) == 0:
                output_list.append(line)

    nhc_file = open(nhc_filename, 'w')
    for line in output_list:
        nhc_file.write(line)
    nhc_file.close()