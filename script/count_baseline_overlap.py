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

from sys import argv

filename1 = argv[1]
filename2 = argv[2]
decomposed_filename = argv[3] # VCF differ from baseline and query

original_id_num = {}

with open(decomposed_filename) as decomposed_file:
    for line in decomposed_file.readlines():
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        vcf_id = int(columns[6])
        if vcf_id in original_id_num:
            original_id_num[vcf_id] += 1
        else:
            original_id_num[vcf_id] = 1

dict1 = {}

matched_id_num = {}

with open(filename1) as file1:
    for line in file1.readlines():
        line = line.strip()
        vcf_id = int(line)
        if vcf_id in matched_id_num:
            matched_id_num[vcf_id] += 1
        else:
            matched_id_num[vcf_id] = 1

for vcf_id in sorted(matched_id_num):
    if vcf_id not in original_id_num:
        print 'Error: vcf id not exist in original VCF file:', vcf_id
        continue
    if original_id_num[vcf_id] == matched_id_num[vcf_id]:
        dict1[vcf_id] = True

matched_id_num_2 = {}

with open(filename2) as file2:
    for line in file2.readlines():
        line = line.strip()
        vcf_id = int(line)
        if vcf_id in matched_id_num_2:
            matched_id_num_2[vcf_id] += 1
        else:
            matched_id_num_2[vcf_id] = 1

dict2 = {}

for vcf_id in sorted(matched_id_num_2):
    if vcf_id not in original_id_num:
        print 'Error: vcf id not exist in original VCF file:', vcf_id
        continue
    if original_id_num[vcf_id] == matched_id_num_2[vcf_id]:
        dict2[vcf_id] = True

match_num = 0

for vcf_id in dict2:
    if vcf_id in dict1:
        match_num += 1

print len(dict1), len(dict2)
print match_num

