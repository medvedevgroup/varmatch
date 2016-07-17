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

if len(argv) < 4:
    print "three parameters required: file.direct baseline/query.vcf 0/1"
    print "0/1 to indicate if it is baseline or query"
    print "\t0 baseline"
    print "\t1 query"
    exit()

decomposed_filename = argv[2] # VCF differ from baseline and query

matching_filename = argv[1] # DIRECT this is the same for baseline and query

base_query_indicator = int(argv[3]) # 0 is baseline, 1 is query

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


matched_id_num = {}
to_ref_id = {}
with open(matching_filename) as matching_file:
    for line in matching_file.readlines():
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        id_index = 2 + base_query_indicator*2
        vcf_id = int(columns[id_index])
        ref_id = int(columns[2])

        if vcf_id in matched_id_num:
            matched_id_num[vcf_id] += 1
            to_ref_id[vcf_id].append(ref_id)
        else:
            matched_id_num[vcf_id] = 1
            to_ref_id[vcf_id] = [ref_id]

matched_variant_num = 0
matched_baseline_num = 0

partial_match_id = {}

baseline_id_filename = "baseline_id.txt"

baseline_id_file = open(baseline_id_filename, 'w')

for vcf_id in sorted(matched_id_num):
    if vcf_id not in original_id_num:
        print 'Error: vcf id not exist in original VCF file:', vcf_id
        continue
    if original_id_num[vcf_id] == matched_id_num[vcf_id]:
        matched_variant_num += 1
        for ref_id in to_ref_id[vcf_id]:
            matched_baseline_num += 1
            baseline_id_file.write(str(ref_id)+'\n')
    elif original_id_num[vcf_id] > matched_id_num[vcf_id]:
        partial_match_id[vcf_id] = True
        #print 'Partial Matching:', vcf_id
    else:
        print 'Error: matched variant number in certain id greater than original variant number:', vcf_id

baseline_id_file.close()

print 'total match number: ', matched_variant_num
print 'partial match number: ', len(partial_match_id)
print 'matched baseline num: ', matched_baseline_num

