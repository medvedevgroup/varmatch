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
    print "three parameters required: file.direct baseline.vcf query.vcf"
    print "0/1 to indicate if it is baseline or query"
    print "\t0 baseline"
    print "\t1 query"
    exit()


matching_filename = argv[1]  # DIRECT this is the same for baseline and query

baseline_final_filename = argv[2]  # baseline.final.vcf
query_final_filename = argv[3]  # query.final.vcf

baseline_original_id_num = {}

with open(baseline_final_filename) as baseline_final_file:
    for line in baseline_final_file.readlines():
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        vcf_id = int(columns[6])
        if vcf_id in baseline_original_id_num:
            baseline_original_id_num[vcf_id] += 1
        else:
            baseline_original_id_num[vcf_id] = 1

query_original_id_num = {}

with open(query_final_filename) as query_final_file:
    for line in query_final_file.readlines():
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        vcf_id = int(columns[6])
        if vcf_id in query_original_id_num:
            query_original_id_num[vcf_id] += 1
        else:
            query_original_id_num[vcf_id] = 1

print len(baseline_original_id_num), len(query_original_id_num)
# by now, only read baseline and query original number for each vcf_id
# vcf_id - num pairs
# only need to check situations in the following code

# record
# baseline_id - matching_num
# query_id - matching_num
# query_id - baseline_id

baseline_matched_id_num = {}
query_matched_id_num = {}
to_ref_id = {}
with open(matching_filename) as matching_file:
    for line in matching_file.readlines():
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        baseline_vcf_id = int(columns[2])
        query_vcf_id = int(columns[4])
        ref_id = int(columns[2])

        if baseline_vcf_id in baseline_matched_id_num:
            baseline_matched_id_num[baseline_vcf_id] += 1
            #to_ref_id[vcf_id].append(ref_id)
        else:
            baseline_matched_id_num[baseline_vcf_id] = 1
            #to_ref_id[vcf_id] = [ref_id]

        if query_vcf_id in query_matched_id_num:
            query_matched_id_num[query_vcf_id] += 1
            to_ref_id[query_vcf_id].append(ref_id)
        else:
            query_matched_id_num[query_vcf_id] = 1
            to_ref_id[query_vcf_id] = [ref_id]


baseline_matched_variant_num = 0
baseline_match_id = {}
baseline_partial_id = {}
# if baseline matched id num == baseline original id num, consider it as an exact match
# if baseline contains complex, this is true
# if baseline contains simples, this is not true, need to process in next function
for vcf_id in sorted(baseline_matched_id_num):
    if vcf_id not in baseline_original_id_num:
        print 'Error: vcf id not exist in original VCF file:', vcf_id
        continue
    if baseline_original_id_num[vcf_id] == baseline_matched_id_num[vcf_id]:
        baseline_matched_variant_num += 1
        baseline_match_id[vcf_id] = True
        #baseline_id_file.write(str(vcf_id)+'\n')
    elif baseline_original_id_num[vcf_id] > baseline_matched_id_num[vcf_id]:
        baseline_partial_id[vcf_id] = True
        #print 'Partial Matching:', vcf_id
    else:
        print 'Error: baseline matched variant number in certain id greater than original variant number:', vcf_id

# for each query matched id
# if corresponding matching baseline id is an exact match, then:
# if query matched id num == query original id num, consider it as an exact match
# if query contains simples, then baseline contains complex, this is true, you do not even have to check equal condition
#   only keeps record of corresponding baseline id set, which is equal to previous baseline match
# if query contains complex, then baseline contains simples:
# equal condition make sure that query is exact match, then query corresponding baseline is exact match,
#    because you do not need to count matching number in baseline

# we do not consider conditions where both baseline and query are complex

query_match_id = {}
query_partial_id = {}

query_matched_variant_num = 0
query_matched_baseline_id = {}



for vcf_id in sorted(query_matched_id_num):
    if vcf_id not in query_original_id_num:
        print 'Error: vcf id not exist in original VCF file:', vcf_id
        continue
    exact_match = True
    for ref_id in to_ref_id[vcf_id]:
        if ref_id not in baseline_match_id:
            exact_match = False
            break
    if not exact_match:
        query_partial_id[vcf_id] = True
        continue
    if query_original_id_num[vcf_id] == query_matched_id_num[vcf_id]:
        query_matched_variant_num += 1
        query_match_id[vcf_id] = True
        for ref_id in to_ref_id[vcf_id]:
            query_matched_baseline_id[ref_id] = True
    elif query_original_id_num[vcf_id] > query_matched_id_num[vcf_id]:
        query_partial_id[vcf_id] = True
        for ref_id in to_ref_id[vcf_id]:
            baseline_partial_id[ref_id] = True
        #print 'Partial Matching:', vcf_id
    else:
        print 'Error: query matched variant number in certain id greater than original variant number:', vcf_id

baseline_id_filename = "baseline_id.txt"
baseline_id_file = open(baseline_id_filename, 'w')
for vcf_id in query_matched_baseline_id:
    baseline_id_file.write(str(vcf_id) + '\n')
baseline_id_file.close()

print 'before query correction:'
print 'baseline match num: ', len(baseline_match_id)
print 'baseline partial num: ', len(baseline_partial_id)
# do not count vcf_id number, count vcf_id set size

print 'after query correction:'
print 'baseline match num: ', len(query_matched_baseline_id)
print 'baseline partial num: ', len(baseline_partial_id)
print 'query match num: ', len(query_match_id)
print 'query partial num: ', len(query_partial_id)


