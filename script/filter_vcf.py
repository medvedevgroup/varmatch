# this script add marker into VCF files FILTER section at columns 6 (0-based)
# the marker is identical for each VCF entry in the VCF file
from sys import argv

input_filename = argv[1]

with open(input_filename) as input_file:
    for line in input_file.readlines():
        line = line.strip()
        if line.startswith('#'):
            print line
            continue
        columns = line.split('\t')
        ref = columns[3]
        alts = columns[4].split(',')
        genotype = columns[9].split(':')[0]
        if genotype == '0/0' or genotype == '0|0':
            continue
        large_var = False
        for alt in alts:
            if abs(len(ref)-len(alt)) > 100:
                large_var = True
                break
        if large_var:
            continue
        print line
