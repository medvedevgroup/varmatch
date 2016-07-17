# this script add marker into VCF files FILTER section at columns 6 (0-based)
# the marker is identical for each VCF entry in the VCF file
from sys import argv

input_filename = argv[1]

unique_id = 1
with open(input_filename) as input_file:
    for line in input_file.readlines():
        line = line.strip()
        if line.startswith('#'):
            print line
            continue
        columns = line.split('\t')
        columns[6] = str(unique_id)
        unique_id += 1
        result = '\t'.join(columns)
        print result
