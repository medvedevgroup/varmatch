# this script add id into VCF files ID section at columns 2 (0-based)
# also output id-marker pair in id_marker.txt
# id is unique, marker is not unique because of decomposition
from sys import argv

input_filename = argv[1]

output_filename = "id_marker.txt"

output = open(output_filename, 'w')

unique_id = 1
with open(input_filename) as input_file:
    for line in input_file.readlines():
        line = line.strip()
        if line.startswith('#'):
            print line
            continue
        columns = line.split('\t')
        columns[2] = str(unique_id)
        marker = columns[6]
        output.write(columns[2] + '\t' + columns[6] + '\n')
        columns[6] = '.'
        unique_id += 1
        result = '\t'.join(columns)
        print result

output.close()

