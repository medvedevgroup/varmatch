# remove id after normalization
# add marker back

from sys import argv

input_filename = argv[1]

id_marker_filename = "id_marker.txt"

id_marker = {}

with open(id_marker_filename) as id_marker_file:
    for line in id_marker_file.readlines():
        line = line.strip()
        columns = line.split('\t')
        id_marker[columns[0]] = columns[1]

with open(input_filename) as input_file:
    for line in input_file.readlines():
        line = line.strip()
        if line.startswith('#'):
            print line
            continue
        columns = line.split('\t')
        unique_id = columns[2]
        columns[6] = id_marker[unique_id]
        columns[2] = '.'
        result = '\t'.join(columns)
        print result
