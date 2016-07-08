from sys import argv

baseline_filename = argv[1]
query_filename = argv[2]

baseline_pos_line = {}
baseline_pos_content = {}

query_pos_line = {}
query_pos_content = {}

def read_file(filename):
    pos_line = {}
    pos_content = {}
    with open (filename) as f:
        for line in f:
            if(line.startswith('#')):
                continue
            line = line.strip()
            columns = line.split('\t')
            content = '\t'.join(columns[1:-1])
            pos = int(columns[1])
            pos_line[pos] = line
            pos_content[pos] = content
    return pos_line, pos_content


(baseline_pos_line, baseline_pos_content)= read_file(baseline_filename)
(query_pos_line, query_pos_content)= read_file(query_filename)

for pos in baseline_pos_content:
    if pos not in query_pos_content:
        print pos, "exist in baseline but not in query"
        print baseline_pos_line[pos]
    else:
        if baseline_pos_content[pos] != query_pos_content[pos]:
            print "same pos but not equal content"
            print baseline_pos_line[pos]
            print query_pos_line[pos]

for pos in query_pos_content:
    if pos not in baseline_pos_content:
        print pos, "in query but not in baseline"
        print query_pos_line[pos]