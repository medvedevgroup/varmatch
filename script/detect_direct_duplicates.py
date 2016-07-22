from sys import argv

filename = argv[1]

match_dict = {}

with open(filename) as file:
    for line in file.readlines():
        line = line.strip()
        if line not in match_dict:
            match_dict[line] = True
        else:
            print "detected duplication"
