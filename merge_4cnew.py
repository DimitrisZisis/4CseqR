import sys
input_file = sys.argv[1]
in_file = open(input_file, "r")

data = []
current = ["", 0, 0]
last    = ["", 0, 0]

for line in in_file:
    attr = line.strip().split()
    current = attr[:3]
    if current == last:
        data += attr[-3:]
    else:
        if data:
            print "\t".join(data)
        data = attr[:3] + attr[-3:]

    last = current


in_file.close()
