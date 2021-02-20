import sys

if __name__ == "__main__":
    fname = sys.argv[1]
    lfile = open(fname, "r")
    labels = dict()
    i = 0
    edges = []
    nodes = dict()
    for line in lfile.readlines():
        tokens = line.strip().split("--")
        x = tokens[0].strip().replace("\"", "")
        y = tokens[1].strip().replace("\"", "")
        if x not in labels:
            i += 1
            labels[x] = i
            nodes[i] = x
        if y not in labels:
            i += 1
            labels[y] = i
            nodes[i] = y
        edges.append([labels[x], labels[y]])
    lfile.close()
    ofile = open(fname+".mtx", "w")
    nfile = open(fname+".labels", "w")
    ofile.write("%%MatrixMarket matrix coordinate pattern symmetric\n")
    ofile.write(str(len(labels)) + " " + str(len(labels)) + " " + str(len(edges)) + "\n")
    for [a,b] in edges:
        ofile.write(str(a) + " " + str(b) + "\n")
    for v in sorted(labels.values()):
        nfile.write(str(nodes[v][:16]) + "\n")
    ofile.close()
    nfile.close()
        
