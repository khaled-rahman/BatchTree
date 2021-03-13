import sys

if __name__ == "__main__":
    
    #python processRawWeight.py Graph_1.txt Graph_2.txt Graph_3.txt Graph_4.txt Graph_5.txt Graph_6.txt Graph_7.txt Graph_8.txt

    arg = 0
    labels = dict()
    i = 0
    nodes = dict()
    edges = []
    taken = dict()
    toplabelval = [550, 500, 450, 400, 350, 300, 250, 200]
    while arg < len(sys.argv) - 1:
        fname = sys.argv[arg+1]
        lfile = open(fname, "r")
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

            if labels[x] > labels[y]:
                x, y = y, x

            if (x,y) not in taken:
                edges.append([labels[x], labels[y], toplabelval[arg]])
                taken[(x,y)] = toplabelval[arg]
        lfile.close()
        arg += 1

    ofile = open(fname+".weighted.mtx", "w")
    nfile = open(fname+".labels", "w")
    ofile.write("%%MatrixMarket matrix coordinate pattern symmetric\n")
    ofile.write(str(len(labels)) + " " + str(len(labels)) + " " + str(len(edges)) + "\n")
    for [a, b, w] in edges:
        ofile.write(str(a) + " " + str(b) + " " + str(w) + "\n")
    for v in sorted(labels.values()):
        nfile.write(str(nodes[v][:16]) + "\n")
    ofile.close()
    nfile.close()
