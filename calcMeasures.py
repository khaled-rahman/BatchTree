import sys, math
from collections import defaultdict  
if __name__ == "__main__":
    #run: python calcMeasures.py testgraph.mtx testgraph.mtxBatch2PrEd64PARAOUT25.txt testgraph.label
    #python calcMeasures.py datasets/raw/Graph_8.txt.mtx datasets/output/Graph_8.txt.mtxBatchPrEL128PARAOUT125.txt datasets/raw/Graph_8.txt.labels datasets/raw/Graph_1.txt.labels datasets/raw/Graph_2.txt.labels datasets/raw/Graph_3.txt.labels datasets/raw/Graph_4.txt.labels datasets/raw/Graph_5.txt.labels datasets/raw/Graph_6.txt.labels datasets/raw/Graph_7.txt.labels datasets/raw/Graph_8.txt.labels

    ename = sys.argv[1]
    fname = sys.argv[2]
    lname = sys.argv[3]
    toplabelval = [550, 500, 450, 400, 350, 300, 250, 200]
    #toplabelval = [200, 250, 300, 350, 400, 450, 500, 550]
    toplabels = defaultdict(list)
    arg = 0
    while arg < len(sys.argv) - 4:
        l1file = open(sys.argv[arg+4], "r")
        for line in l1file.readlines():
            toplabels[arg].append(line.strip())
        arg += 1

    cfile = open(fname, "r")
    lfile = open(lname, "r")
    efile = open(ename, "r")
    labels = dict()
    i = 0
    X = []
    xmin = ymin = float("inf")
    xmax = ymax = float("-inf")
    nodes = dict()
    
    for line in cfile.readlines():
        tokens = line.strip().split()
        x = float(tokens[0])
        y = float(tokens[1])
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
        X.append([x, y])
    cfile.close()
    for line in lfile.readlines():
        labels[i] = len(line)
        nodes[i] = line.strip()
        i += 1
    lfile.close()


    l1nodes = dict()
    for line in l1file.readlines():
        l1nodes[line.strip()] = 1
    l1file.close()
    edges = []
    efile.readline() # first commented line
    efile.readline() # node node edge line
    for line in efile.readlines():
        tokens = line.strip().split()
        x = int(tokens[0])-1
        y = int(tokens[1])-1
        for lv in range(len(toplabelval)):
            if nodes[x] in toplabels[lv] and nodes[y] in toplabels[lv]:
                edges.append([x, y, toplabelval[lv]])
                break
        #if nodes[x] in l1nodes and nodes[y] in l1nodes:
        #    edges.append([x, y])
    efile.close()
    ################Area Coverage################
    labelarea = 0
    H = 2 * 1.15
    axislimit = max([abs(xmin), abs(ymin), abs(xmax), abs(ymax)])
    for n in range(len(X)):
        offset = labels[n] * 0.45;
        W = 2 * offset
        labelarea += W * H
        print("Area of :", n, " is:", W * H)
    boundingArea = (ymax - ymin)*(xmax - xmin)
    #boundingArea = 4 * axislimit * axislimit
    print("Total #Edges:", len(edges))

    print("X length:", len(X))
    print("Xmin:", xmin, "Xmax:", xmax, "Ymin:", ymin, "Ymax:", ymax)
    print("Total Label Area:", labelarea, "Bounding Area:", boundingArea)
    print("Area Coverage:", labelarea / boundingArea)

    ###############Desired Edge Length###########
    totaldiff = 0
    for [a, b, c] in edges:
        dist = math.sqrt((X[a][0] - X[b][0]) * (X[a][0] - X[b][0]) + (X[a][1] - X[b][1]) * (X[a][1] - X[b][1]))
        diff = abs(dist - c)
        totaldiff += math.pow(diff / c, 2)
    averagediff = math.sqrt(totaldiff / len(edges))
    print("Edge Length Preservation:", 1.0 - averagediff)
        
