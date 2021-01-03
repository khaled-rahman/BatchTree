import sys, math
  
if __name__ == "__main__":
    #run: python calcMeasures.py testgraph.mtx testgraph.mtxBatch2PrEd64PARAOUT25.txt testgraph.label
    ename = sys.argv[1]
    fname = sys.argv[2]
    lname = sys.argv[3]
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
        i += 1
    lfile.close()
    edges = []
    efile.readline() # first commented line
    efile.readline() # node node edge line
    for line in efile.readlines():
        tokens = line.strip().split()
        x = int(tokens[0])-1
        y = int(tokens[1])-1
        edges.append([x, y])
    efile.close()

    ################Area Coverage################
    labelarea = 0
    H = 4 * 1.5
    axislimit = max([abs(xmin), abs(ymin), abs(xmax), abs(ymax)])
    for n in range(len(X)):
        offset = labels[n] * axislimit * 1.5 / 800.0;
        W = 2 * offset
        labelarea += W * H
        print("Area of :", n, " is:", W * H)
    boundingArea = (ymax - ymin)*(xmax - xmin)
    print("X length:", len(X))
    print("Xmin:", xmin, "Xmax:", xmax, "Ymin:", ymin, "Ymax:", ymax)
    print("Total Label Area:", labelarea, "Bounding Area:", boundingArea)
    print("Area Coverage:", labelarea / boundingArea)

    ###############Desired Edge Length###########
    ideallength = 200.0
    totaldiff = 0
    for [a,b] in edges:
        dist = math.sqrt((X[a][0] - X[b][0]) * (X[a][0] - X[b][0]) + (X[a][1] - X[b][1]) * (X[a][1] - X[b][1]))
        diff = abs(dist - ideallength)
        totaldiff += math.pow(diff / ideallength, 2)
    averagediff = math.sqrt(totaldiff / len(edges))
    print("Edge Length Preservation:", 1.0 - averagediff)
        
