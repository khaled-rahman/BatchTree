import networkx as nx
import scipy
import scipy.io
import io, sys, random
from random import choice
from string import ascii_lowercase


if len(sys.argv) > 1:
    num = int(sys.argv[1])
else:
    num = 1000

fw = io.BytesIO()
g = nx.random_tree(n=num, seed=0)
mtx = nx.to_scipy_sparse_matrix(g)
scipy.io.mmwrite(fw, mtx)
tree_file = open("synthetictree_"+str(num)+".mtx", "w")
tree_file.write(fw.getvalue().decode('utf-8'))
tree_file.close()

label_file = open("syntheticlabel_"+str(num)+".labels", "w")
for i in range(num):
    length = random.randint(2, 10)
    label = ''.join(choice(ascii_lowercase) for i in range(length))
    label_file.write(label + "\n")
label_file.close()

