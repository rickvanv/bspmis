import networkx as nx
from pyvis.network import Network
import fast_matrix_market as fmm

matrix = fmm.mmread("ca-GrQc-test.mtx")
matnet = Network()
size = 5242

newmatrix = matrix.tocsr()

for k in range(0, size):
    matnet.add_node(k)

for i in range(0, size):
    for j in range(0, size):
        if (newmatrix[i, j] == 1.0):
            matnet.add_edge(i, j)
    print(i)

matnet.show("ca-GrQc.html")
