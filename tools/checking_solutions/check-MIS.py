import networkx as nx
from pyvis.network import Network
import fast_matrix_market as fmm

matrix = fmm.mmread("ca-GrQc-test.mtx")
newmatrix = matrix.tocsr()

size = matrix.shape[0]
print(size)


mis_raw_in = input(
    "Please enter your computed MIS to check if it is correct: ")
mis_list = mis_raw_in.split(', ')

mis_fin = []
for m in mis_list:
    mis_fin.append(int(m)-1)

nodes = {0, 1}
for n in range(0, size):
    nodes.add(n)

print("This is your entered MIS:")
print(mis_fin)


# this checks if there are too many nodes in the entered MIS:
count = 0
for i in mis_fin:
    print("node %d" % (i))
    if (i not in nodes):
        count = count + 1
        print("wrong MIS - node %d has already a neighbouring node, which is part of the MIS, thus this node can't be part of this MIS" % (i))
    else:
        for j in range(0, size):
            if (newmatrix[i, j] == 1.0):
                if (j in nodes):
                    print("remove edge (%d, %d)" % (i, j))
                    nodes.remove(j)
        print("removed node %d" % (i))
        nodes.remove(i)

if (count == 0):
    print("No faults found")
else:
    print("There has been %d nodes found which already have neighbouring nodes in the MIS. Thus the entered MIS is not correct." % (count))

# it should be checked as well, if the entered node set is an actual MIS or only a IS!
