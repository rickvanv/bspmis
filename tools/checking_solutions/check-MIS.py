from io import StringIO
from scipy.io import mmread

input_file = "ca-GrQc.mtx"
matrix = mmread(input_file)
newmatrix = matrix.tocsr()
size = matrix.shape[0]

mis_raw_file = open("../../output.txt", "r")
mis_raw_in = mis_raw_file.read()[:-2]
mis_list = mis_raw_in.split(', ')

mis_fin = []
for m in mis_list:
    mis_fin.append(int(m)-1)

nodes = {0,1}
for n in range(0, size):
    nodes.add(n)

print("Checking if MIS found of %s is an IS" %(input_file))
count = 0
for i in mis_fin:
    if (i not in nodes):
        count = count + 1
        print("wrong MIS - node %d has already a neighbouring node, which is part of the MIS, thus this node can't be part of this MIS" % (i))
        break
    else:
        for j in range(0, size):
            if (newmatrix[i, j] == 1.0):
                if (j in nodes):
                    nodes.remove(j)
        if (i not in nodes):
            count = count + 1
            print("wrong MIS - node %d has already a neighbouring node, which is part of the MIS, thus this node can't be part of this MIS" % (i))
            break
        else:
            nodes.remove(i)

if (count == 0):
    print("No faults found, hence MIS found is an IS")
else:
    print("There has been %d nodes found which already have neighbouring nodes in the MIS. Thus the entered MIS is not correct." % (count))
    exit("Not even an IS")

print("Checking if MIS found is an MIS")
notamis = False
for i in range(0,size):
    if i in mis_fin:
        continue
    else:
        nonzeros = [k for k in range(0,size) if newmatrix[i, k] == 1.0]
        inmis = False
        for j in nonzeros:
            if j in mis_fin:
                inmis = True
                break
        if not inmis:
            notamis = True
            break

if notamis:
    print("Found IS is not a MIS")
else:
    print("Found IS is a MIS!!")





# it should be checked as well, if the entered node set is an actual MIS or only a IS!
