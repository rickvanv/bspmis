import re

name = input(
    "Please enter the filename of your graph (without '.mtx-P4') to delete it's selfedges: ")

file = open(name + ".mtx-P4")
# file = open("ca-GrQc.mtx-P4")

lines = file.readlines()
count = 0

for line in lines:
    if (line[0] != "%"):
        node = re.split(" ", line)
        if ((len(node) > 1) and (len(node) < 4)):
            if (int(node[0]) == int(node[1])):
                print("This node has a self edge: " +
                      node[0] + " - REMOVED!")
                lines.remove(line)
                count = count + 1
name = name + "-without-loops.mtx-P4"
print("%d self edges have been removed. New file withouth selfedges is created: " %
      (count) + name)

new_file = open(name, 'w')
new_file.writelines(lines)
