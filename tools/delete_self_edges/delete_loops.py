import re

#TODO: adapt it to change .mtx file instead of .mtx-Px file
name = input(
    "Please enter the filename of your graph to delete it's selfedges: ")

file = open(name)

lines = file.readlines()
count = 0
reached_firstline = False

firstline_split = []
firstline_idx = 0

for idx, line in enumerate(lines):
    if line[0] != "%":
        if not reached_firstline:
            reached_firstline = True
            firstline_idx = idx
            firstline_split = re.split(" ", line)
            print(firstline_split)
            continue
        node = re.split(" ", line)
        if (len(node) > 1) and (len(node) < 4):
            if int(node[0]) == int(node[1]):
                print("This node has a self edge: " +
                      node[0] + " - REMOVED!")
                lines.remove(line)
                count = count + 1

if not reached_firstline:
    print("Error: File is empty")

firstline_split[2] = str(int(firstline_split[2]) - count) + "\n"
print(' '.join(firstline_split))
lines[firstline_idx] = ' '.join(firstline_split)


print("%d self edges have been removed. New file without selfedges is created: " %
      count + name)

new_file = open(name, 'w')
new_file.writelines(lines)
