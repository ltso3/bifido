import fileinput
import re

with open('protein_mapping.txt', 'r') as file:
    data = file.read().replace('\n', '')

mapping = eval(data)

print(len(mapping))

counter = 0
for line in fileinput.input("protein_hits_blon.faa", inplace=True):
    if re.search("^>", line):
        acc, blon = mapping[counter]
        if re.search("pdb", line):
            tmp = acc.split("_")[0]
            print("-".join(line.replace("pdb|"+tmp+"|", blon+" "+acc+" ").replace(',', '').strip().split()))
        else:
            print("-".join(line.replace(acc, blon+" "+acc).replace(',', '').strip().split()))
        counter += 1
    else:
        print(line.strip())
