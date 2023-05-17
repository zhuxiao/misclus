import sys
import random
with open(sys.argv[1], 'r') as f:
    reg = f.readlines()

print("scaffold\tstartPos\tendPos")
num = 0
tmp = reg[0].split("\t")
vec = []
for i in reg:
    if tmp[0] == i.split("\t")[0]:
        vec.append(i)
    else:
        tmp = i.split("\t")
        print(random.choice(vec),end="")
        vec=[]
        vec.append(i)
        num = num +1
        if num>=100:
            break

