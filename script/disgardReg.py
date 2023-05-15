#!/usr/bin/env python3
import sys
with open(sys.argv[1], 'r') as finreg:
    reg = finreg.read().splitlines()
with open(sys.argv[2], 'r') as fai:
    scaffold = fai.read().splitlines()
reglist = []
scaffoldend = []
for i in reg:
    tmplist = i.split()
    int(tmplist[1])
    int(tmplist[2])
    reglist.append(tmplist)
for i in scaffold:
    tmplist = i.split()
    int(tmplist[1])
    scaffoldend.append(tmplist)
for i in reglist:
    for j in scaffoldend:
        if (j[0]==i[0]):
            tmp = j
            break
    if(int(i[1]) > 1000 and int(i[2]) < int(j[1])-1000):
        print(i[0],":",i[1],'-',i[2],sep='')