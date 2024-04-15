#!/usr/bin/env python3
import sys

def printHelp():
    print("""
Usage:  filterOverlapReg.py [options] <REG> <FAI> <OUTFILE>

Positional options:
    REG
        Region file: CHR|CHR:START-END (required).
    FAI
        fai file: the index of scaffold (required).
    OUTFILE
        The output file name (required).

Options:       
    -h,-help  Show this help message and exit
""")
    exit()

if(sys.argv.count('-h') or sys.argv.count('-help')):
    printHelp()

if(len(sys.argv)!=4):
    printHelp()

with open(sys.argv[1], 'r') as finreg:
    reg = finreg.read().splitlines()
with open(sys.argv[2], 'r') as fai:
    scaffold = fai.read().splitlines()

output = open(sys.argv[3],'w+')
reglist = []
scaffoldend = []
output.write("scaffold\tstartPos\tendPos")
for i in reg:
    tmplist = i.split()
    if len(tmplist) == 3:
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
        output.write("\n"+i[0]+"\t"+i[1]+"\t"+i[2])
output.write("\n")
output.close()
