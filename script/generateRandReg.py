#!/usr/bin/env python3
import sys
import random

def printHelp():
    print("""
Usage:  generateRandReg.py [options] <REG> <OUTFILE> 

Input:
    REG
        Region file: CHR|CHR:START-END (required).
    OUTFILE
        The output file name (required).
        
        
Options:
    -n        The number of randomly selected regions to output [100].
    -h,-help  Show this help message and exit
""")
    exit()
    
if(sys.argv.count('-h') or sys.argv.count('-help') or len(sys.argv)<=1):
    printHelp()

ranNum = 200
if(sys.argv.count("-n")):
    with open(sys.argv[3]) as f:
        reg = f.readlines()
        ranNum = int(sys.argv[2])
        output = open(sys.argv[4],'w+')
else: 
    with open(sys.argv[1], 'r') as f:
        next(f)
        reg = f.readlines()
        output = open(sys.argv[2],'w+')

output.write("scaffold\tstartPos\tendPos\n")

num = 0
tmp = reg[0].split("\t")
vec = []
for i in reg:
    tmplist = i.split()
    if len(tmplist) != 4:
        continue
    if tmp[0] == i.split("\t")[0]:
        vec.append(i)
    else:
        tmp = i.split("\t")
        output.write(random.choice(vec))
        vec=[]
        vec.append(i)
        num = num +1
        if num>=ranNum:
            break

output.close()
