#!/usr/bin/env python3
import sys

def printHelp():
    print("""
Usage:  extractQUASTreg.py [options] <QUASTFILE> <OUTFILE>

Positional options:
    QUASTFILE
        The stdout file reported by QUAST (required).
    OUTFILE
        The output file name (required).

Options:
    -h,-help  Show this help message and exit
""")
    exit()

if(sys.argv.count('-h') or sys.argv.count('-help')):
    printHelp()

if(len(sys.argv)!=3):
    printHelp()



with open(sys.argv[1], 'r') as f:
    tmp = f.read().split("Analyzing contigs...")
    s = tmp[1].split("\n\n")

output = open(sys.argv[2],'w+')
output.write("#scaffold\tstartPos\tendPos")
num=0
for i in s:
    if (i.__contains__("Indel") or i.__contains__("local misassembly") or i.__contains__("Extensive misassembly")):
        lines=i.split("\n")
        num = num+1
        mis = []
        for j in range (len(lines)):
            if lines[j].__contains__("CONTIG"):
                scaffold = lines[j].split(":")[1].split()[0]
            if (lines[j].__contains__("Indel") or lines[j].__contains__("local misassembly") or lines[j].__contains__("Extensive misassembly")):
                align1 = lines[j-1].split("|")[1]
                pos1 = align1.split()
                if  lines[j+1].__contains__("Real Alignment"):
                    align2 = lines[j+1].split("|")[1]
                else:
                    align2 = lines[j+2].split("|")[1]
                pos2 = align2.split()
                mis.append(int(pos1[0]))
                mis.append(int(pos1[1]))
                mis.append(int(pos2[0]))
                mis.append(int(pos2[1]))
                seq = mis
                seq.sort()
                
                output.write("\n"+scaffold+"\t"+str(seq[1])+"\t"+str(seq[2]))
                mis=[]

output.close()