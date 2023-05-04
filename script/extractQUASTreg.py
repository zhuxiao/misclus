import sys
with open(sys.argv[1], 'r') as f:
    tmp = f.read().split("Analyzing contigs...")
    s = tmp[1].split("\n\n")
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
                print(scaffold,"\t",seq[1],"\t",seq[2],sep="")
                mis=[]
