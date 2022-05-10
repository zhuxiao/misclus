#!/usr/bin/env python3
from curses.ascii import isdigit
import os
import sys
import time
def printHelp():
    print("""Program: miseval v0.3.0

    Usage:  miseval.py [options] -r <REF> -b <BAM> -r <REG>

    Description:
        REF     Reference file (required)
        BAM     Coordinate-sorted BAM file (required)
        REG     Regions to be analyzed: CHR|CHR:START-END.(required)
    Inputs:
        -r PATH        reference file
        -b PATH        bam file
        -f PATH        regions file
    General:
        -min-regsize INT    [200]
            the evaluated region will be extended to INT bp symmetrically 
            when it is shorter than INT bp. 
        -extend-regsize-fold    FLOAT   [2]
            auxiliary region is used to set baseline for an evaluated region, 
            it is determined by extending the evaluated region by 
            FLOAT*(the evaluated region size) on both sides. 

    Mode: 
        -cov    on/off    coverage clustering   [on]
            --min-cov-fold  FLOAT   [0.5]
                the coverage is anomalous when it is lower than FLOAT*meancov,
                 where meancov is the average coverage of auxiliary region.
            --max-cov-fold  FOLAT   [2]
                the coverage is anomalous when it is higher than FLOAT*meancov.
        -indels on/off   indels clustering      [on]
            --min-locs-ratio    FLOAT   [0.2]
                mask the ultra-low coverage locations whose coverage is lower 
                than FLOAT*meancov. when indels appear at ultra-low coverage,
                the performance on indels is unreliable.
        -abstrand   on/off   anomalous orientation reads clustering     [0.2]
            --min-abstrand-ratio    FLOAT   [0.3]
                search the anomalous orientation clumps in which the anomalous 
                orientation reads take a higher ratio than FLOAT at each location.
        -abisize    on/off    anomalous insertsize reads clustering     [on]
            --min-abisize-ratio     FLOAT   [0.3]
                search the anomalous insertsize clumps in which the anomalous 
                insertsize reads take a higher ratio than FLOAT at each location.
            --isize-sdev-fold FLOAT     [3]
                the anomalous insertsize is short than (the average insertsize)-FLOAT*sdev bp 
                or longer than  (the average insertsize)+FLOAT*sdev bp, otherwise it is anomalous.
        -abmate     on/off   chromosomal rearrangement reads clustering     [on]
            --min-abmate-ratio  FLOAT   [0.3]
                search the chromosomal rearrangement clumps in which the chromosomal 
                rearrangement reads take a higher ratio than FLOAT at each location.
   
    Options:
        -o PATH             outdir
        -h                  show this help message and exit
    """)
    exit()

def clusterall():
    files = os.listdir("./")
    if 'result' in files:
        files.remove('result')
    else:
        os.mkdir('result')
    for i in files:
        cmd = 'mlpack_kmeans -c 2 -i {0} -o ./result/{1}'.format(i,i)
        os.system(cmd)

def format(i):
    flag = 0
    f = open(i)
    fvec = []
    while(1):
        rline = f.readline()
        if (rline==''):
            break
        itmp = rline.strip().split(',')
        itmp = list(map(float, itmp))
        # print(itmp)
        fvec.append(itmp)
    for j in fvec:
        if (j[0] > fvec[0][0] and j[1] < fvec[0][1]):
            flag = 1
            break
    if (flag==1):
        for k in fvec:
            k[1] = 1 - k[1]
    return fvec

start = time.time()
minRegsize = '200'
exRegFold = '2'
cov = 'on'
minCovFold='0.5'
maxCovFold='2'
indel = 'on'
minLocRatio = '0.2'
abstrand = 'on'
minStrandRatio = '0.3'
abIsize = 'on'
minisizeRatio = '0.3'
IsizeSdevFold='3'
abMate = 'on'
minMateRatio = '0.3'
outdir = "misEval_out"

try:
    if(sys.argv.count('-r') and sys.argv.count('-b') and sys.argv.count('-f')):
        if not (os.path.exists(sys.argv[sys.argv.index('-r')+1])):
            print("reference file is not valid")
            printHelp()
        else:
            reFile = sys.argv[sys.argv.index('-r')+1]
        if not (os.path.exists(sys.argv[sys.argv.index('-b')+1])):
            print("bam file is not valid")
            printHelp()   
        else:
            bamFile =  sys.argv[sys.argv.index('-b')+1]   
        if not (os.path.exists(sys.argv[sys.argv.index('-f')+1])):
            print("regions file is not valid")
            printHelp()
        else:
            regFile = sys.argv[sys.argv.index('-f')+1]
    else:
        printHelp()

    if(sys.argv.count('-h')):
        printHelp()

#region
    if(sys.argv.count('-min-regsize')):
        minRegsize = int(sys.argv[sys.argv.index('-min-regsize')+1])
    if(sys.argv.count('-extend-regsize-fold')):
        exRegFold = float(sys.argv[sys.argv.index('-extend-regsize-fold')+1])    

#coverage
    if(sys.argv.count('-cov')):
        cov = sys.argv[sys.argv.index('-cov')+1]

    if(sys.argv.count('--min-cov-fold')):
        minCovFold = float(sys.argv[sys.argv.index('--min-cov-fold')+1])

    if(sys.argv.count('--max-cov-fold')):
        maxCovFold = float(sys.argv[sys.argv.index('--max-cov-fold')+1])
#indels
    if(sys.argv.count('-indels')):
        indel = sys.argv[sys.argv.index('-indels')+1]

    if(sys.argv.count('--min-locs-ratio')):
        minLocRatio = float(sys.argv[sys.argv.index('--min-locs-ratio')+1])

#abstrand
    if(sys.argv.count('-abstrand')):
        abstrand = sys.argv[sys.argv.index('-abstrand')+1]
    if(sys.argv.count('--min-abstrand-ratio')):
        minStrandRatio = float(sys.argv[sys.argv.index('--min-abstrand-ratio')+1])

#abisize
    if(sys.argv.count('-abisize')):
        abIsize = sys.argv[sys.argv.index('-abisize')+1]
    if(sys.argv.count('--min-abstrand-ratio')):
        minisizeRatio = float(sys.argv[sys.argv.index('--min-abisize-ratio')+1])
    if(sys.argv.count('--isize-sdev-fold')):
        IsizeSdevFold = float(sys.argv[sys.argv.index('--isize-sdev-fold')+1])

#abmate
    if(sys.argv.count('-abmate')):
        abMate = sys.argv[sys.argv.index('-abmate')+1]
    if(sys.argv.count('--min-abmate-ratio')):
        minMateRatio = float(sys.argv[sys.argv.index('--min-abmate-ratio')+1])

    if(sys.argv.count('-o')):
        outdir = sys.argv[sys.argv.index('-o')+1]
    cmd = 'miseval {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17}'.format(reFile, bamFile, regFile, minRegsize, exRegFold, cov, minCovFold, maxCovFold, indel, minLocRatio, abstrand, minStrandRatio, abIsize, minisizeRatio, IsizeSdevFold, abMate, minMateRatio, outdir)
    print(cmd)
    os.system(cmd)
    os.chdir(outdir)
    clusterall()
    os.chdir('result')
#file analysis
    with open('../../{}'.format(regFile), 'r') as reg:
        regions = reg.read().splitlines()
    del(regions[0])
    f = open('clusterResult',mode='w')
    fout = open('regionClassification', mode='w')
    # num=0
    titleline = "#scaffold\tstartPos\tendPos"
    if (cov=='on'):
        covCluster=format('cov.csv')
        titleline += "\tcov"
        # num=len(covCluster)
    if (indel=='on'):
        indelCluster=format('indel.csv')
        titleline += "\tindel"
        # num=len(indelCluster)
    if (abstrand=='on'):
        abstrandCluster=format('strand.csv')
        titleline += "\tabstrand"
        # num=len(abstrandCluster)    
    if (abIsize=='on'):
        abIsizeCluster=format('insert.csv')
        titleline += "\tabisize"
        # num=len(abIsizeCluster)
    if (abMate=='on'):
        abmMateCluster=format('mate.csv')
        titleline += "\tabmate"
        # num=len(abmMateCluster)
# title+mark
    if (cov=='on'):
        titleline += "\tcovMark"
    if (indel=='on'):
        titleline += "\tindelMark"
    if (abstrand=='on'):
        titleline += "\tabstrandMark"
    if (abIsize=='on'):
        titleline += "\tabisizeMark"
    if (abMate=='on'):
        titleline += "\tabmateMark"
    titleline += "\tregionmark"
    # print(titleline)
    fout.write(titleline)
    fout.write('\n')
    last=[]
    last_1=0
    for i in range (len(regions)):
        linei = regions[i]
        if (cov=='on'):
            last_1 += covCluster[i][1]
            linei = linei + '\t' + str(covCluster[i][0])
        if (indel=='on'):
            last_1 += indelCluster[i][1]
            linei = linei + '\t' + str(indelCluster[i][0])
        if (abstrand=='on'):
            last_1 += abstrandCluster[i][1]
            linei = linei + '\t' + str(abstrandCluster[i][0])
        if (abIsize=='on'):
            last_1 += abIsizeCluster[i][1]
            linei = linei + '\t' + str(abIsizeCluster[i][0])
        if (abMate=='on'):
            last_1 += abmMateCluster[i][1]
            linei = linei + '\t' + str(abmMateCluster[i][0])

        if (cov=='on'):
            if not (covCluster[i][1]==0):
                linei = linei + '\t' + 'P'
            else:
                linei = linei + '\t' + 'N'
        if (indel=='on'):
            if not (indelCluster[i][1]==0):
                linei = linei + '\t' + 'P'
            else:
                linei = linei + '\t' + 'N'
        if (abstrand=='on'):
            if not (abstrandCluster[i][1]==0):
                linei = linei + '\t' + 'P'
            else:
                linei = linei + '\t' + 'N'
        if (abIsize=='on'):
            if not (abIsizeCluster[i][1]==0):
                linei = linei + '\t' + 'P'
            else:
                linei = linei + '\t' + 'N'
        if (abMate=='on'):
            if not (abmMateCluster[i][1]==0):
                linei = linei + '\t' + 'P'
            else:
                linei = linei + '\t' + 'N'
        if not (last_1==0):
            linei = linei + '\t' + 'P'
        else:
            linei = linei + '\t' + 'N' 
        # print(linei)      
        fout.write(linei) 
        fout.write('\n')
        last.append(last_1)
        last_1=0
    fout.close()
    for i in range (len(last)):
        # print(i+1,':',last[i],sep='')
        outline = str(i+1) + ':' + str(int(last[i])) + '\n'
        f.writelines(outline)
    f.close()
    with open('regionClassification', 'r') as fsplit:
        regSplit = fsplit.read().splitlines()
    regVec = []
    for i in range (len(regSplit)-1):
        regVec.append(regSplit[i+1].split())
    regPostive = []
    regNegative = []
    for i in range (len(regVec)):
        if (regVec[i][-1]=="P"):
            regPostive.append(i+1)
        else:
            regNegative.append(i+1)
    fpositive = open("misassemblyRegion",mode='w')
    fpositive.write('#scaffold\tstartPos\tendPos\tmark')
    fpositive.write('\n')
    for i in regPostive:
        fpositive.write(regions[i-1])
        fpositive.write("\tmisassembly")
        fpositive.write('\n')
    fpositive.close()
    fnegative = open("normalRegion",mode='w')
    fnegative.write("#scaffold\tstartPos\tendPos\tmark")
    fnegative.write('\n')
    for i in regNegative:
        fnegative.write(regions[i-1])
        fnegative.write("\tnormal")
        fnegative.write('\n')
    fnegative.close()  
    if (cov=='on'):
        os.remove("cov.csv")
    if (indel=='on'):
        os.remove("indel.csv")
    if (abstrand=='on'):
        os.remove("strand.csv")
    if (abIsize=='on'):
        os.remove("insert.csv")
    if (abMate=='on'):
        os.remove("mate.csv")
    os.remove("clusterResult")
    print("""
############# output files #############
The indicators of regions are saved in {0}/result/regionClassification.
The misassembly regions are saved in {1}/result/misassemblyRegion.
The normal regions are saved in {2}/result/normalRegion.""".format(outdir, outdir, outdir))
    end = time.time()
    print('\nRunning time: %s seconds'%(end-start))
except ValueError:
    print("invalid value")
    printHelp()