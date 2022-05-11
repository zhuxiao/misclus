# misEval #
Misassembly Evaluator based on paired-end reads

----------
misEval is a tool that aims to evaluate the misassembly regions using paired-end reads data. It analyzes the signals from the given regions to determine whether the intervals include assembly errors.

## Prerequisites ##
misEval depends on the following libraries and tools:
1. [HTSlib(V1.9 or later)](https://github.com/samtools/htslib)
2. [mlpack(v3.4.2 or later)](https://github.com/mlpack/mlpack)
3. g++ (v4.7 or later which supports c++11).
4. python3
5. autotools

## Quick install ##
* Debian / Ubuntu 
```
$ sudo apt update  # Ensure the package list is up to date
$ sudo apt install mlpack-bin g++ python3 autoconf automake libtool
```
HTSlib needs to be installed from source files.

* RedHat / CentOS
```
$ sudo yum install g++ python3 autoconf automake libtool
```
HTSlib and mplack need to be installed from source files.

## Building misEval ##

The miseval can be builded by typing:
```sh
$ git clone https://github.com/zhuxiao/miseval.git
$ cd miseval/
$ ./autogen.sh
```

## misEval command ##

Simply，misEval can be run by the command:
```sh
$ miseval.py -r referenceSequence.fa -b bamFile.bam -f regionFile
```
Then, the assembly error regions and normal regions will be output in `misassemblyRegion` file and `normalRegion` file, respectively.

The help information can be shown:
```sh
Program: miseval
Version: 0.4.0

Usage:  miseval.py [options] -r <REF> -b <BAM> -r <REG>

Description:
    REF     Reference file (required)
    BAM     Coordinate-sorted BAM file (required)
    REG     Regions to be analyzed: CHR|CHR:START-END.(required)

Inputs:
    -r PATH    Reference file
    -b PATH    Bam file
    -f PATH    Regions file
General:
    -min-regsize INT [200]
        the evaluated region will be extended to INT bp symmetrically 
        when it is shorter than INT bp. 
    -extend-regsize-fold FLOAT [2]
        auxiliary region is used to set baseline for an evaluated region, 
        it is determined by extending the evaluated region by 
        FLOAT*(the evaluated region size) on both sides. 

Mode: 
    -cov on/off coverage clustering [on]
        --min-cov-fold FLOAT [0.5]
            the coverage is anomalous when it is lower than FLOAT*meancov,
            where meancov is the average coverage of auxiliary region.
        --max-cov-fold FOLAT [2]
            the coverage is anomalous when it is higher than FLOAT*meancov.
    -indels on/off indels clustering [on]
        --min-locs-ratio FLOAT [0.2]
            mask the ultra-low coverage locations whose coverage is lower 
            than FLOAT*meancov. When indels appear at ultra-low coverage location,
            the performance on indels is unreliable.
    -abstrand on/off anomalous orientation reads clustering [0.2]
        --min-abstrand-ratio FLOAT [0.3]
            search the anomalous orientation clumps in which the anomalous 
            orientation reads take a higher ratio than FLOAT at each location.
    -abisize on/off anomalous insertsize reads clustering [on]
        --min-abisize-ratio FLOAT [0.3]
            search the anomalous insertsize clumps in which the anomalous 
            insertsize reads take a higher ratio than FLOAT at each location.
        --isize-sdev-fold FLOAT [3]
            the anomalous insertsize is short than (average_isize)-FLOAT*sdev bp 
            or longer than (average_isize)+FLOAT*sdev bp, otherwise it is anomalous.
    -abmate on/off chromosomal rearrangement reads clustering [on]
        --min-abmate-ratio FLOAT [0.3]
            search the chromosomal rearrangement clumps in which the chromosomal 
            rearrangement reads take a higher ratio than FLOAT at each location.
   
Options:
    -o PATH    Outdir
    -h         Show this help message and exit
```

## Contact ##

If you have problems or some suggestions, please contact: [zhuxiao_hit@yeah.net](zhuxiao_hit@yeah.net) without hesitation. 

---- Enjoy !!! -----

