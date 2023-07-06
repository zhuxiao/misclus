# misClas #
Misassembly classification based on paired-end reads

----------
misClas is a tool that aims to classify the misassembly regions using paired-end reads data. It analyzes the signals from the given regions to determine whether the intervals include assembly errors.

## Prerequisites ##
misClas depends on the following libraries and tools:
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
HTSlib and mlpack need to be installed from source files.

## Building misClas ##

misClas can be built by typing:
```sh
$ git clone https://github.com/zhuxiao/misclas.git
$ cd misclas/
$ ./autogen.sh
```
Then, the program `misclas` will be generated in directory `misclas`.

Alternatively, misClas can also be built using the released package `misclas_${VERSION}.tar.gz` from `https://github.com/zhuxiao/misclas/releases`:
```sh
$ tar zxvf misclas_${VERSION}.tar.gz
$ cd misclas_${VERSION}/
$ ./autogen.sh
```
Then, the program `misclas` will be generated in directory `misclas_${VERSION}`.

It is recommended to create the symbolic link of the program `misclas` to the `PATH` directories on the machine.

## misClas command ##

Simply，misClas can be run by the command:
```sh
$ misclas -r referenceSequence.fa -b bamFile.bam -f regionFile
```
Then, the assembly error regions and normal regions will be output in `misassemblyRegion` file and `normalRegion` file, respectively.

The help information can be shown:
```sh
Program: misclas
Version: 0.4.2

Usage:  misclas [options] -r <ASM> -b <BAM> -f <REG>

Description:
    ASM     Assembled sequences file (required)
    BAM     Coordinate-sorted BAM file (required)
    REG     Regions to be analyzed: CHR|CHR:START-END.(required)

Inputs:
    -r PATH    Assembled sequences file
    -b PATH    Bam file
    -f PATH    Regions file
General:
    -ranNorRegCoef [0.2]
        set the ratio of normal random regions (generated by misClas randomly) to candidate regions
    -min-regsize INT [200]
        the candidate misassembly region will be extended to INT bp symmetrically 
        when it is shorter than INT bp. 
    -extend-regsize-fold FLOAT [2]
        auxiliary region is used to set baseline for an candidate misassembly region, 
        it is determined by extending the candidate misassembly region by 
        FLOAT*(the candidate misassembly region size) on both sides. 

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

## Misassembly region extraction scripts ##
The scripts in directory `script` can be used to extract regions used for misClas, the function of these scripts as below:
* `extractQuastReg.py`: extract misassembly regions from the results of QUAST;
* `filterOverlapReg.py`: filter out the overlapped regions in extracted region file;
* `generateRandReg.py`: randomly select certain number of regions from the extracted region file;

Furthermore, the examples for the misassembly detection results for QUAST, Misasm, Pilon, REAPR and misFinder can be downloaded at "https://github.com/zhuxiao/results_misclas/tree/origin/Ecoli_real/misassembly_detection", and the misassembly region extract commands for these tools are as below:

* QUAST
```sh
$ extractQUASTreg.py misassemblies_QUAST misReg_QUAST
```
Then, the extracted misassembly regions for QUAST will be saved in the `misReg_QUAST` file.

* Misasm
```sh
$ cat genome_Indel genome_Misjoin | awk '{split($1,a,":"); split(a[2],b,"-"); print a[1]"\t"b[1]"\t"b[2]}' > misReg_Misasm
```
Then, the extracted misassembly regions for Misasm will be saved in the `misReg_Misasm` file.

* Pilon
```sh
$ cat misassemblies_Pilon | grep "fix" | grep ":" | cut -d ":" -f 2,3 | awk '{split($1,a,":");split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > misReg_Pilon
```

Then, the extracted misassembly regions for Pilon will be saved in the `misReg_Pilon` file.

* REAPR
```sh
$ cat misassemblies_REAPR | grep "Error" | awk '{print $1"\t"$4"\t"$5}' > misReg_REAPR
```
Then, the extracted misassembly regions for REAPR will be saved in the `misReg_REAPR` file.

* misFinder
```sh
$ cat misassemblies_misFinder | awk '{print $1"\t"$3"\t"$4}' > misReg_misFinder
```
Then, the extracted misassembly regions for misFinder will be saved in the `misReg_misFinder` file.

## Extracted misassembly region file format ##
There are three fields in extracted misassembly region file: (1) name of scaffold/contig, (2) start position of the extracted misassembly region, (3) end position of the extracted misassembly region. And these fields are separated by tabulators.
```sh
#scaffold	startPos	endPos
scaffold_1	43698	44050
scaffold_3	14720	14886
scaffold_6	171522	171824
...
```

## Contact ##

If you have problems or some suggestions, please contact: [xzhu@ytu.edu.cn](xzhu@ytu.edu.cn) without hesitation. 

---- Enjoy !!! -----

