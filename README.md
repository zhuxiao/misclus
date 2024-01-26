# misClus #
Misassembly clustering based on paired-end reads

----------
misClus is a tool that aims to cluster the misassembly regions using paired-end reads data. It analyzes the signals from the given regions to determine whether the intervals include assembly errors.

For more detailed misclus experimental information, please refer to [misclus-experiments](https://github.com/zhuxiao/misclus-experiments)

## Prerequisites ##
misClus depends on the following libraries and tools:
1. [HTSlib(V1.9 or later)](https://github.com/samtools/htslib)
2. [mlpack(v3.4.2 or later)](https://github.com/mlpack/mlpack)
3. g++ (v4.7 or later which supports c++11).
4. python3
5. autotools

## Quick install ##
* Debian / Ubuntu 
```bash
$ sudo apt update  # Ensure the package list is up to date
$ sudo apt install mlpack-bin g++ python3 autoconf automake libtool
```
HTSlib needs to be installed from source files.

* RedHat / CentOS
```bash
$ sudo yum install g++ python3 autoconf automake libtool
```
HTSlib and mlpack need to be installed from source files.

## Building misClus ##

misclus can be built using the released package `misclus_0.5.0.tar.gz` from `https://github.com/zhuxiao/misclus/releases`:

```sh
$ wget -c https://github.com/zhuxiao/misclus/releases/download/0.5.0/misclus_0.5.0.tar.xz
$ tar zxvf misclus_0.5.0.tar.xz
$ cd misclus_0.5.0/
$ ./autogen.sh
```

Then, the program `misclus` will be generated in directory `bin`.

Alternatively, misClus can also  be built by typing:

```sh
$ git clone https://github.com/zhuxiao/misclus.git
$ cd misclus/
$ ./autogen.sh
```
Then, the program `misclus` will be generated in directory `bin`.

In addition,misclus can also be built by docker.

```bash
$ docker pull xzhu2/misclus:0.5.0
$ docker run -it --name your_contaniner_name misclus_iamges:0.5.0
$ cd /opt/misclus/
$ ./autogen.sh
```

Then, the program `misclus` will be generated in directory `bin`

It is recommended to create the symbolic link of the program `misclus` to the `PATH` directories on the machine.

## misClus command ##

Simply，misClus can be run by the command:
```sh
$ misclus regionFile referenceSequence.fa  bamFile.bam  
```
Then, the assembly error regions and normal regions will be output in `misassemblyRegion` file and `normalRegion` file, respectively.

The help information can be shown:
```sh
Program: misclus
Version: 0.5.0 (using htslib 1.14)

Usage:  misclus [options] <REG_FILE> <REF_FILE> <BAM_FILE>

Description:
   REG_FILE   Region file that need to be processed (required).
              Each region in REG_FILE should be specified in the format
              per line: CHR|CHR:START-END.
   REF_FILE   Reference file (required)
   BAM_FILE   Coordinate-sorted BAM file (required)

Options: 
   -t INT     Number of threads [0]. 0 for the maximal number
              of threads in machine
   -o DIR     Output directory [output]
   --norm-reg-percent DOUBLE
              The percentage of total regions to generate random normal regions[0.2]
              These randomly generated normal regions are used for comparison with
              candidate misassembly regions during clustering.
   --min-regsize  INT
              The minimal region size for feature extraction. The region will be extended to
              INT bp on both sides if it is shorter than INT bp.[200]
   --extend-regsize-fold DOUBLE
              Auxiliary region is used for an candidate misassembly region.
              The auxiliary region is determined by extending the candidate region by
              DOUBLE fold of the size of candidate region on both sides.[2]
   --cov-off  disable the clustering on coverage feature
      --min-cov-fold DOUBLE
              Coverage is anomalous when it is lower than DOUBLE*meancov for an auxiliary region,
              where meancov is the average coverage of the auxiliary region.[0.5]
              This option takes effect only if '--cov-off' option is not specified.
      --max-cov-fold DOUBLE
              Coverage is anomalous when it is higher than DOUBLE*meancov for an auxiliary region.[2]
              This option takes effect only if '--cov-off' option is not specified.
   --indel-off  disable the clustering on indel feature
      --min-locs-ratio DOUBLE
              Mask the ultra-low coverage locations whose coverage is lower
              than DOUBLE*meancov. When indels appear at ultra-low coverage location,
              the performance on indels is unreliable.[0.2]
              This option takes effect only if '--indel-off' option is not specified.
   --ab-strand-off  disable the clustering on anomalous orientation feature
      --min-abstrand-ratio DOUBLE
              Search the anomalous orientation clumps in which the anomalous
              orientation reads take a higher ratio than DOUBLE at each location.[0.3]
              This option takes effect only if '--ab-strand-off' option is not specified.
   --ab-isize-off  disable the clustering on anomalous insert size feature
      --min-abisize-ratio DOUBLE
              Search the anomalous insert size clumps in which the anomalous
              insert size reads take a higher ratio than DOUBLE at each location.[0.3]
              This option takes effect only if '--ab-isize-off' option is not specified.
      --isize-sdev-fold DOUBLE
              Normal insert size is usually in the range of 
              [isize-DOUBLE*sdev, isize+DOUBLE*sdev] bp, otherwise it is anomalous.[3]
              This option takes effect only if '--ab-isize-off' option is not specified.
   --ab-mate-off  disable the clustering on abnormal mate read pair feature
      --min-abmate-ratio DOUBLE
              Search the mate read pair clumps in which the abnormal mate read pairs
              take a higher ratio than DOUBLE at each location.[0.5]
              This option takes effect only if '--ab-mate-off' option is not specified.
   -v         show version information
   -h         show this help message and exit
```

## Misassembly region extraction scripts ##
The scripts in directory `script` can be used to extract regions used for misClus, the function of these scripts as below:
* `extractQuastReg.py`: extract misassembly regions from the results of QUAST;
* `filterOverlapReg.py`: filter out the overlapped regions in extracted region file;
* `generateRandReg.py`: randomly select certain number of regions from the extracted region file;

Furthermore, the examples for the misassembly detection results for QUAST, Misasm, Pilon, REAPR and misFinder can be downloaded at "https://github.com/zhuxiao/misclus-experiments/tree/origin/Ecoli_real/misassembly_detection", and the misassembly region extraction commands for these tools are as below:

* QUAST
```sh
$ wget -c https://github.com/zhuxiao/misclus-experiments/blob/origin/Ecoli_real/misassembly_detection/QUAST/misassemblies_QUAST
$ extractQUASTreg.py misassemblies_QUAST misReg_QUAST
```
Then, the extracted misassembly regions for QUAST will be saved in the `misReg_QUAST` file.

* Misasm
```sh
$ wget -c https://github.com/zhuxiao/misclus-experiments/blob/origin/Ecoli_real/misassembly_detection/Misasm/genome_Indel
$ wget -c https://github.com/zhuxiao/misclus-experiments/blob/origin/Ecoli_real/misassembly_detection/Misasm/genome_Misjoin
$ cat genome_Indel genome_Misjoin | awk '{split($1,a,":"); split(a[2],b,"-"); print a[1]"\t"b[1]"\t"b[2]}' > misReg_Misasm
```
Then, the extracted misassembly regions for Misasm will be saved in the `misReg_Misasm` file.

* Pilon
```sh
$ wget -c https://github.com/zhuxiao/misclus-experiments/blob/origin/Ecoli_real/misassembly_detection/Pilon/misassemblies_Pilon
$ cat misassemblies_Pilon | grep "fix" | grep ":" | cut -d ":" -f 2,3 | awk '{split($1,a,":");split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > misReg_Pilon
```

Then, the extracted misassembly regions for Pilon will be saved in the `misReg_Pilon` file.

* REAPR
```sh
$ wget -c https://github.com/zhuxiao/misclus-experiments/blob/origin/Ecoli_real/misassembly_detection/REAPR/misassemblies_REAPR
$ cat misassemblies_REAPR | grep "Error" | awk '{print $1"\t"$4"\t"$5}' > misReg_REAPR
```
Then, the extracted misassembly regions for REAPR will be saved in the `misReg_REAPR` file.

* misFinder
```sh
$ wget -c https://github.com/zhuxiao/misclus-experiments/blob/origin/Ecoli_real/misassembly_detection/REAPR/misassemblies_REAPR
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