# Running Novoalign mapper

Starting with sequence files that have already been inspected with md5sum and fastqc and have been trimmed using trimmomatic (See other files for running trimmomatic)

Novoalign link tutorial: http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/basic-short-read-mapping/

1) Need to create directory for project (project_dir) and to house mapping outputs

ex. mkdir novoalign (project), and mkdir novo_dir (mapping outputs)

2) Find path to call Novoalign

ex. /usr/local/novoalign


3) Make a scripts directory to house scripts to run 

ex. mkdir novo_scripts

4) Novoindex reference

Need index dir: mkdir novo_index

Example Script:
```
#! /bin/bash

#Create variable for location of reference genome
ref_genome=/home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta

#Variables for project, novoalign, and output directory
project_dir=/home/paul/episodicData/novoalign
novoalign=/usr/local/novoalign
novo_index=${project_dir}/novo_index

#Index the reference
${novoalign}/novoindex ${novo_index}/dmel-all-chromosome-r5.57_2.nix  ${ref_genome}
```

5) Run Novoalign
Note: Compressed read files are not supported in unlicensed versions.

Need to unzip (in trimmomatic output dir)
```
gunzip *.gz
```

Novoalign Flags
-d == Full pathname of indexed reference sequence from novoindex
-f == Files containing the read sequences to be aligned
-o == Specifies output report format and options (SAM)

-i ###,## ==
Sets fragment orientation and approximate fragment length for proper pairs.
ex. -i 250 50  Defaults to paired end Illumina or Mate Pair ABI with 250bp insert and 50bp standard deviation
- Using 400,100 based on initial mapping with novoalign first run through
See Kofler *et al.* 2016 for some idea of different flags 

The script
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=/home/paul/episodicData/trim_dir

#Path to output directory
novo_dir=${project_dir}/novo_dir


files=(${trim_dir}/*_R1_PE.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq`
${novoalign}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 400,100 -o SAM > ${novo_dir}/${base}_novo.sam

done
```
Takes a long time: Only uses 1 thread (100% computer)
From novoalign reference manual: -c 99 Sets the number of threads to be used. On licensed versions it defaults 
to the number of CPUs as reported by sysinfo(). On free version the 
option is disabled

To avoid this problem, run scripts for each in parallel

1) make the script to make multiple scripts
Make dir for all output scripts:
./split_mappingScripts/
Run in Scripts dir ( not split_mappingScripts)
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=/home/paul/episodicData/trim_dir

#Path to output directory
novo_dir=${project_dir}/novo_dir

files=(${trim_dir}/*_R1_PE.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq`
echo "${novoalign}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 400,100 -o SAM > ${novo_dir}/${base}_novo.sam" > ./split_mappingScripts/${base}.sh

done
```
2) Create script to call all and run in parallel (us & which puts job in background
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Variable for each script location
map_scripts=${project_dir}/novo_scripts/split_mappingScripts

#Variable for location of whole script (novo_scripts)
scripts=${project_dir}/novo_scripts


files=(${map_scripts}/*.sh)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sh`
echo "${map_scripts}/${base}.sh &" > ${scripts}/novo_parallel_map.sh

done
```
Change permissions and run novo_parallel_map.sh to run all files at once (can change input parametes to run subsets on different days)

Rezip files in trim_dir (saves space)
```
gzip *.fastq
```

sam-bam files: need bam directory (mkdir novo_bam)
Flags:
- b == output = .bam, -S = input .sam, -q 20 = quality mapping score of 20 (standard throughout all experiments)
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_dir=${project_dir}/novo_dir

#Path to output directory
novo_bam=${project_dir}/novo_bam

files=(${novo_dir}/*.sam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sam`
samtools view -b -S -q 20 ${novo_dir}/${base}.sam | samtools sort -o ${novo_bam}/${base}.bam
done
```


Next Steps:
> Merge; merge the base multiply
    > no flags
    > mkdir novo_merge


```    
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_merge=${project_dir}/novo_merge


files=(${novo_bam}/*_L001_novo.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_novo.bam`
samtools merge ${novo_merge}/${base}_novo_merge.bam ${novo_bam}/${base}_L001_novo.bam ${novo_bam}/${base}_L002_novo.bam
done

```
  
> Picard Sort
    > java -Xmx2g -Djava.io.tmpdir=${tmp} -jar ${pic} SortSam I= ${merged}/${base}.bam O= ${sort_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${tmp}

> Remove Duplicates

> More QC

> Create mpileup

> Create .sync file
