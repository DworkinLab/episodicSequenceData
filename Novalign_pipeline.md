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

Flags
-d == Full pathname of indexed reference sequence from novoindex
-f == Files containing the read sequences to be aligned
-o == Specifies output report format and options (SAM)

-i ###,## ==
Sets fragment orientation and approximate fragment length for proper pairs.
ex. -i 250 50  Defaults to paired end Illumina or Mate Pair ABI with 250bp insert and 50bp standard deviation

```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/Novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=${project_dir}/../trim_dir

#Path to output directory
novo_dir=${project_dir}/novo_dir


files=(${trim_dir}/*_R1_PE.fastq.gz)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
${novoalign}/novoalign -d ${novo_index} -f  ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz -i 200,50 -o SAM > ${novo_dir}/${base}_novo.sam > log.txt

done
```
