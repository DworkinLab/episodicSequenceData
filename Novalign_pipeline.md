# Running Novoalign mapper

Starting with sequence files that have already been inspected with md5sum and fastqc and have been trimmed using trimmomatic (See other files for running trimmomatic)

Novoalign link tutorial: http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/basic-short-read-mapping/

1) Need to create directory for project (project_dir) and to house mapping outputs

ex. mkdir novoalign (project), and mkdir novo_dir (mapping outputs)

2) Find path to call Novoalign

ex. /usr/local/novoalign

3) Novoindex referenece





Modify example script

```
novoalign -d hg18.nix -f SRR040810_1.fastq.gz SRR040810_2.fastq.gz -i 200,50 -o SAM > alignments.sam > log.txt
```
```
#! /bin/bash

project_dir=/home/paul/episodicData/Novoalign
index_dir=${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir}/trim_dir
novoalign=/usr/local/novoalign/novoalign
novo_dir=${project_dir}/novo_dir

```
