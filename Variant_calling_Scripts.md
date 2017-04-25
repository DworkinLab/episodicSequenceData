# Running Variant Calling packages: scripts

CRISP>
Varscan>
SnpEff>
LoFreq>
Snape>

### CRISP

Requirements:
-- fasta reference genome (with .fai indexed)
-- input .bam files (better if used with GATK indel realigner)
-- CRISP available (download or on ssh already)

FLAGS:
```
--bams/bam: bams == textfile with list of bam file paths (one for each pool), bam == bam file for one pool, specify filename for each pool using --bam pool1.bam --bam pool2.bam .... --bam pooln.bam

--ref: Indexed Reference Sequence file (fasta)

--poolsize: poolsize (number of haploid genomes in each pool), for diploid genomes: 2 x # individuals [120]

--perms: maximum number of permutations for calculating contingency table p-value, default 20,000

--filterreads: filter reads with excessive number of mismatches/gaps compared to the reference sequence, Default==1, 0 to disable

--regions: region(s) in which variants will be called, e.g chr1:654432-763332. BAM files should be indexed for this option with the index file pooln.bam.bai 

--qvoffset: quality value offset, 33 for Sanger format, 64 for Illumina 1.3+ format

--mbq: minimum base quality to consider a base for variant calling, default 10

--mmq: minimum read mapping quality to consider a read for variant calling, default 20

--minc: minimum number of reads with alternate allele required for calling a variant, default 4

--ctpval: threshold on the contingency table p-value for calling position as variant (specified as log10), default is -3.5, increase this threshold as number of pools increases

--qvpval: threshold on the quality values based p-value for calling position as variant (specified as log10), default is -5

--bed: bedfile with intervals in which variants will be called

--VCF: VCF file to which the variant calls will be output 
```

The Code:

```
#! /bin/bash

#Variable for project name (file name)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to CRISP
crisp=/home/paul/CRISP-122713/CRISP

#Variable for reference genome (non-zipped) ## .faidx indexed with samtools? .fai present
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

#Path to .bam files from GATK
input=${project_dir}/novo_GATK

#Output
output=${project_dir}/novo_crisp


files=(${input}/*_realigned.bam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _realigned.bam`
echo "${input}/${base}_realigned.bam" >> ${input}/${project_name}_BAMlist.txt

done

${crisp} --bams ${input}/${project_name}_BAMlist.txt \
			--ref ${ref_genome} \
 				--poolsize 120 \
 				--perms 1000 \
 				--filterreads 0 \
 				--qvoffset 33 \
 				--mbq 10 \
 				--mmq 10 \
 				--minc 4 \
 				--VCF ${output}/${project_name}.vcf > ${output}/${project_name}_variantcalls.log
```
