# The analysis of pooled sequence data after the final .bam files and .sync files have been created. 

- The analaysis here is for data mapped and finalized with three mappers: **bwa mem**, **bowtie2** and **novoalign** 

- See associated scripts for steps up to final bams and creating .sync files
_______________________________________________________________________________________

## Outline of analysis:

### 1) Tajima's Pi of non-overlapping windows for each sequence
/
/
### 2) Run Fst on windows for each pairwise comparision of sequenced data
/
/
### 3) per SNP logistic regression for each treatment by generation
/
/
### 4) estimates of selection coefficient at each position for selection and control lineages
/
/
### 5) Combine outputs of three mappers for Fst, model, and selection coefficients
/
/

_______________________________________________________________________________________

## Notes:

- For one mapper: Have a diretory with all the .final.bam files created and .mpileup /.sync files created using these .bam files

- The analysis (up to step 5) will be generic for one mapper (Novoalign) but completed similarily for other chosen mappers


_______________________________________________________________________________________

## 1) Tajima's Pi of non-overlapping windows for each sequence
1) Create Directories

Will need two new directories, one for pileup files, and one for the output Pi files: Be sure they are in project dir

```
mkdir novo_pileup
mkdir novo_pi
```

2) Create pileup files of every .bam file

Each generated .bam file needs to be in mpileup format for the use in Popoolation Scripts:

Flags:

- B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment

- Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)

- f -- path to reference sequence
        
Script: novo_Pi_pileups.sh
```
#! /bin/bash

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_GATK

# Path to output novoalign pileup files
output=${project_dir}/novo_pileup

index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

files=(${input}/*_merge_novo_final_realigned.bam)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _merge_novo_final_realigned.bam`

samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_merge_novo_final_realigned.bam > ${output}/${base}.pileup

done
```

3) Run script to calcualte Tajima's Pi

Using the Variance-sliding.pl script from Popoolation1

Flags:

- input -- input pileup file

- output -- output file with Tajima's Pi calculated

- measure [pi] -- Options include Tajima's Pi or Wattersons Theta or Tajima's D along chromosomes using a sliding window approach

- window-size [10000] -- size of the sliding window 

- step-size [10000] -- how far to move along with chromosome (if step size smaller, windows will overlap)
 
- min-count [2] -- minimum allele count 

- min-coverage [4] -- minimum coverage (not important if subsampling done..)

- max-coverage [400] --maximum coverage

- min-qual [20] -- minimum base quality (already filtered for 20 multiple times)

- pool-size [120] -- number of chromosomes (So double the number of individuals per pool)

- fastq-type [sanger] -- depending on the encoding of the fastq files

- min-covered-fraction [0.5] -- minimum percentage of sites having sufficient coverage in the given window -- 0.5 from example


Script: novo_tajima_pi.sh

```
#! /bin/bash

# Path to PoPoolation1 (Currently in Paul's Home directory)
popoolation=/home/paul/popoolation_1.2.2

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_pileup

# Path to output Tajima Pi files
output=${project_dir}/novo_pi

files=(${input}/*.pileup)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} .pileup`

perl ${popoolation}/Variance-sliding.pl \
	--input ${input}/${base}.pileup \
	--output ${output}/${base}.pi \
	--measure pi \
	--window-size 10000 \
	--step-size 10000 \
	--min-count 2 \
	--min-coverage 4 \
	--max-coverage 400 \
	--min-qual 20 \
	--pool-size 120 \
	--fastq-type sanger \
	--snp-output ${output}/${base}.snps \
	--min-covered-fraction 0.5
done
```
Need to check the files were created properly (check that there are no empty files).


4) Bring to local machine
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_pi/*.pi /Users/paulknoops/Bioinformatics/R-projects_git/episodicSequenceData/R_scripts/Pi_Analysis_Novo
```

5) Function to create plots of tajima Pi data

R Function to run and create plots [FUNTION](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Pi_PlotFunction.R) 

See: R_scripts/Pi_Analysis_Novo/Pi_PlotFunction.R


6) Run the function for each .pi file

```
Pi_PlotFunction('F115ConR1_TAGCTT_novo.pi', "Novoalign")
```





_______________________________________________________________________________________

## 2) Run Fst on windows for each pairwise comparision of sequenced data




_______________________________________________________________________________________

## 3) per SNP logistic regression for each treatment by generation




_______________________________________________________________________________________

## 4) estimates of selection coefficient at each position for selection and control lineages




_______________________________________________________________________________________

## 5) Combine outputs of three mappers for Fst, model, and selection coefficients


