# The analysis of pooled sequence data

## After the final .bam files and .sync files have been created. The analaysis here is for data mapped and finalized with three mappers: bwa mem, bowtie2 and novoalign 

See associated scripts for steps up to final bams and creating .sync files
_______________________________________________________________________________________

## Outline of analysis:

### 1) Tajima's Pi of non-overlapping windows for each sequence

### 2) Run Fst on windows for each pairwise comparision of sequenced data

### 3) per SNP logistic regression for each treatment by generation

### 4) estimates of selection coefficient at each position for selection and control lineages
