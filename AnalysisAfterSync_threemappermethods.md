# Details of procedure for each mapper used

Starting with each mapper having been used and created the .bam files and the large .sync file:

Each .sync files (main chromosomes all together) and the .bam files after GATK relaignment

**Novoalign**
```
novo_episodic_main.sync
```
```
F115ConR1_TAGCTT_novo_merge_novo_final_realigned.bam
F115ConR2_GGCTAC_novo_merge_novo_final_realigned.bam
F115SelR1_GTTTCG_novo_merge_novo_final_realigned.bam
F115SelR2_GTGGCC_novo_merge_novo_final_realigned.bam
F38ConR1_ATCACG_novo_merge_novo_final_realigned.bam
F38ConR2_TTAGGC_novo_merge_novo_final_realigned.bam
F38SelR1_ACTTGA_novo_merge_novo_final_realigned.bam
F38SelR2_GATCAG_novo_merge_novo_final_realigned.bam
F77ConR1_ATGTCA_novo_merge_novo_final_realigned.bam
F77ConR2_ATTCCT_novo_merge_novo_final_realigned.bam
F77SelR1_TTAGGC_novo_merge_novo_final_realigned.bam
F77SelR2_GATCAG_novo_merge_novo_final_realigned.bam
MGD3_SO_CAGATC_novo_merge_novo_final_realigned.bam
```

**Bwa -mem**
```
episodic_data_main.gatk.sync
```
```
F115ConR1_TAGCTT_merged_aligned_pe.final_realigned.bam
F115ConR2_GGCTAC_merged_aligned_pe.final_realigned.bam
F115SelR1_GTTTCG_merged_aligned_pe.final_realigned.bam
F115SelR2_GTGGCC_merged_aligned_pe.final_realigned.bam
F38ConR1_ATCACG_merged_aligned_pe.final_realigned.bam
F38ConR2_TTAGGC_merged_aligned_pe.final_realigned.bam
F38SelR1_ACTTGA_merged_aligned_pe.final_realigned.bam
F38SelR2_GATCAG_merged_aligned_pe.final_realigned.bam
F77ConR1_ATGTCA_merged_aligned_pe.final_realigned.bam
F77ConR2_ATTCCT_merged_aligned_pe.final_realigned.bam
F77SelR1_TTAGGC_merged_aligned_pe.final_realigned.bam
F77SelR2_GATCAG_merged_aligned_pe.final_realigned.bam
MGD3_SO_CAGATC_merged_aligned_pe.final_realigned.bam
```

**Bowtie**

```
episodic_data_bowtie_main.gatk.sync 
```
```
F115ConR1_TAGCTT_merged_bowtie_pe.final_realigned.bam
F115ConR2_GGCTAC_merged_bowtie_pe.final_realigned.bam
F115SelR1_GTTTCG_merged_bowtie_pe.final_realigned.bam
F115SelR2_GTGGCC_merged_bowtie_pe.final_realigned.bam
F38ConR1_ATCACG_merged_bowtie_pe.final_realigned.bam
F38ConR2_TTAGGC_merged_bowtie_pe.final_realigned.bam
F38SelR1_ACTTGA_merged_bowtie_pe.final_realigned.bam
F38SelR2_GATCAG_merged_bowtie_pe.final_realigned.bam
F77ConR1_ATGTCA_merged_bowtie_pe.final_realigned.bam
F77ConR2_ATTCCT_merged_bowtie_pe.final_realigned.bam
F77SelR1_TTAGGC_merged_bowtie_pe.final_realigned.bam
F77SelR2_GATCAG_merged_bowtie_pe.final_realigned.bam
MGD3_SO_CAGATC_merged_bowtie_pe.final_realigned.bam
```

Split the sync file by chromosome:

**Script:** [novo_split_sync2chromosomes.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_split_sync2chromosomes.sh)

Will show the procedure for each mapper to generate files for [outline of methods](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018.md)

_______________________________________________________________________________________

## 1) Tajima's Pi of non-overlapping windows for each sequence

### [Create](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Pi_pileups.sh) pileup files for every .bam file

ex.
```
samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_merge_novo_final_realigned.bam > ${output}/${base}.pileup
```

**Novoalign**

**Bwa -mem**

**Bowtie**


### Run [script](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_tajima_pi.sh) to calcualte Tajima's Pi using the Variance-sliding.pl script from Popoolation1

ex. 
```
perl ${popoolation}/Variance-sliding.pl --input ${input}/${base}.pileup --output ${output}/${base}.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output ${output}/${base}.snps --min-covered-fraction 0.5
```


### [Create plots](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Pi_PlotFunction.R) of tajima Pi data

On local machine, this R function can run each .pi file to output a plot for each mapper

ex. 
```
Pi_PlotFunction('FILE.pi', "Plot_Title_Details-(i.e mapper)")
```

Not done for bwa-mem (worth it??)
_______________________________________________________________________________________

## 2) Run Fst on windows for each pairwise comparision of sequenced data

### Running [Fst](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Novo_Fst.sh)

ex.
```
perl ${fst} --input ${novo_mpileup}/novo_episodic_main.sync --output ${novo_fst}/novo_episodic_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```
**Novoalign**

```
#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
novo_mpileup=${project_dir}/novo_mpileup

#Path and variable for script from PoPoolation to create .sync files
fst=/usr/local/popoolation/fst-sliding.pl

mkdir ${novo_mpileup}/novo_fst
novo_fst=${novo_mpileup}/novo_fst

perl ${fst} --input ${novo_mpileup}/novo_episodic_main.sync --output ${novo_fst}/novo_episodic_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```

**Bwa -mem**
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData


#Path to .sync files
input=${project_dir}/R_dir

#Path and variable for script from PoPoolation to create .sync files
fst=/usr/local/popoolation/fst-sliding.pl

#mkdir ${project_dir}/bwa_fst
#output=${project_dir}/bwa_fst
output=/home/paul

perl ${fst} --input ${input}/episodic_data_main.gatk.sync  --output ${output}/episodic_data_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```

**Bowtie**
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/bowtie

#Path to .sync files
input=${project_dir}/R_bowtie

#Path and variable for script from PoPoolation to create .sync files
fst=/usr/local/popoolation/fst-sliding.pl

#mkdir ${project_dir}/bowtie_fst
#output=${project_dir}/bowtie_fst
output=/home/paul

perl ${fst} --input ${input}/episodic_data_bowtie_main.gatk.sync  --output ${output}/episodic_data_bowtie_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```


Combine here??

### In R, [split](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Fst_Split_Comparisons.R) the file into each compasison

**Works for all three just by changing the working directory and the file read in**
```
# Script to read in one .fst file and split this into many .csv files based on comparisons
# .fst files generated from Popoolation2 fst-sliding.pl script
# Need to change the working directory, the input, and the number of comparisons present (i.e 6:ncol for .fst file)

### Packages Required (tidyverse, but more specifically tidyr)
  require(data.table)
  require(tidyverse)

######### NOVOALIGN ##########
### Set working directory to location of .fst file
  setwd("~/Bioinformatics/episodic_practice/FST/novo_fst")
### Read in the .fst file into R (requires data.table)
  Fst <- fread('novo_episodic_main.fst')
 
######### Bowtie2 ##########
### Set working directory to location of .fst file
  #setwd("~/Bioinformatics/episodic_practice/FST/bowtie_fst")
### Read in the .fst file into R (requires data.table)
  #Fst <- fread('episodic_data_bowtie_main.fst')

######### BWA-MEM ##########
### Set working directory to location of .fst file
  #setwd("~/Bioinformatics/episodic_practice/FST/bwa_fst") 
### Read in the .fst file into R (requires data.table)
  #Fst <- fread('episodic_data_main.fst')

### Make into long format
  XCC  <- gather(Fst, Comparison, Fst_measure, 6:83, factor_key=TRUE)

### Remove intermediate:
  rm(Fst_novo)

### Seperate the Fst (ex. 1:2=na) into a comparison column and a fst column
  novo_Fst <- XCC %>%
    separate(Fst_measure, "=", into = c('Comp', 'Fst'))

### Remove intermediate:
  rm(XCC)

### Remove unnecessary column (column 6 has no value)
  novo_Fst <- novo_Fst[,c(1,2,3,4,5,7,8)]

### Rename columns:
  colnames(novo_Fst) <- c('chr', 'window', "num", 'frac', 'meanCov','Comp', 'Fst')

### Create list of all unique comparisons:
  X_compLIST <- unique(novo_Fst$Comp)

### For loop that will create a .csv file for every comparison:
  for (vxcx in X_compLIST){

    CXV_Comp <- novo_Fst[which(novo_Fst$Comp==vxcx),]
    
    write.csv(CXV_Comp, file=paste("Novo_fst_", vxcx, '.csv', sep = ""), row.names = FALSE)
    }
```

R script that will split the .fst file into many .csv files with each comparison of interest (can choose the necessary ones from here)

### Plotting Fst files for comparisons of interest




_______________________________________________________________________________________

## 3) per SNP logistic regression for each treatment by generation

**Long Script:*** [novo_regression_model_LONGSCRIPT.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_regression_model_LONGSCRIPT.sh)

This script will break the chromosomal .sync files into smaller managable pieces and run through multiple R scripts while removing intermediates:

**R script to covert sync to Count data:** [Sync_to_counts.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Sync_to_counts.R)

ex.
```
Rscript Sync_to_counts.R '${Sync Directory}'
```

**Running the model for each position along the chromosome:** [Counts_to_model.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Counts_to_model.R)

In long script: this is set up to work in parallel, having each chromosome running at the same time (6 instances running over 11 sections)

NOTE: Script needs to be changed to run faster/ more efficiently
ex.
```
Rscript Counts_to_model_2.R 'DIRECTORY'
```

**Combine all the split chromosome pieces back into one chromosome:** [Combine_chromo.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Combine_chromo.R)

ex.
```
Rscript Combine_chromo.R 'DIRECTORY' 'OutputDIRECTORY'
```



_______________________________________________________________________________________

## 4) estimates of selection coefficient at each position for selection and control lineages




_______________________________________________________________________________________

## 5) Combine outputs of three mappers for Fst, model, and selection coefficients

