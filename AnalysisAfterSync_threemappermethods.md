Starting with each mapper having been used and created the .bam files and the large .sync file:

Each Mapper files and layout:

**Novoalign**

**Bwa -mem**

**Bowtie**

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

**Bowtie**
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



### In R, [split](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Fst_Split_Comparisons.R) the file into each compasison


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

