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

## Notes / additional set up:

- For one mapper: Have a diretory with all the .final.bam files created and .mpileup /.sync files created using these .bam files

- The analysis (up to step 5) will be generic for one mapper (Novoalign) but completed similarily for other chosen mappers

- need to know the order of the .sync files: will be based on the order of the .bam files read into .sync 

- Some scripts require the .sync file split into chromosomes: script below shows method of splitting

**Script:** [novo_split_sync2chromosomes.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_split_sync2chromosomes.sh)
_______________________________________________________________________________________

## 1) Tajima's Pi of non-overlapping windows for each sequence

### Create pileup files for every .bam file

**Script:** [novo_Pi_pileups.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Pi_pileups.sh)

ex.
```
samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_merge_novo_final_realigned.bam > ${output}/${base}.pileup
```

Flags:

- B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment

- Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)

- f -- path to reference sequence

### Run script to calcualte Tajima's Pi using the Variance-sliding.pl script from Popoolation1

**Script:** [novo_tajima_pi.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_tajima_pi.sh)

ex. 
```
perl ${popoolation}/Variance-sliding.pl --input ${input}/${base}.pileup --output ${output}/${base}.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output ${output}/${base}.snps --min-covered-fraction 0.5
```

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

### Create plots of tajima Pi data

On local machine, this R function can run each .pi file to output a plot

**Script:** [Pi_PlotFunction.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Pi_PlotFunction.R)

### In R, run the function for each .pi file

ex. 
```
Pi_PlotFunction('F115ConR1_TAGCTT_novo.pi', "Novoalign")
```


_______________________________________________________________________________________

## 2) Run Fst on windows for each pairwise comparision of sequenced data

### Running Fst

**Script:** [Novo_Fst.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Novo_Fst.sh)

ex.
```
perl ${fst} --input ${novo_mpileup}/novo_episodic_main.sync --output ${novo_fst}/novo_episodic_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```
Flags:

- input -- input sync file
- output -- output file with Fst calculated 
- window-size [500] -- size of the window 
- step-size [500] -- distance to move along chromosome
- min-count [6] -- minimum allele count 
- min-coverage [10] -- minimum coverage
- max-coverage [250] --maximum coverage
- pool-size [120] -- double pooled size (diploid)
- min-covered-fraction [1] -- minimum percentage of sites having sufficient coverage in the given window


### In R, split the file into each compasison

**Script:** [novo_Fst_Split_Comparisons.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Fst_Split_Comparisons.R)

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

3) **Combine all the split chromosome pieces back into one chromosome:** [Combine_chromo.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Combine_chromo.R)

ex.
```
Rscript Combine_chromo.R 'DIRECTORY' 'OutputDIRECTORY'
```



_______________________________________________________________________________________

## 4) estimates of selection coefficient at each position for selection and control lineages




_______________________________________________________________________________________

## 5) Combine outputs of three mappers for Fst, model, and selection coefficients


