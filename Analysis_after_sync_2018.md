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

### Create pileup files for every .bam file

**Script:** [novo_Pi_pileups.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Pi_pileups.sh)

Flags:

- B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment

- Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)

- f -- path to reference sequence

### Run script to calcualte Tajima's Pi using the Variance-sliding.pl script from Popoolation1

**Script:** [novo_tajima_pi.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_tajima_pi.sh)

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

### Bring to local machine

### Create plots of tajima Pi data

R function that will run for each .pi file to output a plot

**Script:** [Pi_PlotFunction.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Pi_PlotFunction.R)

### Run the function for each .pi file

Ex. 
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


