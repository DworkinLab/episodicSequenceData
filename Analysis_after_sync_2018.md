# The analysis of pooled sequence data after the final .bam files and .sync files have been created. 

- The analaysis here is for data mapped and finalized with three mappers: **bwa mem**, **bowtie2** and **novoalign** 

- See associated scripts for steps up to final bams and creating .sync files
_______________________________________________________________________________________

## Outline of analysis:

### 1) Tajima's Pi of non-overlapping windows for each sequence
/
/
### 2) Run Fst on windows for each pairwise comparision of sequenced data and calculate average Fst across three mappers
/
/
### 3) per SNP logistic regression for each treatment by generation averaged for Novoalign, Bwa-mem and bowtie2
/
/
### 4) Comparison of Fst and model output
/
/
### 5) Average estimates of selection coefficient at each position for selection and control lineages for thee mappers
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
Pi_PlotFunction('FILE.pi', "Plot_Title_Details")
```

### Example Output: Ancestral Pi for Novoalign

![Ancestral Pi Plot for Novoalign](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Ancestral_Pi.png)


Currently, Pi files for Novoalign and Bowtie2 generated.

- Should these be averaged (with bwa-mem that can be generated easily enough)

- is one mapper enough for this

- Do we want to show trends in nucleotide diversity (other Pi plots) or just focus on ancestoral Pi

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

### Combining three mappers:

**Script**[FST_combine3mappers.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/FST_combine3mappers.R)

Will take the split comparisons, and combine those specified into one FST file with the average Fst between the three mappers


### Plotting Fst files for comparisons of interest

 - average Fst of three mappers **and** average between replicates
 
 - comparison betweeen control and selection lines
 
Generation 38:
![meanFst for F38](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F38_meanFstPlot.png)

Generation 77: 
![meanFst for F77](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F77_meanFstPlot.png)

Generation 115: 
![meanFst for F115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_meanFstPlot.png)




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
Basic Model at each positon (tmp2):
```
modlist_2[[i]] <- 
        glm(cbind(Major_count, Minor_count) ~ Treatment*Generation, 
            data = tmp2, family = "binomial")
```

**Combine all the split chromosome pieces back into one chromosome:** [Combine_chromo.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Combine_chromo.R)

ex.
```
Rscript Combine_chromo.R 'DIRECTORY' 'OutputDIRECTORY'
```

**Combine three mappers into one file** [model_combine3mappers.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/model_combine3mappers.R)
```
Rscript combine_threeMappers.R 'OutputDirectory
```

**Write files with coeffefficent of interest** [model_3mappersTxG.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/model_3mappersTxG.R)

Keeps positions that are present in all three files (position needs to be mapped three times)
```
Rscript model_3mappersTxG.R 'Input/OutputDirectory'
```

**Plots**
Which size is better????

Treatment x Generation -log10(meanP-value) for model output: 
![FullGenomeTxGPlot](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/CHROMO_meanP.png)

Treatment x Generation -log10(meanP-value) for model output: colours dictate the qunatiles (the top 1% in dark grey, green being the top 0.01% of positions)
![QuantilePlot](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/CHROMOs_Qunatiles_D1.png)

Just the top 1% of positions (~ 32000 positions)
![top1%](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/model_top1percent.png)

## 4) Comparison of Fst (of generation 115, control:Selection) and Model Outputs:

**First draft of plots** 

**Full Chromosome:**

 -- Numbers/ peaks do not line up perfectly (product of how they are plotted/ measured (windows vs. positions)
 
![FstV.model_full](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/log10p_fst_comboFUll.png)

 -- Individual chromosome plots line up per position, but the fst have 1/500th (ish) data points (windows of 500 used)
 
**X Chromosome**
![XChromo](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fst_pvalue_X.png)

**2L Chromosome**
![2LChromo](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fst_pvalue_2L.png)

**2R Chromosome**
![2RChromo](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fst_pvalue_2R.png)

**3L Chromosome**
![3LChromo](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fst_pvalue_3L.png)

**3R Chromosome**
![3RChromo](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fst_pvalue_3R.png)

**4 Chromosome**
![4Chromo](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fst_pvalue_4.png)


**The 2L peak seems consistant and present in both Fst and Model output:**

2L:4991444..5304468 in [flybase.org](http://flybase.org/cgi-bin/gbrowse2/dmel/?id=8a57f1c89e508abfc4b4a39afdac905e&snapname=snap_-onDblclick&snapcode=bcfe8ef3c884aad6ffaaecd5b8d853e5&source=dmel)


_______________________________________________________________________________________

## 5) estimates of selection coefficient at each position for selection and control lineages using [poolSeq](https://github.com/ThomasTaus/poolSeq) R package:

Not completed currently:

- Machine was full and current script takes a long time

- should this be done on the full chromosome or focus on regions of interest (say the top 1% above)

**Script:** [poolseq_SelectionCoefficientEstimate.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/poolseq_SelectionCoefficientEstimate.sh)

This script will break the file into two treatment .sync files, break apart these .sync files (smaller sized files), and run through a R script to run poolSeq Package. 

**Rscript:** [poolSeq_selectionCoeff.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/poolSeq_selectionCoeff.R)

Combine the .csv into one file:
**Rscript:** [combinePoolseqCSV.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/combinePoolseqCSV.R)

Notes: 

- need to edit to make more efficient to run (not taking 3 days per section)

- poolSeq may not work on certain versions of R; but can bring in Taus poolSeq scripts in manually and source (done in this script)

- may need to ensure updated packages and install if necessary (i.e matrixStats_0.53.0 installed for this reason)

- breaking the .sync file causes changes in structure, so a modified read_sync function is used (in R script; Taus_ReadSync.R))

_______________________________________________________________________________________



