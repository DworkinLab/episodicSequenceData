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

### Combining the data files from [three mappers](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/FST_combine3mappers.R)

```
## Packages:
  require(data.table)
  require(tidyverse)

## Main directory holding fst directories for mappers:
  setwd("/Users/paulknoops/Bioinformatics/episodic_practice/FST")
  
## List directories in working directories:
  workdir <- getwd()
  mydirs <- list.dirs(path=workdir, recursive = FALSE)
  
## Make list of compasisons to combine (no use doing them all?)
  patty <- c('_fst_1:3.csv', '_fst_2:4.csv', '_fst_5:7.csv', '_fst_6:8.csv', '_fst_9:11.csv', '_fst_10:12.csv')

#For each comparison:
for (patt in patty){
  mtdf <- NULL
  # For each directory holding the comparisons for each mapper (should be in working directory)
  for (dir in mydirs){
    #All files that end with the comparison of interest in the dir (i.e. the one comparison for each mapper)
    mycomp <- list.files(path = dir, pattern=patt, full.names = T)
    
    for (file in mycomp){
      Xc <- fread(file)
      
      Xc <- Xc[-which(Xc$Fst=='na'),]
      
      Xc$Fst <- as.numeric(Xc$Fst)
    
    }
    
    mtdf <- rbind(mtdf, Xc)
  }
  
  # Get a count of the number of occurances of a position (window)
  comp_final <- mtdf %>%
    group_by(chr, window, Comp) %>%
    mutate(count = n())
  
  # Keep those with a count of 3 (three times mapped)
  compcomb <- comp_final[which(comp_final$count==3),]
  
  #Calcualte mean Fst:
  comp_final_2 <- mtdf %>%
    group_by(chr, window, Comp) %>%
    summarise(meanFst = (mean(Fst)))
    
  # Writes many files with the mean FST value of the three mappers for the comparisons of interest
  
  write.csv(comp_final_2, file=paste("combined", patt, sep = ""), row.names = FALSE)
  }
  
```

### Cheap, dirty and quick plot:
```
#Plot data?

require(data.table)
require(tidyverse)


setwd("/Users/paulknoops/Bioinformatics/episodic_practice/FST")
XC2 <- fread('combined_fst_1:3.csv')
CX2 <- fread('combined_fst_2:4.csv')
tttle <- 'F115'

ddat <- rbind(CX2, XC2)
rm(CX2)
rm(XC2)
#ddat
ddat2 <- ddat %>%
  group_by(chr, window) %>%
  summarise(meanFst = (mean(meanFst)))

rm(ddat)

xcs <- mean(ddat2$meanFst)
#median(ddat2$meanFst)

g <- nrow(ddat2[which(ddat2$chr=='2L'),])
h <- nrow(ddat2[which(ddat2$chr=='2R'),])
i <- nrow(ddat2[which(ddat2$chr=='3L'),])
j <- nrow(ddat2[which(ddat2$chr=='3R'),])
k <- nrow(ddat2[which(ddat2$chr=='4'),])
l <- nrow(ddat2[which(ddat2$chr=='X'),])

#To change the order for X to be first:
ddat2$number <- c((l+1):(l+g), 
                  (l+g+1):(l+g+h), 
                  (l+g+h+1):(l+g+h+i),
                  (l+g+h+i+1):(l+g+h+i+j),
                  (l+g+h+i+j+1):(l+g+h+i+j+k), 
                  (1:l))



ggxcv <-  ggplot(data = ddat2, aes(x=number, y=meanFst, color=chr))
ggxcv2 <- ggxcv + 
  geom_point(size=0.2, show.legend = F, alpha = 0.75) + 
  theme(panel.background = element_blank()) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
  geom_hline(yintercept = xcs) + 
  xlab("Chromosome") +
  ggtitle(tttle) +
  scale_x_discrete(limits=c(l/2, l+(g/2), (l+g+(h/2)), (l+g+h+(i/2)), (l+g+h+i+(j/2)), (l+g+h+i+j+(k/2))), labels = c("X","2L", "2R", '3L', '3R', "4")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))

print(ggxcv2)
```
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

```
#! /bin/bash

### Set all variables (need to make an output directory):

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup
	
mkdir ${SyncFiles}/splitsync_dir
splitSync=${SyncFiles}/splitsync_dir

#Output dir:
poolSeq=${project_dir}/novo_PoolSeq

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts
	
# The seperated .sync files
sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 
sync[5]=${SyncFiles}/novo_episodic_4.sync 

##-----------------------------------------------##

### Split into treatment vs. control

for file in ${sync[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$6,$7,$10, $11, $14, $15, $16, $16}' > ${splitSync}/${base}_Sel.sync
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$4,$5,$8, $9, $12, $13, $16, $16}' > ${splitSync}/${base}_Con.sync

done


##------------------------------------------------##

### Split the sync files into many sized files (12):

files=(${splitSync}/*.sync)

for file in ${files[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	
	mkdir ${splitSync}/${base}_Split
	split_sync=${splitSync}/${base}_Split
	
	length=($(wc -l ${splitSync}/${base}.sync))
		
	#Split length into 12 segements (12th == length) (can extend this if to large)
	cut=$((${length}/12))
	cut_2=$((${cut}*2))
	cut_3=$((${cut}*3))
	cut_4=$((${cut}*4))
	cut_5=$((${cut}*5))
	cut_6=$((${cut}*6))
	cut_7=$((${cut}*7))
	cut_8=$((${cut}*8))
	cut_9=$((${cut}*9))
	cut_10=$((${cut}*10))
	cut_11=$((${cut}*11))
	
	sed -n " 1, ${cut} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_1.sync

	sed -n " $((${cut} + 1)), ${cut_2} p"  ${splitSync}/${base}.sync >  ${split_sync}/${base}_2.sync

	sed -n " $((${cut_2} + 1)), ${cut_3} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_3.sync
	
	sed -n " $((${cut_3} + 1)), ${cut_4} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_4.sync

	sed -n " $((${cut_4} + 1)), ${cut_5} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_5.sync

	sed -n " $((${cut_5} + 1)), ${cut_6} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_6.sync

	sed -n " $((${cut_6} + 1)), ${cut_7} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_7.sync
	
	sed -n " $((${cut_7} + 1)), ${cut_8} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_8.sync
	
	sed -n " $((${cut_8} + 1)), ${cut_9} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_9.sync
	
	sed -n " $((${cut_9} + 1)), ${cut_10} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_10.sync
	
	sed -n " $((${cut_10} + 1)), ${cut_11} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_11.sync
	
	sed -n " $((${cut_11} + 1)), ${length} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_12.sync
	
	syncs=(${split_sync}/*.sync)
 	
	for file in ${syncs[@]}
	  	do
	  	(Chromo=$(cat ${file} | awk '{print $1; exit}')
	  	Rscript ${Rscripts}/PoolSeq_SelCoeff.R ${file} ${Chromo} ${split_sync}) &
	done 
	wait
	rm -f ${split_sync}/*.sync
done
wait

Rscript ${Rscripts}/combinePoolseqCSV.R ${splitSync}

##------------------------------------------------##
```


```
#RscriptTest:
args <- commandArgs(trailingOnly = TRUE)

require(methods)
require(data.table)
require(foreach)
require(stringi)
require(matrixStats)

### Source the scripts (Copied) for Pool-Seq (only one fails and is not needed)
source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/loadaf.R')  

source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estsh.R')

source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/idsel.R')

source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/simaf.R')

### Possibly need custom function to read in manipulated .sync files:
### With fiddling with the .sync file, a personal read.sync function is needed

source("/home/paul/episodicData/novoalign/novo_Rscripts/Taus_ReadSync.R")

### Read in the data file for args[1]

setwd(args[3])

mySync <- read.sync_Personal(file=args[1], gen=c(115, 115, 38, 38, 77, 77, 0, 0), repl=c(1,2,1,2,1,2,1,2), polarization = "rising")

# Turn alleles to data frame:
ff <- as.data.frame(mySync@alleles)

# Find length (number of positons)
len <- round(length(ff$posID))

# Keep only positions:
pst <- as.numeric(ff$pos)
pst2 <- sort(pst)

# Generations:
ccc <- c(0,38,77,115)

rm(pst)
rm(ff)

pst2  <- sample(pst2[1]:pst2[len], 1000)

### Create empty matrix to read into for estiamting S:
pbj <- matrix(NA,length(pst2), 3)


  for (i in 1:length(pst2)) {
    b_b <- pst2[i]
    TrajTEST <- af.traj(mySync, args[2], repl=c(1,2), pos=b_b)
    BfsfTEST <- estimateSH(TrajTEST, Ne=150, t=ccc, h=0.5, haploid = FALSE, simulate.p.value=TRUE)
    pbj[i,] <- c(BfsfTEST$s, BfsfTEST$p.value, b_b)
    rm(TrajTEST)
    rm(BfsfTEST)
    rm(b_b)
  }


x2 <- args[1]
x3 <- gsub("\\..*","", x2)
write.csv(pbj, file=paste(args[3], "/", x3, ".csv", sep=""), row.names=FALSE)

rm(pbj)
rm(mySync)
rm(ccc)
rm(pst2)
```

### BWA-mem Pool-Seq

- Rscripts should be the same:
- need to split into Chromosomes first:
```
#!/bin/bash

#Name of the .sync file for spliting:
project_name=episodic_data

#Variable for project:
project_dir=/home/paul/episodicData

#Path to .sync files
SyncFiles=${project_dir}/mpileup_dir

mkdir ${SyncFiles}/splitsync_dir
splitSync=${SyncFiles}/splitsync_dir

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novoalign/novo_Rscripts

grep '3R' ${SyncFiles}/${project_name}_main.gatk.sync > ${SyncFiles}/${project_name}_3R.sync &
grep '2R' ${SyncFiles}/${project_name}_main.gatk.sync > ${SyncFiles}/${project_name}_2R.sync &
grep '3L' ${SyncFiles}/${project_name}_main.gatk.sync > ${SyncFiles}/${project_name}_3L.sync &
grep '2L' ${SyncFiles}/${project_name}_main.gatk.sync > ${SyncFiles}/${project_name}_2L.sync &
grep '^4' ${SyncFiles}/${project_name}_main.gatk.sync > ${SyncFiles}/${project_name}_4.sync &
grep 'X' ${SyncFiles}/${project_name}_main.gatk.sync > ${SyncFiles}/${project_name}_X.sync 

wait

# The seperated .sync files
sync[0]=${SyncFiles}/${project_name}_3R.sync
sync[1]=${SyncFiles}/${project_name}_2R.sync
sync[2]=${SyncFiles}/${project_name}_3L.sync
sync[3]=${SyncFiles}/${project_name}_2L.sync
sync[4]=${SyncFiles}/${project_name}_X.sync 
sync[5]=${SyncFiles}/${project_name}_4.sync 

##-----------------------------------------------##

### Split into treatment vs. control

for file in ${sync[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$6,$7,$10, $11, $14, $15, $16, $16}' > ${splitSync}/${base}_Sel.sync
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$4,$5,$8, $9, $12, $13, $16, $16}' > ${splitSync}/${base}_Con.sync

done


##------------------------------------------------##

### Split the sync files into many sized files (15):

files=(${splitSync}/*.sync)

for file in ${files[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	
	mkdir ${splitSync}/${base}_Split
	split_sync=${splitSync}/${base}_Split
	
	length=($(wc -l ${splitSync}/${base}.sync))
		
	#Split length into 15 segements (15th == length) (can extend this if to large)
	cut=$((${length}/15))
	cut_2=$((${cut}*2))
	cut_3=$((${cut}*3))
	cut_4=$((${cut}*4))
	cut_5=$((${cut}*5))
	cut_6=$((${cut}*6))
	cut_7=$((${cut}*7))
	cut_8=$((${cut}*8))
	cut_9=$((${cut}*9))
	cut_10=$((${cut}*10))
	cut_11=$((${cut}*11))
	cut_12=$((${cut}*12))
	cut_13=$((${cut}*13))
	cut_14=$((${cut}*14))
	
	sed -n " 1, ${cut} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_1.sync

	sed -n " $((${cut} + 1)), ${cut_2} p"  ${splitSync}/${base}.sync >  ${split_sync}/${base}_2.sync

	sed -n " $((${cut_2} + 1)), ${cut_3} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_3.sync
	
	sed -n " $((${cut_3} + 1)), ${cut_4} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_4.sync

	sed -n " $((${cut_4} + 1)), ${cut_5} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_5.sync

	sed -n " $((${cut_5} + 1)), ${cut_6} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_6.sync

	sed -n " $((${cut_6} + 1)), ${cut_7} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_7.sync
	
	sed -n " $((${cut_7} + 1)), ${cut_8} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_8.sync
	
	sed -n " $((${cut_8} + 1)), ${cut_9} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_9.sync
	
	sed -n " $((${cut_9} + 1)), ${cut_10} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_10.sync
	
	sed -n " $((${cut_10} + 1)), ${cut_11} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_11.sync
	
	sed -n " $((${cut_11} + 1)), ${cut_12} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_12.sync
	
	sed -n " $((${cut_12} + 1)), ${cut_13} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_13.sync
	
	sed -n " $((${cut_13} + 1)), ${cut_14} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_14.sync
	
	sed -n " $((${cut_14} + 1)), ${length} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_15.sync
	
	syncs=(${split_sync}/*.sync)
 	
	for file in ${syncs[@]}
	  	do
	  	(Chromo=$(cat ${file} | awk '{print $1; exit}')
	  	Rscript ${Rscripts}/poolSeq_selectionCoeff.R ${file} ${Chromo} ${split_sync}) &
	done 
	wait
	rm -f ${split_sync}/*.sync
done
wait

# Did not work in the other: will do manually later
#Rscript ${Rscripts}/combinePoolseqCSV.R ${splitSync}

##------------------------------------------------##
```

### Bowtie2 



_______________________________________________________________________________________

## 5) Combine outputs of three mappers for Fst, model, and selection coefficients

