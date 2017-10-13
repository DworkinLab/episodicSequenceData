# Analysis with custom R scripts
________________________________________________________________________________________________________________

Layout of work done with R on generated .sync files

## Before R: Setting Up .sync Files

The .sync files are generated using Popoolation2 depending on the data cleaning done to generate .bam files 

The .sync files used can have the unused regions removed for R analysis: use grep
```
#Remove Heterochromatin regions ('Het'), the Unmapped Regions ('U'), and the Mitochondria genome ('dmel_mitochondrion_genome'):
# -v == does not contain (keeps those without the 'Het' etc.)
# ${sync} == path to the location of sync files
# ${project_name} == name of the .sync file that uniquely identifies it 

grep -v 'Het' ${sync}/${project_name}.sync > ${sync}/${project_name}_less_het.sync
grep -v 'U' ${sync}/${project_name}_less_het.sync > ${sync}/${project_name}_removed_U_Het.sync
grep -v 'dmel_mitochondrion_genome' ${sync}/${project_name}_removed_U_Het.sync > ${sync}/${project_name}_main.sync

#Remove intermediates (-f == force remove) to have main.sync files
rm -f ${sync}/${project_name}_less_het.sync
rm -f ${sync}/${project_name}_removed_U_Het.sync
```

Left with a .sync file with only the 3R, 3L, 2R, 2L, X and 4 chromosomal arms

Split these more: split into individual files based on chromosome using grep
```
# ${sync} == path to the location of sync files
# ${project_name} == name of the .sync file that uniquely identifies it 
# '&' runs each in parallel: all will run together not one at a time
# For 4th chromo: using ^4 which only keeps rows that begin with 4 (layout of .sync begins with chromosome)

grep '3R' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_3R.sync &
grep '2R' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_2R.sync &
grep '3L' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_3L.sync &
grep '2L' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_2L.sync &
grep '^4' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_4.sync &
grep 'X' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_X.sync
```

.sync files for each different chromsome

Depending on the size, these may be able to run individually (likely possible with small 4th) but R only works if using smaller data sets: may need to break up further

Each will be in own directory (base_dir) for each chromosome to keep track of easier

This can be adjusted based on the length of each file: aim for ~2,000,000 positions per section: Example below for 10 sections

Note: If one file is much longer (i.e generally 3R) it may be more beneficial to do each file seperatly to make more or less sections based on each length (i.e 10 sections for 2R, 13 for 3R).

First need a directory for the subset base_dir directories: change based on chosen path etc.

```
mkdir /home/paul/episodicData/bowtie/R_bowtie/bowtie_subset
```

```
#! /bin/bash

#Variable location (location of the split sync files)
SyncFiles=/home/paul/episodicData/bowtie/R_bowtie

#Variable for output for subsets
base_dir=/home/paul/episodicData/bowtie/R_bowtie/bowtie_subset

#Variable for each different .sync chromosome:

sync[0]=${SyncFiles}/episodic_data_bowtie_3R.gatk.sync
sync[1]=${SyncFiles}/episodic_data_bowtie_2R.gatk.sync
sync[2]=${SyncFiles}/episodic_data_bowtie_3L.gatk.sync
sync[3]=${SyncFiles}/episodic_data_bowtie_2L.gatk.sync
sync[4]=${SyncFiles}/episodic_data_bowtie_X.gatk.sync

#Don't need to do 4th chromo: Small enough to run on own
#sync[5]=${SyncFiles}/episodic_data_bowtie_4.gatk.sync


# Loop for each sync file into 10 (or more) sections depending on the ultimate length of the file to have ~2,000,000
# Length will go to the end: so ensure that the final section is not to long (alternative for more sections hashed out below)

for file in ${sync[@]}
do
name=${file}
base=`basename ${name} .gatk.sync`

mkdir ${base_dir}/${base}_dir
basedir=${base_dir}/${base}_dir

length=($(wc -l ${SyncFiles}/${base}.gatk.sync))

sed -n ' 1, 2054752 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_1.sync
sed -n ' 2054753, 4109504 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_2.sync
sed -n ' 4109505, 6164256 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_3.sync
sed -n ' 6164257, 8219008 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_4.sync
sed -n ' 8219009, 10273760 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_5.sync
sed -n ' 10273761, 12328512 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_6.sync
sed -n ' 12328513, 14383264 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_7.sync
sed -n ' 14383265, 16438016 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_8.sync
sed -n ' 16438017, 18492768 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_9.sync
sed -n " 18492769, ${length} p" ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_10.sync

#Adding an 11th section: Hash out _10 above and unhash below (can continue trend past these for more or less sections)
#sed -n ' 18492769, 20547500 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_10.sync
#sed -n " 20547501, ${length} p" ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_11.sync 

done
``` 
Need to make a directory for the 4th chromosome (same location of the basedir's above) and move the single file for 4th chromosome there



Now have many .sync files in sized able to be used, each in seperate directory within one larger directory (and should be only directories within large directory): this is important for later!

## Running In R

### Packages:

Packages Needed through analysis: be sure to have installed on either local or remote machine based on place of analysis

```
#install.packages('dplyr')
library(dplyr)

#install.packages('tidyr')
library(tidyr)

#install.packages('ggplot2')
library(ggplot2)
```

### Converting Sync file to counts

Method to create files with 

This process is quite slow: will do each file (2,000,000 positions) seperatly: look into parallizing script

To run: open R (> R) and source this script (be sure to edit based on file used). 

Convert a .sync file into long format, filter somewhat, and have only position, treatment, Cage, Generation and Maj/Min counts

```

# For loop sync to counts

# To run: open R (> R) and source this script (be sure to edit based on file used). 

## Convert a .sync file into long format, filter somewhat, and have only position, treatment, Cage, Generation and Maj/Min counts

### Packages source code: only need these two for this script
require('dplyr')

#tidyr may not work: require on own before running
require('tidyr')

#1) Need to change details as needed above and below string of #####

#2) Needs a .sync file made by popoolation2

#3) Need to change most importantly for analysis the read in and read out names 

# Read in Data: Big Data Sets
# Pwd a directory containing only the directories of interest (made with other sed -n script)

mydirs <- list.dirs(path = "/home/paul/episodicData/R_dir/bwa_subsetDirectories", recursive = FALSE)

#to not include that actual dir, recursive = FALSE
#Will go through each directory:
for (dir in mydirs){

    setwd(dir)
    
  #within each directory: go through each .sync subset:
  
  mysyncs <- list.files(pattern=".sync")
  
  for (sync in mysyncs){
  
      #Read in the date:
      episodic_counts <- read.table(sync)
    
    #adjust colnames
    print("data read in")
    name.Columns <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR2_115", "SelR1_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "SelR1_0")
    
    colnames(episodic_counts) <- name.Columns
    
    #Add "replicates" of ancestor -- all are equal
    episodic_counts$SelR2_0 <- episodic_counts$SelR1_0
    episodic_counts$ConR1_0 <- episodic_counts$SelR1_0
    episodic_counts$ConR2_0 <- episodic_counts$SelR1_0
      
    #Need the ancestor to stay (after making long) to call major/minor alleles later
    episodic_counts$Ancestor <- episodic_counts$SelR1_0
    
    # Make long by bring populations down
    print("making long")
    long_episodic <- gather(episodic_counts, Population, Allele_Freq , ConR1_115:ConR2_0, factor_key=TRUE)
    
    #Remove intermediates to save storage:
    
    rm(episodic_counts)
    print("removed counts")
    
    # All geneneric below for sync files (only real issue through file is population naming convention)
    
    ###################################################
    
    #Seperate the allele counts into independent columns for each base
    print("splitting allele freq")
    Episodic_split_2 <- long_episodic %>% 
      separate(Allele_Freq, c("A","T","C","G","N","del"), ":")
    
    rm(long_episodic)
    
    print("removed long")
    
    #Seperate the ancestor to base later things on
    Episodic_split_2 <- Episodic_split_2 %>% 
      separate(Ancestor, c("A_0","T_0","C_0","G_0","N_0","del_0"), ":")
    
    # as.numeric to multiple columns:
    cols.num <- c("A_0", "T_0", "C_0", "G_0", "N_0", "del_0", "A", "T", "C", "G", "N", "del")
    
    #Seems to take a long time for this step?
    Episodic_split_2[cols.num] <- sapply(Episodic_split_2[cols.num],as.numeric)
    
    #Get the sum of all the rows (all the different bases) for each population position:
    
    print("getting row sums")
    Episodic_split_2$sum <- (rowSums(Episodic_split_2[,11:16]))
    
    #Ancestor Major_Allele and minor allele:
    
    # Major allele of ancestor == the maximum positional count
    Episodic_split_2$anc_max <- apply(Episodic_split_2[,4:9], 1, max)
    # Minor is the ancestor second highest count
    Episodic_split_2$anc_min <- apply(Episodic_split_2[,4:9], 1, 
                                      function(x)max(x[x!=max(x)]))
    
    #Major / Minor Base name: match the number of anc_max with the column to call the correct base:
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MajorAllele = ifelse(anc_max== Episodic_split_2[,4], "A", ifelse(anc_max== Episodic_split_2[,5],  "T", ifelse(anc_max== Episodic_split_2[,6],  "C",ifelse(anc_max== Episodic_split_2[,7],  "G", ifelse(anc_max== Episodic_split_2[,8],  "N", ifelse(anc_max== Episodic_split_2[,9],  "del", "N/A" ))))))})
    
    #Major Allele Count of evolved populations; match the Major allele with the count of certain columns for each population 
    
    Episodic_split_2 <- within(Episodic_split_2, {
      Maj_count = ifelse (MajorAllele == "A", Episodic_split_2[,11], ifelse (MajorAllele == "T", Episodic_split_2[,12], ifelse (MajorAllele == "C", Episodic_split_2[,13], ifelse (MajorAllele == "G", Episodic_split_2[,14], ifelse (MajorAllele == "N", Episodic_split_2[,15], ifelse (MajorAllele == "del", Episodic_split_2[,16], "N/A"))))))})
    
    
    # Same thing for minor allele: first ensure that if the sum of all counts == the Major coutn and the ancestor had no minor allele, their is no minor allele (N/A), then follow the same match of anc_min to a certain base
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MinorAllele = ifelse(Maj_count==Episodic_split_2[,17] & anc_min==0, "N/A", ifelse(anc_min== Episodic_split_2[,4], "A", ifelse(anc_min== Episodic_split_2[,5],  "T", ifelse(anc_min== Episodic_split_2[,6],  "C",ifelse(anc_min== Episodic_split_2[,7],  "G", ifelse(anc_min== Episodic_split_2[,8],  "N", ifelse(anc_min== Episodic_split_2[,9],  "del", "Z") ))))))})
    
    
    #Minor Allele Count of the ancestreal minor allele count
    Episodic_split_2 <- within(Episodic_split_2, {
      Min_count = ifelse (MinorAllele == "A", Episodic_split_2[,11], ifelse (MinorAllele == "T", Episodic_split_2[,12], ifelse (MinorAllele == "C", Episodic_split_2[,13], ifelse (MinorAllele == "G", Episodic_split_2[,14], ifelse (MinorAllele == "N", Episodic_split_2[,15],ifelse (MinorAllele == "del", Episodic_split_2[,16],"N/A"))))))})
    
    print("called major and minor alleles and counts")
    # To determine the minor allele base if not specified by the ancestor (new allele brough up etc.)
    
    #max for the population (could be the minor allele)
    Episodic_split_2$maj_all <- apply(Episodic_split_2[,11:16], 1, max)
    
    #alt== second highest count for populations
    Episodic_split_2$alt_allele <- apply(Episodic_split_2[,11:16], 1, 
                                         function(x)max(x[x!=max(x)]))
    
    print("define unknown alleles")
    Episodic_split_2 <- within(Episodic_split_2, {
      Min_count_2 = ifelse (Maj_count == sum, 0, ifelse(Maj_count==maj_all, alt_allele, maj_all))})
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MinorAllele_base = ifelse(Min_count_2==0, "N/A", ifelse(Min_count_2== Episodic_split_2[,11], "A", ifelse(Min_count_2== Episodic_split_2[,12],  "T", ifelse(Min_count_2== Episodic_split_2[,13],  "C",ifelse(Min_count_2== Episodic_split_2[,14],  "G", ifelse(Min_count_2== Episodic_split_2[,15],  "N", ifelse(Min_count_2== Episodic_split_2[,16],  "del", "Z") ))))))})
    
    # Remove unneeded columns (6,7,8,9,10,11,13,14,15)
    Episodic_split_2 <- subset(Episodic_split_2, select = -c(A_0,T_0,C_0,G_0,N_0,del_0,A,T,C,G,N,del,anc_max,anc_min, MinorAllele, Min_count, maj_all, alt_allele))
    
    print("removed unneeded columns")
    nam.col <- c("chr", "pos", "ref", "Population", "sum", "MajorAllele", "Major_count", "Minor_count", "MinorAllele")
    colnames(Episodic_split_2) <- nam.col
    
    
    #Remove unneccessary Columns (as needed)
    #Keep them all for now (except sum) as may be needed later
    #Episodic_split_2 <- subset( Episodic_split_2, select = -ref )
    #Episodic_split_2 <- subset( Episodic_split_2, select = -chr)
    #Episodic_split_2 <- subset( Episodic_split_2, select = -MajorAllele )
    #Episodic_split_2 <- subset( Episodic_split_2, select = -MinorAllele)
    Episodic_split_2<- subset( Episodic_split_2, select = -sum)
    
 #%%%%%%%%%%%%%%%%%%%%%%
 
    ## Depends on the filter method:
    print("begin filtering")
    
    #Filter method: take the sum of each position, and must have at least 5 counts called (i.e over the 16 populations, the total of alleles called for the minor allele must be over 5)
    
    grp <- Episodic_split_2 %>%
      group_by(pos) %>%
      summarise(sum=sum(Minor_count))
    grp2 <- grp[which(grp$sum<=5),]
    Episodic_split_2 <- Episodic_split_2[!(Episodic_split_2$pos %in% grp2$pos),]
    
    #check that the number of obs for episodic_long2 == obs for those without 0's sum (*16 for number of "populations") (needed later as well == grp3)
    
    #grp3 <- grp[-which(grp$sum<=5),]
    rm(grp)
    rm(grp2)
    
    print("remove filter inermediates")
    
 #%%%%%%%%%%%%%%%%%%%%%%%%
 
    #################################################
    #Should be all genetic above (from start specificed)
    
    ## Below depends on the population name layout etc. made above
    
    #Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)
    
    print("seperate population to Treatment, Generation and Cage")
    episodic_long <- Episodic_split_2 %>%
      separate(Population, c("Treatment", "Generation"), "_")
    
    rm(Episodic_split_2)
    
    episodic_long <- episodic_long %>%
      separate(Treatment, c("Treatment", "Cage"), "R")
    
    cols.num <- c("Cage", "Generation", "Major_count", "Minor_count")
    episodic_long[cols.num] <- sapply(episodic_long[cols.num],as.numeric) 
    
    print("Have final episodic long; now write a csv")
  
    #will need to rename .csv files  
    write.csv(episodic_long, file=paste(sync, ".csv", sep=""))
  
    
    print("wrote csv and now done this .sync file")
  }
}

```

Have a .csv file that can be used to running the model:

