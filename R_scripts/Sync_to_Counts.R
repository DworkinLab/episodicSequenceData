#Test Dir == /home/paul/episodicData/subsetting

# To run: open R (> R) and source this script (be sure to edit based on file used). 

## Convert a .sync file into long format, filter somewhat, and have only position, treatment, Cage, Generation and Maj/Min counts

#Source a packages script with necessary packages needed
#source("packages.R")

### Packages source code: only need these two for this script
require('dplyr')

#tidyr may not work: require on own before running
require('tidyr')


#1) Need to change details as needed above and below string of #####

#2) Needs a .sync file made by popoolation2

#3) Need to change most importantly for analysis the read in and read out names 


#Read in data: subset testing:

#episodic_counts <- read.table("/home/paul/episodicData/subsetting/episodic_data_bowtie_2R_subset.gatk.sync")

#episodic_counts <- read.table("/home/paul/episodicData/subsetting/episodic_data_2R_subset.gatk.sync")

# Read in Data: Big Data Sets

#BWA:
episodic_counts <- read.table("/home/paul/episodicData/R_dir/episodic_data_main.gatk.sync")

#Bowtie:

#episodic_counts <- read.table("/home/paul/episodicData/bowtie/R_bowtie/episodic_data_bowtie_main.gatk.sync")



#adjust colnames

name.Columns <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR2_115", "SelR1_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "SelR1_0")
colnames(episodic_counts) <- name.Columns

#Add "replicates" of ancestor -- all are equal
episodic_counts$SelR2_0 <- episodic_counts$SelR1_0
episodic_counts$ConR1_0 <- episodic_counts$SelR1_0
episodic_counts$ConR2_0 <- episodic_counts$SelR1_0


#Need the ancestor to stay (after making long) to call major/minor alleles later
episodic_counts$Ancestor <- episodic_counts$SelR1_0

# Make long by bring populations down
long_episodic <- gather(episodic_counts, Population, Allele_Freq , ConR1_115:ConR2_0, factor_key=TRUE)

#Error???

# All geneneric below for sync files (only real issue through file is population naming convention)
###################################################

#Seperate the allele counts into independent columns for each base
Episodic_seperate <- long_episodic %>% 
  separate(Allele_Freq, c("A","T","C","G","N","del"), ":")

#Seperate the ancestor to base later things on
Episodic_seperate <- Episodic_seperate %>% 
  separate(Ancestor, c("A_0","T_0","C_0","G_0","N_0","del_0"), ":")

# as.numeric to multiple columns:
cols.num <- c("A_0", "T_0", "C_0", "G_0", "N_0", "del_0", "A", "T", "C", "G", "N", "del")

Episodic_seperate[cols.num] <- sapply(Episodic_seperate[cols.num],as.numeric)

#Get the sum of all the rows (all the different bases) for each population position:

Episodic_seperate$sum <- (rowSums(Episodic_seperate[,11:16]))

#Renamed incase of error:
Episodic_split_2 <- Episodic_seperate

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


# To determine the minor allele base if not specified by the ancestor (new allele brough up etc.)

#max for the population (could be the minor allele)
Episodic_split_2$maj_all <- apply(Episodic_split_2[,11:16], 1, max)

#alt== second highest count for populations
Episodic_split_2$alt_allele <- apply(Episodic_split_2[,11:16], 1, 
                                     function(x)max(x[x!=max(x)]))

Episodic_split_2 <- within(Episodic_split_2, {
  Min_count_2 = ifelse (Maj_count == sum, 0, ifelse(Maj_count==maj_all, alt_allele, maj_all))})

Episodic_split_2 <- within(Episodic_split_2, {
  MinorAllele_base = ifelse(Min_count_2==0, "N/A", ifelse(Min_count_2== Episodic_split_2[,11], "A", ifelse(Min_count_2== Episodic_split_2[,12],  "T", ifelse(Min_count_2== Episodic_split_2[,13],  "C",ifelse(Min_count_2== Episodic_split_2[,14],  "G", ifelse(Min_count_2== Episodic_split_2[,15],  "N", ifelse(Min_count_2== Episodic_split_2[,16],  "del", "Z") ))))))})

# Remove unneeded columns (6,7,8,9,10,11,13,14,15)
Episodic_split_2 <- subset(Episodic_split_2, select = -c(A_0,T_0,C_0,G_0,N_0,del_0,A,T,C,G,N,del,anc_max,anc_min, MinorAllele, Min_count, maj_all, alt_allele))

nam.col <- c("chr", "pos", "ref", "Population", "sum", "MajorAllele", "Major_count", "Minor_count", "MinorAllele")
colnames(Episodic_split_2) <- nam.col


#Remove unneccessary Columns (as needed)
#Keep them all for now (except sum) as may be needed later
#Episodic_split_2 <- subset( Episodic_split_2, select = -ref )
#Episodic_split_2 <- subset( Episodic_split_2, select = -chr)
#Episodic_split_2 <- subset( Episodic_split_2, select = -MajorAllele )
#Episodic_split_2 <- subset( Episodic_split_2, select = -MinorAllele)
Episodic_split_2<- subset( Episodic_split_2, select = -sum)


Episodic_split_3 <- Episodic_split_2


## Depends on the filter method:



#Filter method: take the sum of each position, and must have at least 5 counts called (i.e over the 16 populations, the total of alleles called for the minor allele must be over 5)
grp <- Episodic_split_3 %>%
  group_by(pos) %>%
  summarise(sum=sum(Minor_count))
grp2 <- grp[which(grp$sum<=5),]
Episodic_split_3 <- Episodic_split_3[!(Episodic_split_3$pos %in% grp2$pos),]

#check that the number of obs for episodic_long2 == obs for those without 0's sum (*16 for number of "populations") (needed later as well == grp3)

#grp3 <- grp[-which(grp$sum<=5),]


#################################################
#Should be all genetic above (from start specificed)

## Below depends on the population name layout etc. made above


#Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)

episodic_long <- Episodic_split_3 %>%
  separate(Population, c("Treatment", "Generation"), "_")

episodic_long <- episodic_long %>%
  separate(Treatment, c("Treatment", "Cage"), "R")

cols.num <- c("Cage", "Generation", "Major_count", "Minor_count")
episodic_long[cols.num] <- sapply(episodic_long[cols.num],as.numeric) 


write.csv(episodic_long, file="episodic_bwa_main_counts.csv", row.names = FALSE)
#write.csv(episodic_long, file="episodic_bowtie_main_counts.csv", row.names = FALSE)