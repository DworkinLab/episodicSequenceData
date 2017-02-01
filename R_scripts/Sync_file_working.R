library(tidyr)
library(ggplot2)
library(dplyr)


#Data from only the Right arm of the 2nd Chromosome
#Data_subset <- read.table("~/Bioinformatics/Sequence_analysis_2016/episodic_data_2R.sync")
# File size == 3.6 Gb --> to large, use a subset

#Random Region near middle of chromosome == 10000 bp
# == 1.8 MB -- Note; have a smaller file (184kb)

Data_subset <- read.table("episodic_data_2R_subset.sync")

#Format of populations ex. 27:0:0:30:0:0 == A-count:T-count:C-count:G-count:N-count:deletion-count

#adjust colnames

colnames(Data_subset) <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "Base")

print(head(Data_subset))

# Use tidyr --> make long
Episodic_long <- gather(Data_subset, Population, Allele_Freq , ConR1_115:Base, factor_key=TRUE)

print(head(Episodic_long))

#Seperate the allele counts into independent columns for each base
Episodic_seperate <- Episodic_long %>% 
  separate(Allele_Freq, c("A","T","C","G","N","del"), ":")

#Order by position
Episodic_seperate <- arrange(Episodic_seperate, Position)

head(Episodic_seperate)

#Make Longer??
longer_Episodic <- gather(Episodic_seperate, Base, Count, A:del, factor_key = TRUE)

head(longer_Episodic)


