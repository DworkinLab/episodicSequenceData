library(tidyr)
library(ggplot2)
library(dplyr)

#To Large
#Data_subset <- read.table("~/Bioinformatics/Sequence_analysis_2016/episodic_data_2R.sync")
# == 3.6 Gb

#Random Region near middle of chromosome == 10000 bp
Data_subset <- read.table("episodic_data_2R_subset.sync")
# == 1.8 MB

head(Data_subset)
colnames(Data_subset) <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "Base")
head(Data_subset)

#Use tidyr --> make long?
Episodic_long <- gather(Data_subset, Population, Allele_Freq , ConR1_115:Base, factor_key=TRUE)
Episodic_long
