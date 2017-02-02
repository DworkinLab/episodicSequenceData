library(tidyr)
library(ggplot2)
library(dplyr)


#Data from only the Right arm of the 2nd Chromosome
#Data_subset <- read.table("~/Bioinformatics/Sequence_analysis_2016/episodic_data_2R.sync")
# File size == 3.6 Gb --> to large, use a subset

#Random Region near middle of chromosome == 10000 bp
# == 1.8 MB -- Note; have a smaller file (184kb)


##Use smaller file to play with data for now!
Data_subset <- read.table("episodic_data_2R_1000_subset.sync")

#Format of populations ex. 27:0:0:30:0:0 == A-count:T-count:C-count:G-count:N-count:deletion-count

#adjust colnames

colnames(Data_subset) <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "Gen0")

print(head(Data_subset))
# Use tidyr --> make long
Episodic_long <- gather(Data_subset, Population, Allele_Freq , ConR1_115:Gen0, factor_key=TRUE)

head(Episodic_long)
#Seperate the allele counts into independent columns for each base
Episodic_seperate <- Episodic_long %>% 
  separate(Allele_Freq, c("A","T","C","G","N","del"), ":")
head(Episodic_seperate)
#Make Longer
longer_Episodic <- gather(Episodic_seperate, Base, Count, A:del, factor_key = TRUE)
head(longer_Episodic)
#numeric counts and arrange by position
longer_Episodic$Count <- as.numeric(longer_Episodic$Count)
longer_Episodic <- arrange(longer_Episodic, Position)
h
#Remove all those with count == 0
Epi_rem <- subset(longer_Episodic, Count>0)
head(Epi_rem)

#Plot for base generation: first subset

Epi_Base <- subset(Epi_rem, Population =="Gen0")
Epi_Base <- arrange(Epi_Base, Position)
head(Epi_Base)


p1 <- ggplot(data = Epi_Base, 
             aes(x = Position, y=Count, colour=Base))
p2 <- ggplot(data = Epi_Base, 
                   aes(x = Position, y=ref, colour=Base))
p3 <- (p1 + geom_point(size = 2, alpha=0.5))
p4 <- (p2 + geom_point(alpha=0.5))
print(p3)
print(p4)
