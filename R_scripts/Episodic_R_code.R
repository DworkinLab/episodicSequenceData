library(tidyr)
library(ggplot2)
library(dplyr)

#Data from subset of the Right arm of the 2nd Chromosome
#Random Region near middle of chromosome == 10000 bp
# == 1.8 MB -- Note; have a smaller file (184kb)

##Using smaller file to play with data for now (1000bp)

Data_subset <- read.table("episodic_data_2R_subset.sync")

#Format of populations ex. 27:0:0:30:0:0 == A-count:T-count:C-count:G-count:N-count:deletion-count

#adjust colnames

colnames(Data_subset) <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "AncestorR1_0")

print(head(Data_subset))


# Use tidyr --> make long
Episodic_long <- gather(Data_subset, Population, Allele_Freq , ConR1_115:AncestorR1_0, factor_key=TRUE)

#Error???

#Seperate the allele counts into independent columns for each base
Episodic_seperate <- Episodic_long %>% 
  separate(Allele_Freq, c("A","T","C","G","N","del"), ":")

#Want to rearrange data:

#Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)

Episodic_split <- Episodic_seperate %>%
  separate(Population, c("Treatment", "Generation"), "_")

Episodic_split <- Episodic_split %>%
  separate(Treatment, c("Treatment", "Replicate"), "R")


#as.numeric to multiple columns:
cols.num <- c("Replicate","Generation", "A", "T", "C", "G", "N", "del")
Episodic_split[cols.num] <- sapply(Episodic_split[cols.num],as.numeric)

#Drop Chromosome (All will be 2R)

Episodic_split <- subset(Episodic_split, select = -c(Chromosome) )

# Two columns for counts: Ref-allele counts and Alternative counts 

#Ref_Allele:

Episodic_split$Ref_Allele <- 0

Episodic_split <- within(Episodic_split, {
  Ref_Allele = ifelse (ref == "A", Episodic_split[,6], ifelse (ref == "T", Episodic_split[,7], ifelse (ref == "C", Episodic_split[,8], ifelse (ref == "G", Episodic_split[,9], Episodic_split$Ref_Allele)  )))})


#Alternative Allele:
#First make columns of Alternative Alleles (to make it easier than a bunh of ifelse statements later)

Episodic_split <- within(Episodic_split, {
  Alt_Allele1 = ifelse (ref == "A", (Episodic_split[,8]), Episodic_split[,6])
  Alt_Allele2 = ifelse (ref == "T", (Episodic_split[,8]), Episodic_split[,7])
  Alt_Allele3 = ifelse (ref == "G", (Episodic_split[,8]), Episodic_split[,9])
})

#Now max for each row into Alt_Allele Column
Episodic_split$Alt_Allele <- apply(Episodic_split[,13:15], 1, max)

#Column to say what the alternative Allele is:

Episodic_split <- within(Episodic_split, {
  Alt_Allele_base = ifelse (Alt_Allele == 0, "N/A", ifelse (Alt_Allele == Episodic_split[,6], "A", ifelse (Alt_Allele == (Episodic_split[,7]), "T", ifelse (Alt_Allele == (Episodic_split[,8]), "C", ifelse (Alt_Allele == (Episodic_split[,9]), "G", "N/A" )))))})


# Remove unneeded columns (6,7,8,9,10,11,13,14,15)
Episodic_split <- subset(Episodic_split, select = -c(A,T,C,G,N,del,Alt_Allele3,Alt_Allele2, Alt_Allele1) )


print(head(Episodic_split))
print(tail(Episodic_split))

#Simple plot showing alternative allele counts by postion

p1 <- ggplot(data = Episodic_split, 
             aes(x = Position, y=Alt_Allele, colour=Alt_Allele_base))

p2 <- (p1 + geom_point(size = 2, alpha = 0.5))

print(p2)




p3 <- ggplot(data = Episodic_split, 
             aes(x=Position, y = Alt_Allele, colour= Treatment))


p4 <- (p3 + geom_point(size = 2, alpha = 0.5))

print(p4)

#Allele Frequencies:
Episodic_split$Ref_Allele_freq <- with(Episodic_split, Ref_Allele/(Ref_Allele+Alt_Allele))
head(Episodic_split)
Episodic_split

p5 <- ggplot(data = Episodic_split, 
             aes(x=Position, y = Ref_Allele_freq, colour = Treatment))
p6 <- (p5 + geom_point(size = 2, alpha = 0.5))

print(p6)
