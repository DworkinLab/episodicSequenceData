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

colnames(Data_subset) <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "BaseR1_0")

print(head(Data_subset))
# Use tidyr --> make long
Episodic_long <- gather(Data_subset, Population, Allele_Freq , ConR1_115:BaseR1_0, factor_key=TRUE)

print(head(Episodic_long))
#Seperate the allele counts into independent columns for each base
Episodic_seperate <- Episodic_long %>% 
  separate(Allele_Freq, c("A","T","C","G","N","del"), ":")

### Do not need next bit, keeping just in case
##   Make Longer --

#longer_Episodic <- gather(Episodic_seperate, Base, Count, A:del, factor_key = TRUE)

##   numeric counts and arrange by position

#longer_Episodic$Count <- as.numeric(longer_Episodic$Count)
#longer_Episodic <- arrange(longer_Episodic, Position)

##    Remove all those with count == 0
#Epi_rem <- subset(longer_Episodic, Count>0)
#   Plot for base generation: first subset
#Epi_Base <- subset(Epi_rem, Population =="BaseR1_0")
#Epi_Base <- arrange(Epi_Base, Position)
#p1 <- ggplot(data = Epi_Base, 
#             aes(x = Position, y=Count, colour=Base))
#p2 <- ggplot(data = Epi_Base, 
#                   aes(x = Position, y=ref, colour=Base))
#p3 <- (p1 + geom_point(size = 2, alpha=0.5))
#p4 <- (p2 + geom_point(alpha=0.5))
#print(p3)
#print(p4)

### Re-start work from here

#Want to rearrange data: start with Episodic_seperate

print(head(Episodic_seperate))

#Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)
Episodic_split <- Episodic_seperate %>%
  separate(Population, c("Treatment", "Generation"), "_")

Episodic_split <- Episodic_split %>%
  separate(Treatment, c("Treatment", "Replicate"), "R")


#as.numeric to multiple columns:
cols.num <- c("Replicate","Generation", "A", "T", "C", "G", "N", "del")
Episodic_split[cols.num] <- sapply(Episodic_split[cols.num],as.numeric)

#Drop Chromosome
Episodic_split <- subset(Episodic_split, select = -c(Chromosome) )

print(head(Episodic_split))


# Two columns for counts: Ref-allele counts and Alternative counts 

# Call only the columns with counts
#Episodic_split[,6:11]

# Need to find method that if the column name == ref; put count into column, and then find 2nd allele (2nd highest)
    # what if two alleles are both not ref????

# Make long (TidyR), if ref == Base, counts = Epidosic_split$Ref_allele??
    # how to call other allele?
#longer_Episodic <- gather(Episodic_split, Base, Count, A:del, factor_key = TRUE)

#Change values in ref to 6,7,8,9 (the column values of ATCG in Epidosic_split)
#Episodic_split <- within(Episodic_split, ref <- factor(ref, labels = c(6, 8, 9, 7)))
#Episodic_split$ref <- as.integer(Episodic_split$ref)
#head(Episodic_split)

# Take max from Episodic_split[,6:9]
#Episodic_split <- for (i in 1:nrow(Episodic_split)) +

#Closer but NO
#i = 1
#while (i <= nrow(Episodic_split)) {
#Episodic_split$Major_Allele = max(Episodic_split[i,6:9])
#i = i+1}
#Episodic_split$Ref_Allele <- 0
#for (i in 1:nrows(Episodic_split)) {
#  if (Episodic_split$ref == "A" {Episodic_split$Ref_allele <- Episodic_split[i,6]}) else
#  if (Episodic_split$ref == "T"{Episodic_split$Ref_allele <- Episodic_split[i,7]}) else
#  if (Episodic_split$ref == "C" {Episodic_split$Ref_allele <- Episodic_split[i,8]}) else
#  if (Episodic_split$ref == "G" {Episodic_split$Ref_allele <- Episodic_split[i,9]}) 
#}

# THIS WORKS FOR REF_ALLELE! -- Long version
#Episodic_split$Ref_Allele <- 0

#Episodic_split <- within(Episodic_split, {
# Ref_Allele = ifelse (Episodic_split$ref == "A", Episodic_split[,6], Episodic_split$Ref_Allele) } )
#Episodic_split <- within(Episodic_split, {
#  Ref_Allele = ifelse (Episodic_split$ref == "T", Episodic_split[,7], Episodic_split$Ref_Allele) } )
#Episodic_split <- within(Episodic_split, {
#  Ref_Allele = ifelse (Episodic_split$ref == "C", Episodic_split[,8], Episodic_split$Ref_Allele) } )
#Episodic_split <- within(Episodic_split, {
#  Ref_Allele = ifelse (Episodic_split$ref == "G", Episodic_split[,9], Episodic_split$Ref_Allele) } )


#Short version !

Episodic_split$Ref_Allele <- 0

Episodic_split <- within(Episodic_split, {
  Ref_Allele = ifelse (ref == "A", Episodic_split[,6], ifelse (ref == "T", Episodic_split[,7], ifelse (ref == "C", Episodic_split[,8], ifelse (ref == "G", Episodic_split[,9], Episodic_split$Ref_Allele)  )))})

head(Episodic_split)
tail(Episodic_split)

within(Episodic_split, max(Episodic_split[,6:9]))
#
# if Epidsodic_split[,6:9] != Ref_Allele or 0, 
# But what if Ref_Allele and alternative = same thing
# Make 
#ref_A <- c(7,8,9)
#if ref = "A", max(Episodic_split[ref_A])

ref_A <- c(7,8,9)
ref_T <- c(6,8,9)
ref_C <- c(6,7,9)
ref_G <- c(7,8,9)

Episodic_split <- within(Episodic_split, {
  Alt_Allele = ifelse (ref == "A", max(Episodic_split[ref_A]), ifelse (ref == "T", max(Episodic_split[ref_T]), ifelse (ref == "C", max(Episodic_split[ref_C]), ifelse (ref == "G", max(Episodic_split[ref_G]), 0))))}) 
