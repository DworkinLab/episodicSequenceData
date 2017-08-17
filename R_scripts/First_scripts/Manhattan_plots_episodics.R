# Using information that worked from the episodics_trial_Manhattan.R

#Using GWASTools package
#source("https://bioconductor.org/biocLite.R")
#biocLite("GWASTools")
#n
library(GWASTools)

#For help on Manhattan plot
?manhattanPlot

#Can change to locations of files (i.e Sequence_Analysis_2016)
getwd()
#setwd("/Users/paulknoops/Bioinformatics/Sequence_Analysis_2016/Data")
#Give each a unique name to work with and bring in Data
Bowtie_13_24_2 <- read.table("episodic_data_bowtie_2_1-3,2-4.cmh.csv",h=T)
# Can view data
str(Bowtie_13_24_2)
head(Bowtie_13_24_2)
tail(Bowtie_13_24_2)
as.data.frame(table(Bowtie_13_24_2$CHR))

#Subset Data to only include 2L, X, 3L, 4, 2R, and 3R
# Use information on lengths to determin
summary(Bowtie_13_24_2$CHR)
#Add all numbers up to 2L (order = "YHet", "2RHet", "2LHet", "3LHet", "3RHet", "U", "XHet", "2L", "X", "3L", "4", "2R", "3R", "Uextra"
#In this case:
#3706+23300+22376+25310+61698+2031+329

#Add all through to 3R (include value above)
#138750+301250+232060+336708+284406+18085+422736

#The Range of first addition (+1) and second will give range of 2L -- 3R
#Create the subset of these and name as a unique name (similar as first)
sub_bow <- Bowtie_13_24_2[138751:1733995, ]
summary(sub_bow)
#Create Manhattan Plot:
manhattanPlot(sub_bow$P, sub_bow$CHR, main = "bowtie", ylim = c(0, (log10(length(sub_bow$P)) +15)))

#Same for BWA
BWA_13_24_2 <- read.table("episodic_data_2_1-3,2-4.cmh.gwas", h=T)
summary(BWA_13_24_2$CHR)
#2636+22909+20609+25047+81515+1631+743
#155090
#155090+171613+129285+196201+156692+12966+314762
#1136609
#Alterntive = lenght(BWA_13_24_2$CHR) -Uextra
#length(BWA_13_24_2$CHR)
#1236435-99826
sub_BWA <- BWA_13_24_2[155091:1136609, ]
head(sub_BWA)
head(sub_bow)
length(sub_BWA$P)
length(sub_bow$P)

#create plot:
manhattanPlot(sub_BWA$P, sub_BWA$CHR, main = "bwa", signif = 5e-8, ylim = c(0, (log10(length(sub_BWA$P)) +35)))
#Set a comparison of SNP to P-value; if different (by say a factor of .1??), remove
#may not want to remove... but should be neglected to base line..?

## Combining into one data frame?
#sub_BWA$Test <- ifelse(sub_BWA$P <= 100000000, "A", "B")
#sub_bow$Test <- ifelse(sub_bow$P <= 1000000000000, "B", "C")
#combined <- rbind(sub_BWA, sub_bow)
#summary(combined)
#manhattanPlot(combined$P, combined$CHR, main = "combined", ylim = c(0, log10(length(combined$P)) + 35))
