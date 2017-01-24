# install.packages("qqman")
library(qqman)
vignette("qqman")
#str(gwasResults)
getwd()
setwd("/Users/paulknoops/Bioinformatics")
Bowtie_13_24_2 <- read.table("episodic_data_bowtie_2_1-3,2-4.cmh.csv",h=T)
str(Bowtie_13_24_2)
head(Bowtie_13_24_2)
tail(Bowtie_13_24_2)
as.data.frame(table(Bowtie_13_24_2$CHR))
manhattan(Bowtie_13_24_2)
summary(Bowtie_13_24_2$CHR)
#Need to change each chromo to a number
#1 = left, 2 = right (i.e 21 = 2L, 22 = 2R)
# U = 5
# X = 6
#Y = 7
# Het=8

#Trying to rename all that are the same.
#Bowtie_13_24_2$CHR <- c("YHet", "XHet", "X", "Uextra", "U", "4", "3RHet", "3R", "3LHet", "3L", "2RHet", "2R", "2LHet", "2L")
#Order = YHet, 2RHet, 2LHet, 3LHet, 3RHet, U, XHet, 2L, X, 3L, 4, 2R, 3R, Uextra
#When Changes = 78, 228, 218, 318, 328, 5, 68, 21, 6, 31, 4, 22, 32, 58
summary(Bowtie_13_24_2$CHR)
#YHet
#Bowtie_13_24_2$CHR[0:329]
#2RHet: 329+23300
#Bowtie_13_24_2$CHR[330:23629]
#2LHet: 23629+3706
#Bowtie_13_24_2$CHR[23630:27335]
#3LHet: 27335+22376
#Bowtie_13_24_2$CHR[27336:49711]
#3RHet: 49711+25310
#Bowtie_13_24_2$CHR[49712:75021]
#U: 75021+61698
#Bowtie_13_24_2$CHR[75022:136719]
#XHet: 136719+2031
#Bowtie_13_24_2$CHR[136720:138750]
#2L: 138750+301250
#Bowtie_13_24_2$CHR[138751:440000]
#X: 440000+422736
#Bowtie_13_24_2$CHR[440001:862736]
#3L: 862736+336708
#Bowtie_13_24_2$CHR[862737:1199444]
#4: 1199444+18085
#Bowtie_13_24_2$CHR[1199445:1217529]
#2R: 1217529+232060
#Bowtie_13_24_2$CHR[1217530:1449589]
#3R: 1449589+284406
#Bowtie_13_24_2$CHR[1449590:1733995]
#Uextra: 1733995+36692
#Bowtie_13_24_2$CHR[1733996:1770687]
#78
#summary(Bowtie_13_24_2$CHR[0:329])
#228
#summary(Bowtie_13_24_2$CHR[330:23629])
#218
#summary(Bowtie_13_24_2$CHR[23630:27335])
#318
#summary(Bowtie_13_24_2$CHR[27336:49711])
#328
#summary(Bowtie_13_24_2$CHR[49712:75021])
#5
#summary(Bowtie_13_24_2$CHR[75022:136719])
#68
#summary(Bowtie_13_24_2$CHR[136720:138750])
#21
#summary(Bowtie_13_24_2$CHR[138751:440000])
#6
#summary(Bowtie_13_24_2$CHR[440001:862736])
#31
#summary(Bowtie_13_24_2$CHR[862737:1199444])
#4
#summary(Bowtie_13_24_2$CHR[1199445:1217529])
#22
#summary(Bowtie_13_24_2$CHR[1217530:1449589])
#32
#summary(Bowtie_13_24_2$CHR[1449590:1733995])
#Below = 58
#summary(Bowtie_13_24_2$CHR[1733996:1770687])
#summary(Bowtie_13_24_2$CHR)
#names(Bowtie_13_24_2$CHR[1733996:1770687])<- c("58")
#within(df, Bowtie_13_24_2$CHR <- factor(Bowtie_13_24_2$CHR, labels = c(78, 228, 218, 318, 328, 5, 68, 21, 6, 31, 4, 22, 32, 58)))
class(Bowtie_13_24_2)
class(Bowtie_13_24_2$CHR)
#YHet, 2RHet, 2LHet, 3LHet, 3RHet, U, XHet, 2L, X, 3L, 4, 2R, 3R, Uextra
#sapply(Bowtie_13_24_2$CHR,switch,'YHet'=78, '2RHet'=228, '2LHet'=218, '3lHet'=318, '3RHet'=328, 'U'=5, 'XHet'=68, '2L'=21, 'X'=6, '3L'=31, '4'=4, '2R'=22, '3R'=32, 'Uextra'=58)
?manhattan

#Redefine manahattan for chr with symbols; Check saved file with this...

#From forum for code
#manhattan(Bowtie_13_24_2)
#c("YHet", "2RHet", "2LHet", "3LHet", "3RHet", "U", "XHet", "2L", "X", "3L", "4", "2R", "3R", "Uextra")

### Does not work!
#source("https://bioconductor.org/biocLite.R")
#biocLite("quantsmooth")
#library(quantsmooth)
#numericCHR(Bowtie_13_24_2$CHR)
#summary(Bowtie_13_24_2$CHR)


#manhattan(Bowtie_13_24_2, chrlabs = c("YHet", "2RHet", "2LHet", "3LHet", "3RHet", "U", "XHet", "2L", "X", "3L", "4", "2R", "3R", "Uextra"))
summary(Bowtie_13_24_2)
summary(Bowtie_13_24_2$CHR)

#YHet = 78

#Bowtie_13_24_2$CHR[0:329] <- 78
#Bowtie_13_24_2$CHR <- as.character(Bowtie_13_24_2$CHR)
#Bowtie_13_24_2$CHR <- as.vector(Bowtie_13_24_2$CHR)

#replace(Bowtie_13_24_2$CHR, 0:329, 78)

##
#manhattan(subset(Bowtie_13_24_2, CHR == 4))
#manahattan(gwasResults)

###
#source("https://bioconductor.org/biocLite.R")
#biocLite("GWASTools")
#n
library(GWASTools)
#manhattanPlot(Bowtie_13_24_2$P, Bowtie_13_24_2$CHR)
#Sub_bowtie1324_2 <- subset(Bowtie_13_24_2$CHR[138751:1733995])
sub_bow <- Bowtie_13_24_2[138751:1733995, ]
#sub_bow$CHR
#manhattanPlot(sub_bow$P, sub_bow$CHR)
?manhattanPlot
#manhattanPlot(sub_bow$P, sub_bow$CHR, thinThreshold = NULL)
#--> 5e-8 sig level
manhattanPlot(sub_bow$P, sub_bow$CHR, ylim = c(0, (log10(length(sub_bow$P)) +6)))
manhattanPlot(sub_bow$P, sub_bow$CHR, ylim = c(0, (log10(length(sub_bow$P)) +10)))
#length(sub_bow$P)
#log10(1595245)

