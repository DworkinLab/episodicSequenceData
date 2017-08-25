#Combine bwa and Bowtie Together and find p-values:
#Combine the 4th chromosome

require(ggplot2)
require(dplyr)
library(GWASTools)
#install.packages("manhattanly")
library(manhattanly)



#First make into one data frame per chromosome (will need with bigger ones split) and add a signifyer to mapper:

#Change the working Directory per Chromosome:

#BWA:
setwd('/Users/paulknoops/Bioinformatics/episodic_practice/Chromo4Coeffs/bwa')
mycsvs <- list.files(pattern=".coeffs.csv")
coeffs_bwa <- NULL
for (file in mycsvs){
  print(file)
  coeffs1 <- read.csv(file, h=T)
  coeffs1$mapper <- "bwa"
  coeffs_bwa <- rbind(coeffs_bwa, coeffs1)
}


#Change the working Directory per Chromosome:

#BOWITE:

setwd('/Users/paulknoops/Bioinformatics/episodic_practice/Chromo4Coeffs/bowtie')
mycsvs <- list.files(pattern=".coeffs.csv")
coeffs_bowtie <- NULL
for (file in mycsvs){
  print(file)
  coeffs2 <- read.csv(file, h=T)
  coeffs2$mapper <- "bowtie"
  coeffs_bowtie <- rbind(coeffs_bowtie, coeffs2)
}


#combine to one:
X <- rbind(coeffs_bowtie, coeffs_bwa)
head(X)

#Split based on chosed Effect: may want to write a csv per effect?

X2 <- X[which(X$Effects=="TreatmentSel"),]
#X2 <- X[which(X$Effects=="TreatmentSel:Generation"),]
#X2 <- X[which(X$Effects=="Generation"),]
#X2 <- X[which(X$Effects=="Interept"),]

#Title for ggplot based on effect
Title <- as.character(X2$Effects[1])


#Fine positions that are mapped by both bwa and bowtie
DF <- X2 %>%
  group_by(position) %>%
  summarise(mapper_1 = mapper[1], mapper_2=mapper[2])

DF2 <- DF[which(DF$mapper_2=="bwa" & DF$mapper_1=="bowtie"),]

#remove those from the X2 (merged and effect) not with both mappers
coeffs_treatment <- X2[(X2$position %in% DF2$position),]

#Can likely remove all coefficients other than p-values

#Using the mean values between the two mappers: note mean(log_p) != -log10(mean(p-value))

mean_coeffs <- coeffs_treatment %>%
  group_by(position, chr, Effects) %>%
  summarise(Estimate = mean(Estimate), SE = mean(Standard_error), Z = mean(z.value), pvalue = mean(p.value), meanlog_p = mean(log_p), log_p = -log10(mean(p.value)))

#Plot to ggplot:

ggmean <- ggplot(data = mean_coeffs, aes(x=position, y=log_p))
ggmean2 <- ggmean + geom_point(size = 0.5) + geom_hline(yintercept = -log10(0.05)) +
  ggtitle(Title) + 
  ylab("-log10(p-value)")

print(ggmean2)

#Plot with manhattan plot (less manipulability)
manhattanPlot(mean_coeffs$pvalue, mean_coeffs$chr, ylim= c(0,2.25), signif = 0.05)

#Plot with manhattanly: nice plots but more useful as an interctive plot (which not useful here)
mean_coeffs2 <- mean_coeffs[c("position", "chr", "pvalue")]
colnames(mean_coeffs2) <- c("BP", "CHR", "P")
manhattanly(mean_coeffs2)



#May want to do smallest value of the two (recomended form Kofler et al. 2016).

min_coeffs <- coeffs_treatment %>%
  group_by(position, chr, Effects) %>%
  summarise(Estimate = min(Estimate), SE = min(Standard_error), Z = min(z.value), pvalue = min(p.value), minlog_p = min(log_p), log_p = -log10(min(p.value)))

ggmin <- ggplot(data = min_coeffs, aes(x=position, y=log_p))
ggmin2 <- ggmin + geom_point(size = 0.5) + geom_hline(yintercept = -log10(0.05)) +
  ggtitle(Title)
print(ggmin2)


manhattanPlot(min_coeffs$pvalue, min_coeffs$chr, ylim= c(0,4), signif = 0.05)


#Just to see them side by side
head(coeffs_treatment)
manhattanPlot(coeffs_treatment$p.value, coeffs_treatment$mapper, ylim= c(0,4), signif = 0.05)


