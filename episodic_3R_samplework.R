setwd('/Users/paulknoops/Sequence_analysis_2016')
ThirdR<- read.table('episodic_data_Ian_subset.sync', h=F)
summary(ThirdR)
Pop_13_24 <- ThirdR[,1:7]
summary(Pop_13_24)
head(Pop_13_24)



