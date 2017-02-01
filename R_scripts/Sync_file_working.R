

#To Large
#Data_subset <- read.table("~/Bioinformatics/Sequence_analysis_2016/episodic_data_2R.sync")
# == 3.6 Gb

#Random Region near middle of chromosome == 10000 bp
Data_subset <- read.table("episodic_data_2R_subset.sync")
# == 1.8 MB


head(Data_subset)
colnames(Data_subset) <- c("Chromosome", "Position", "ref", "115ConR1", "115ConR2", "115SelR1", "115SelR2", "38ConR1", "38ConR2", "38SelR1", "38SelR2", "77ConR1", "77ConR2", "77SelR1", "77SelR2", "Base")
head(Data_subset)



#Use tidyr --> make long?
