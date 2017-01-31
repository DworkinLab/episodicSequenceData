

#To Large
#Data_subset <- read.table("~/Bioinformatics/Sequence_analysis_2016/episodic_data_2R.sync")
# == 3.6 Gb

Data_subset <- read.table("episodic_data_Ian_subset.sync")
# == 194 KB
#old data (incorrect base generation?)
# Same method into Git_journal to get a random region of one chromosome.....

head(Data_subset)
colnames(Data_subset) <- c("Chromosome", "Position", "ref", "115ConR1", "115ConR2", "115SelR1", "115SelR2", "38ConR1", "38ConR2", "38SelR1", "38SelR2", "77ConR1", "77ConR2", "77SelR1", "77SelR2", "MDG2", "MGD")
head(Data_subset)
#Use tidyr --> make long?
