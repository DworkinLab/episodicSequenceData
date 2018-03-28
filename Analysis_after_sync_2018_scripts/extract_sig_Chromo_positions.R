# Extract positions of interest:
## Assuming a list or data frame is available with all positions:

require('data.table')

name.Columns <- c("Chromosome", "Position", "ref", 
	"ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", 
	"ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", 
	"ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", 
	"SelR1_0")

### X Chromosome:
## Read in data: 

	episodic_data <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_X.sync')
	colnames(episodic_data) <- name.Columns

## Keep rows with Position in list (position == column 2)

	siglist <- fread('/home/paul/episodicData/positions/X_positions.csv', h=T)
	siglist_2 <- as.list(siglist$x)
	episodic_counts <- episodic_data[episodic_data$Position %in% siglist_2 ,]
	
## write csv file with only positions
	
	write.csv(episodic_counts, file="/home/paul/episodicData/positions/X_positions.sync")

## remove:
	rm(episodic_data)
	rm(siglist)
	rm(episodic_counts)
	
### 2R Chromosome:
	episodic_data <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_2R.sync')
	colnames(episodic_data) <- name.Columns
	siglist <- fread('/home/paul/episodicData/positions/2R_positions.csv', h=T)
	siglist_2 <- as.list(siglist$x)
	episodic_counts <- episodic_data[episodic_data$Position %in% siglist_2 ,]
	write.csv(episodic_counts, file="/home/paul/episodicData/positions/2R_positions.sync")
	rm(episodic_data)
	rm(siglist)
	rm(episodic_counts)
	
### 2L Chromosome:
	episodic_data <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_2L.sync')
	colnames(episodic_data) <- name.Columns
	siglist <- fread('/home/paul/episodicData/positions/2L_positions.csv', h=T)
	siglist_2 <- as.list(siglist$x)
	episodic_counts <- episodic_data[episodic_data$Position %in% siglist_2 ,]
	write.csv(episodic_counts, file="/home/paul/episodicData/positions/2L_positions.sync")
	
	rm(episodic_data)
	rm(siglist)
	rm(episodic_counts)
	
### 3R Chromosome:
	episodic_data <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_3R.sync')
	colnames(episodic_data) <- name.Columns
	siglist <- fread('/home/paul/episodicData/positions/3R_positions.csv', h=T)
	siglist_2 <- as.list(siglist$x)
	episodic_counts <- episodic_data[episodic_data$Position %in% siglist_2 ,]
	write.csv(episodic_counts, file="/home/paul/episodicData/positions/3R_positions.sync")
	rm(episodic_data)
	rm(siglist)
	rm(episodic_counts)
	
### 3L Chromosome:
	episodic_data <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_3L.sync')
	colnames(episodic_data) <- name.Columns
	siglist <- fread('/home/paul/episodicData/positions/3L_positions.csv', h=T)
	siglist_2 <- as.list(siglist$x)
	episodic_counts <- episodic_data[episodic_data$Position %in% siglist_2 ,]
	write.csv(episodic_counts, file="/home/paul/episodicData/positions/3L_positions.sync")
	rm(episodic_data)
	rm(siglist)
	rm(episodic_counts)
	
### 4 Chromosome:
	episodic_data <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_4.sync')
	colnames(episodic_data) <- name.Columns
	siglist <- fread('/home/paul/episodicData/positions/4_positions.csv', h=T)
	siglist_2 <- as.list(siglist$x)
	episodic_counts <- episodic_data[episodic_data$Position %in% siglist_2 ,]
	write.csv(episodic_counts, file="/home/paul/episodicData/positions/4_positions.sync")
	rm(episodic_data)
	rm(siglist)
	rm(episodic_counts)
