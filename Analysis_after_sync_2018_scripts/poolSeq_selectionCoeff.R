### Running the PoolSeq Package to run on .sync files for estimates of selection coefficeints per position
### Requires: R (>= 3.3.1), data.table (>= 1.9.4), foreach (>= 1.4.2), stringi (>= 0.4-1), matrixStats (>= 0.14.2)

  args <- commandArgs(trailingOnly = TRUE)
  
### Required Packages:

  #install.packages("/home/paul/poolSeq_0.3.2.tar.gz", repos=NULL, type="source")
  #install.packages("/home/paul/matrixStats_0.53.0.tar.gz", repos=NULL, type="source")

### Not available: so source seperate:
  #require(poolSeq)
  
### These are part of the dependencies for poolSeq
  
  require(methods)
  require(data.table)
  require(foreach)
  require(stringi)
  require(matrixStats)
  
### Source the scripts (Copied) for Pool-Seq (only one fails and is not needed)
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/loadaf.R')  
  #estne.R Fails.
  #source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estne.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estsh.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/idsel.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/simaf.R')

### Possibly need custom function to read in manipulated .sync files:
	### Needed for manipulated .sync files (one basic change labeled at top of changed script:
  
  source("/home/paul/episodicData/novoalign/novo_Rscripts/Taus_ReadSync.R")

### Read in the data file for args[1]

  setwd(args[3])
  
  mySync <- read.sync_Personal(file=args[1], gen=c(115, 115, 38, 38, 77, 77, 0, 0), repl=c(1,2,1,2,1,2,1,2), polarization = "rising")

### Make data.frame of just alleles information to sort out relevent positions:

  ff <- as.data.frame(mySync@alleles)
  pst <- as.numeric(ff$pos)
  pst2 <- sort(pst)
  rm(pst)
  rm(ff)

### Create empty data frame to read into for estiamting S:
  
  DF <- data.frame(NULL)
  ccc <- c(0,38,77,115)

  for (i in pst2) {
  	Traj115 <- af.traj(mySync, args[2], repl=c(1,2), pos=i)
  	Bfsf <- estimateSH(Traj115, Ne=150, t=ccc, h=0.5, haploid = FALSE, simulate.p.value=TRUE)
  	Fd <- data.frame(Bfsf$s, Bfsf$p0, Bfsf$p.value)
 	 Fd$pos <- i
  	DF <- rbind(DF, Fd)
  	DF <- na.omit(DF)
  	#print(paste("Running entity:", i, "of", END))
  	rm(i)
	
	}
 
  x2 <- args[1]
  x3 <- gsub("\\..*","", x2)
  write.csv(DF, file=paste(args[3], "/", x3, ".csv", sep=""), row.names=FALSE)
  
  rm(DF)
  rm(mySync)
  rm(ccc)
  rm(pst2)
