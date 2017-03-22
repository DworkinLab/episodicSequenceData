
#packages:
source('episodic_packages.R')
#??haplotype
#Working Dir ==Scripts
episodic_freq <- sync_to_frequencies(file="../R_Data/episodic_data_2R_subset.sync", base.pops=c(rep(FALSE,12), TRUE), header= FALSE)
colnames(episodic_freq) <- c("chr", "pos", "ref","minallele","majallele", 'basePops' ,"ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "AncR1_0") 

#Remove those under certain threshold: 25% minor allele frequency in sum of sequences
episodic_freq$sum <- rowSums(episodic_freq[,6:19])
episodic_freq_sub <- episodic_freq[ -which(episodic_freq$sum<0.25),]
episodic_freq_sub <- subset(episodic_freq_sub, select = -c(sum) )
episodic_freq_sub <- na.omit(episodic_freq_sub)
head(episodic_freq_sub)
# filter replicated time series data for informative SNPs
#dat_filtered=initialize_SNP_time_series(chr=ex_dat$chr, pos=ex_dat$pos, base.freq=ex_dat$basePops, lib.freqs=ex_dat[,7:ncol(ex_dat), with=FALSE], pop.ident=rep(1:5,each=4), pop.generation=rep(c(0:3)*20,times = 5), use.libs=rep(TRUE,20))

episodic_filtered <- with(episodic_freq, initialize_SNP_time_series(chr = chr, pos= pos, base.freq = basePops, lib.freqs = episodic_freq[,7:19], pop.ident = c(1,2,1,2,1,2,1,2,1,2,1,2,1), pop.generation = c(115,115,115,115,38,38,38,38,77,77,77,77,0), use.libs = rep(TRUE,13)))
