source('episodic_packages.R')
#??haplotype
#Working Dir ==Scripts
episodic_freq <- sync_to_frequencies(file="../R_Data/episodic_data_2R_subset.sync", base.pops=c(rep(FALSE,12), TRUE), header= FALSE)
colnames(episodic_freq) <- c("chr", "pos", "ref","minallele","majallele", 'basePops' ,"ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "AncR1_0") 


head(episodic_freq)

episodic_freq$max <- do.call(pmax, episodic_freq[,6:19])
episodic_freq_sub <- episodic_freq[ -which(episodic_freq$max<0.05),]
episodic_freq_sub <- subset(episodic_freq_sub, select = -c(max) )


head(episodic_freq_sub)
dim(episodic_freq_sub)

#ep_long <- gather(episodic_freq, population, min_all_freq, basePops:AncR1_0)
#ep_long <- ep_long[ -which(ep_long$population=="basePops"),]
#head(ep_long)
#ep_long2 <- ep_long[ which(ep_long$min_all_freq>0.9),]
#ep_long2
