#Read episodics as frequencies

#packages:
source('episodic_packages.R')
#??haplotype
#Working Dir ==Scripts
episodic_freq <- sync_to_frequencies(file="../R_Data/episodic_data_2R_subset.sync", base.pops=c(rep(FALSE,12), TRUE), header= FALSE)
colnames(episodic_freq) <- c("chr", "pos", "ref","minallele","majallele", 'basePops' ,"ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "AncR1_0") 

summary(episodic_freq)
#Remove any that have all populations minor allele frequency of 0
episodic_freq$sum <- rowSums(episodic_freq[,6:19])
head(episodic_freq)
episodic_freq_sub <- episodic_freq[ -which(episodic_freq$sum<0.05),]

summary(episodic_freq_sub)
dim(episodic_freq_sub)
#171 sites of 10000 that have allele frequncies over 25%

episodic_freq_sub <- subset(episodic_freq_sub, select = -c(sum) )
##Copied

#Why some with N/A
episodic_freq[pos==10635699]

#Remove na with na.omit

dat <- na.omit(episodic_freq_sub)


#summary((dat))
dat_long <- gather(dat, population, min_all_freq, basePops:AncR1_0)
print(dat_long)

#dat_plot <- ggplot(dat_long, aes(x= pos, y = min_all_freq, colour = population))
#dat_plot2 <- dat_plot+geom_point() + ggtitle("all") 
#+ guides(colour=FALSE)
#print(dat_plot2)

#Change per position:
#Group by position before making long: take the minor allele frequency and subtract the base == Change

ddiff <- dat %>%
  group_by(pos) %>%
  summarise(ref=ref, minallele=minallele, majallele=majallele, ConR1_115=ConR1_115-AncR1_0, ConR2_115=ConR2_115-AncR1_0, SelR1_115=SelR1_115-AncR1_0, SelR2_115=SelR2_115-AncR1_0, ConR1_38=ConR1_38-AncR1_0, ConR2_38=ConR2_38-AncR1_0, SelR1_38=SelR1_38-AncR1_0, SelR2_38=SelR2_38-AncR1_0, ConR1_77=ConR1_77-AncR1_0, ConR2_77=ConR2_77-AncR1_0, SelR1_77=SelR1_77-AncR1_0, SelR2_77=SelR2_77-AncR1_0, AncR1_0=AncR1_0-AncR1_0)

#summary(ddiff)
#Make long
ddiff_long <- gather(ddiff, population_diff, minallele_freq, ConR1_115:AncR1_0)

#Removes 0 values-- plots have a long strip at 0 b/c many alleles are not changing.
ddiff_long <- ddiff_long[ -which(ddiff_long$minallele_freq==0),]


#head(ddiff_long)
#ddiff_long$pos <- as.character(ddiff_long$pos)
#plot of every position, and the change from the base generation (0 == base generation)
ddiff_plot <- ggplot(ddiff_long, aes(x= pos, y = minallele_freq, colour = population_diff))
ddiff_plot2 <- ddiff_plot+geom_point() + ggtitle("All sequences change from base generation") 
#+ guides(colour=FALSE)
#print(ddiff_plot2)


#Split population into Treatment and Generation
ddiff_long2 <- ddiff_long %>% 
  separate(population_diff, c('Treatment', "Generation"), "_")


#head(ddiff_long2)

#split by generation: Plot each generation in own plot by postion

#Spliting each generation:
ddiff_38 <- ddiff_long2[ which(ddiff_long2$Generation==38),]
ddiff_77 <- ddiff_long2[ which(ddiff_long2$Generation==77),]
ddiff_115 <- ddiff_long2[ which(ddiff_long2$Generation==115),]

#head(ddiff_115)
plot_115 <- ggplot(ddiff_115, aes(x= pos, y = minallele_freq, colour = Treatment))
plot2_115 <- plot_115+geom_point() + ggtitle("delta115")
#print(plot2_115)


#head(ddiff_38)
plot_38 <- ggplot(ddiff_38, aes(x= pos, y = minallele_freq, colour = Treatment))
plot2_38 <- plot_38+geom_point() + ggtitle("dalta38")
#print(plot2_38)

#head(ddiff_77)
plot_77 <- ggplot(ddiff_77, aes(x= pos, y = minallele_freq, colour = Treatment))
plot2_77 <- plot_77+geom_point() + ggtitle("delta77")
#print(plot2_77)

#Combine plots together to show in one window:


#multiplot(ddiff_plot2, plot2_77, plot2_38, plot2_115,cols=2)


#Split by Treatment:
ddiff_long3 <- ddiff_long2 %>% 
  separate(Treatment, c('Treatment', "Replicate"), "R")
ddiff_sel <- ddiff_long3[ which(ddiff_long3$Treatment=='Sel'),]
#print(head(ddiff_sel))
plot_sel <- ggplot(ddiff_sel, aes(x= pos, y = minallele_freq, colour = Generation))
plot2_sel <- plot_sel+geom_point() + ggtitle("Selection")
#print(plot2_sel)

ddiff_con <- ddiff_long3[ which(ddiff_long3$Treatment=='Con'),]
#print(head(ddiff_con))
plot_con <- ggplot(ddiff_con, aes(x= pos, y = minallele_freq, colour = Generation))
plot2_con <- plot_con+geom_point() + ggtitle("Control")
#print(plot2_con)

#multiplot(plot2_sel, plot2_con, cols =2)




###### Changes from previous generation
genDiff <- dat %>%
  group_by(pos) %>%
  summarise(ref=ref, minallele=minallele, majallele=majallele, ConR1_115=ConR1_115-ConR1_77, ConR2_115=ConR2_115-ConR2_77, SelR1_115=SelR1_115-SelR1_77, SelR2_115=SelR2_115-SelR2_77, ConR1_38=ConR1_38-AncR1_0, ConR2_38=ConR2_38-AncR1_0, SelR1_38=SelR1_38-AncR1_0, SelR2_38=SelR2_38-AncR1_0, ConR1_77=ConR1_77-ConR1_38, ConR2_77=ConR2_77-ConR2_38, SelR1_77=SelR1_77-SelR1_38, SelR2_77=SelR2_77-SelR2_38)

genDiff_long <- gather(genDiff, population_diff, minallele_freq, ConR1_115:SelR2_77)

#Removes 0 values 
genDiff_long <- genDiff_long[ -which(genDiff_long$minallele_freq==0),]

#head(genDiff_long)
gendiff_plot <- ggplot(genDiff_long, aes(x= pos, y = minallele_freq, colour = population_diff))
gendiff_plot2 <- gendiff_plot+geom_point() + ggtitle("Generation Differences; 115-77, 77-38, 38-base")
#print(gendiff_plot2)

#gendiff_long2
gendiff_long2 <- genDiff_long %>% 
  separate(population_diff, c('Treatment', "Generation"), "_")


#Needed to order 
gendiff_long2$Generation <- as.numeric(gendiff_long2$Generation)
gendiff_long2 <- gendiff_long2[order(gendiff_long2$Generation),]
gendiff_long2$Generation <- as.factor(gendiff_long2$Generation)



gendiff_plot55 <- ggplot(gendiff_long2, aes(x=Generation, y=minallele_freq, colour=Treatment))
gendiff_plot56 <- gendiff_plot55 + geom_jitter() + ggtitle('Generational change in allele frequency') 
#+ guides(colour=T)
#print(gendiff_plot56)

gendiff_plot44 <- ggplot(gendiff_long2, aes(x=Treatment, y=minallele_freq, colour=Generation))
gendiff_plot45 <- gendiff_plot44 + geom_jitter() + ggtitle('Generation Spread')
#+ guides(colour=T)
#print(gendiff_plot45)

#gendiff_long2$Generation <- as.numeric(gendiff_long2$Generation)
#gendiff_plot66 <- ggplot(gendiff_long2, aes(x=Generation, y=minallele_freq, colour=pos))
#gendiff_plot65 <- gendiff_plot66 + geom_point() + ggtitle('Generation') + geom_
#+ guides(colour=T)
#print(gendiff_plot65)
#dat



##Correlations:

#with(DGRP_sub, cor(prop_spider[Sex == "Female"],prop_spider[Sex == "Male"] ))
#with(DGRP_sub, plot(x = prop_spider[Sex == "Female"],y = prop_spider[Sex == "Male"], xlab = "Females", ylab = "Males", abline(line2), main = "Male Female correlation: proportion with spider"))
dat_long <- dat_long[ -which(dat_long$population=="basePops"),]

dat_long2 <- dat_long %>% 
  separate(population, c('Treatment', "Generation"), "_")
head(dat_long2)
with(dat_long2, cor(min_all_freq[Treatment =="SelR1"], min_all_freq[Treatment=="SelR2"]))
#with(dat_long2, plot(x=min_all_freq[Treatment =="SelR1"],y= min_all_freq[Treatment=="SelR2"]))

#head(dat_long)
#par(mfrow=c(3,2))
#with(dat_long, plot(x=min_all_freq[population=="SelR1_115"],y= min_all_freq[population=="SelR2_77"], xlab = "SelR1_115", ylab= "SelR2_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
#with(dat_long, plot(x=min_all_freq[population=="SelR2_115"],y= min_all_freq[population=="SelR2_77"], xlab = "SelR2_115", ylab= "SelR2_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
#with(dat_long, plot(x=min_all_freq[population=="SelR1_115"],y= min_all_freq[population=="SelR1_77"], xlab = "SelR1_115", ylab= "SelR1_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))

#with(dat_long, plot(x=min_all_freq[population=="SelR1_77"],y= min_all_freq[population=="SelR2_38"], xlab = "SelR1_77", ylab= "SelR2_38",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
#with(dat_long, plot(x=min_all_freq[population=="SelR2_77"],y= min_all_freq[population=="SelR2_38"], xlab = "SelR2_77", ylab= "SelR2_38",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
#with(dat_long, plot(x=min_all_freq[population=="SelR1_77"],y= min_all_freq[population=="SelR1_38"], xlab = "SelR1_77", ylab= "SelR1_38",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))



#par(mfrow=c(3,1))
#with(ddiff_115,plot(x=minallele_freq[Treatment=="SelR1"], y= minallele_freq[Treatment=="SelR2"], main = "Generation 115 Selection Lineages"))
#with(ddiff_77,plot(x=minallele_freq[Treatment=="SelR1"], y= minallele_freq[Treatment=="SelR2"], main = "Generation 77 Selection Lineages"))
#with(ddiff_38,plot(x=minallele_freq[Treatment=="SelR1"], y= minallele_freq[Treatment=="SelR2"], main = "Generation 38 Selection Lineages"))


par(mfrow=c(1,1))


d_plot <- ggplot(dat_long, aes(x= pos, y = min_all_freq, colour = population))
d_plot2 <- d_plot+geom_point() + ggtitle("all") 
#print(d_plot2)

#Spliting each generation:
d_38 <- dat_long2[ which(dat_long2$Generation==38),]
d_77 <- dat_long2[ which(dat_long2$Generation==77),]
d_115 <- dat_long2[ which(dat_long2$Generation==115),]

#head(ddiff_115)
pl_115 <- ggplot(d_115, aes(x= pos, y = min_all_freq, colour = Treatment))
pl2_115 <- pl_115+geom_point() + ggtitle(115)
#print(pl2_115)


#head(ddiff_38)
pl_38 <- ggplot(d_38, aes(x= pos, y = min_all_freq, colour = Treatment))
pl2_38 <- pl_38+geom_point() + ggtitle(38)
#print(pl2_38)

#head(ddiff_77)
pl_77 <- ggplot(d_77, aes(x= pos, y = min_all_freq, colour = Treatment))
pl2_77 <- pl_77+geom_point() + ggtitle(77)
#print(pl2_77)

#Combine plots together to show in one window:
#multiplot(d_plot2, pl2_77, pl2_38, pl2_115,cols=2)


#Plots for analysis
plot_sel_2 <- ggplot(ddiff_sel, aes(x= pos, y = minallele_freq, colour = Generation, shape = Replicate))
plot2_sel_2 <- plot_sel_2+geom_point() + ggtitle("Selection difference from base generation")

#Plots for analysis:

#Selection
print(plot2_sel_2)

print(gendiff_plot56)

#Plot of allele frequencies
multiplot(d_plot2, pl2_77, pl2_38, pl2_115,cols=2)

#Plot of differences from base generation
multiplot(ddiff_plot2, plot2_77, plot2_38, plot2_115,cols=2)

#Correlation Plots
par(mfrow=c(3,2))
with(dat_long, plot(x=min_all_freq[population=="SelR1_115"],y= min_all_freq[population=="SelR2_77"], xlab = "SelR1_115", ylab= "SelR2_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
with(dat_long, plot(x=min_all_freq[population=="SelR2_115"],y= min_all_freq[population=="SelR2_77"], xlab = "SelR2_115", ylab= "SelR2_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
with(dat_long, plot(x=min_all_freq[population=="SelR1_115"],y= min_all_freq[population=="SelR1_77"], xlab = "SelR1_115", ylab= "SelR1_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))

with(dat_long, plot(x=min_all_freq[population=="SelR1_77"],y= min_all_freq[population=="SelR2_38"], xlab = "SelR1_77", ylab= "SelR2_38",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
with(dat_long, plot(x=min_all_freq[population=="SelR2_77"],y= min_all_freq[population=="SelR2_38"], xlab = "SelR2_77", ylab= "SelR2_38",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
with(dat_long, plot(x=min_all_freq[population=="SelR1_77"],y= min_all_freq[population=="SelR1_38"], xlab = "SelR1_77", ylab= "SelR1_38",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))

par(mfrow=c(2,2))

with(dat_long, plot(x=min_all_freq[population=="SelR2_115"],y= min_all_freq[population=="SelR1_77"], xlab = "SelR2_115", ylab= "SelR1_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
with(dat_long, plot(x=min_all_freq[population=="SelR1_115"],y= min_all_freq[population=="SelR2_77"], xlab = "SelR1_115", ylab= "SelR2_77",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))

par(mfrow=c(3,1))
with(d_115,plot(x=min_all_freq[Treatment=="SelR1"], y= min_all_freq[Treatment=="SelR2"], main = "Generation 115 Selection Lineages"))
with(d_77,plot(x=min_all_freq[Treatment=="SelR1"], y= min_all_freq[Treatment=="SelR2"], main = "Generation 77 Selection Lineages"))
with(d_38,plot(x=min_all_freq[Treatment=="SelR1"], y= min_all_freq[Treatment=="SelR2"], main = "Generation 38 Selection Lineages"))



#par(mfrow=c(3,1))
#with(d_115,plot(x=min_all_freq[Treatment=="SelR1"], y= min_all_freq[Treatment=="SelR2"], main = "Generation 115 Selection Lineages minor allele frequencies"))
#with(d_77,plot(x=min_all_freq[Treatment=="SelR1"], y= min_all_freq[Treatment=="SelR2"], main = "Generation 77 Selection Lineages minor allele frequencies"))
#with(d_38,plot(x=min_all_freq[Treatment=="SelR1"], y= min_all_freq[Treatment=="SelR2"], main = "Generation 38 Selection Lineages minor allele frequencies"))
