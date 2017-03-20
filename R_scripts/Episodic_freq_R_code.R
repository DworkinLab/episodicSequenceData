#Read episodics as frequencies
#install.packages('dplyr')
library(dplyr)

#install.packages('tidyr')
library(tidyr)

#install.packages('ggplot2')
library(ggplot2)

#install.packages('haploReconstruct')
library(haploReconstruct)

#Multiplot:
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#Working Dir ==Scripts
episodic_freq <- sync_to_frequencies(file="../R_Data/episodic_data_2R_subset.sync", base.pops=c(rep(FALSE,12), TRUE), header= FALSE)
colnames(episodic_freq) <- c("chr", "pos", "ref","minallele","majallele", 'basePops' ,"ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "AncR1_0") 

#Remove any that have all populations minor allele frequency of 0
episodic_freq$sum <- rowSums(episodic_freq[,6:19])
head(episodic_freq)
episodic_freq_sub <- episodic_freq[ -which(episodic_freq$sum==0),]

head(episodic_freq_sub)

episodic_freq_sub <- subset(episodic_freq_sub, select = -c(sum) )
##Copied


dat <- episodic_freq_sub
print(head(dat))
dat_long <- gather(dat, population, min_all_freq, basePops:AncR1_0)
print(dat_long)


#dat_long



#Hashed out all below: fixed and made neat in epidosic_plots.R 

#Long based on each population average minor allele frequency:

#dat2 <- dat_long %>%
#  group_by(population) %>%
#  summarise(avg_allele=mean(min_all_freq))

#dat2 <- dat2[-c(2),]
#dat2
#dat2 <- dat2 %>% 
#  separate(population, c('Treatment', "Generation"), "_")

#dat2$Generation <- as.numeric(dat2$Generation)
#dat2 <- dat2[order(dat2$Generation),]
#dat2$Generation <- as.character(dat2$Generation)

#dat2 <- as.data.frame(dat2)
#Plot of average minor allele frequencies across full genome:
#g1 <- ggplot(dat2, aes(x=Generation, y=avg_allele, colour=Treatment))
#g2 <- g1+geom_point(size = 5)
#g2

#Long with grouping population and the minor allele, for the average minor allele frequencies in each population based on the allele.

#dat3 <- dat_long %>%
#  group_by(population, minallele) %>%
#  summarise(avg_allele=mean(min_all_freq))
#dat3
#dat3 <- dat3[-c(5:8),]

#dat3 <- dat3 %>% 
#  separate(population, c('Treatment', "Generation"), "_")

#g3 <- ggplot(dat3, aes(x=population, y = avg_allele, colour = minallele))
#g4 <- g3 + geom_point(size = 5)
#g4

#Change from base generation -- average for all (from dat2) -- most imformative -- better version in episodic_plots.R
#dat2$change <- dat2$avg_allele - dat2[1,3]
#dat2
#g5 <- ggplot(dat2, aes(x=Generation, y = change, colour = Treatment))
#g6 <- g5 + geom_point(size = 8)
#g6 + geom_line()

#Per base change; using dat3, change for each minor allele from the base pop.


#dat3 <- as.data.frame(dat3)
#dat3$change <- 0
#dat3
#dat3$change <- ifelse(dat3$minallele=="A", dat3$avg_allele - dat3[1,3], ifelse(dat3$minallele=="C", dat3$avg_allele - dat3[2,3], ifelse(dat3$minallele=="G", dat3$avg_allele - dat3[3,3],ifelse(dat3$minallele=="T", dat3$avg_allele - dat3[4,3], dat3$change))))


#Split up the population into Treatment and Generation for plot of minor allele change at each generation:
#BAD
#dat4 <- dat3 %>% 
#  separate(population, c('Treatment', "Generation"), "_")


#dat4$Generation <- as.numeric(dat4$Generation)
#g7 <- ggplot(dat4, aes(x=Generation, y = change, colour = minallele))
#g8 <- g7 + geom_jitter(size = 8)
#g8






#Plot of each position and the accompaning minor allele frequency for all populations: Better version in episodic_plots.R
#dat_long$pos <- as.character(dat_long$pos)
#p1 <- ggplot(dat_long, aes(x= pos, y = min_all_freq, colour = population))
#p2 <- p1+geom_point()
#print(p2)








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
#plot of every position, and the change from the base generation (0 == base generation)
ddiff_plot <- ggplot(ddiff_long, aes(x= pos, y = minallele_freq, colour = population_diff))
ddiff_plot2 <- ddiff_plot+geom_point() + ggtitle("all") 
#+ guides(colour=FALSE)
print(ddiff_plot2)


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
plot2_115 <- plot_115+geom_point() + ggtitle(115)
print(plot2_115)


#head(ddiff_38)
plot_38 <- ggplot(ddiff_38, aes(x= pos, y = minallele_freq, colour = Treatment))
plot2_38 <- plot_38+geom_point() + ggtitle(38)
print(plot2_38)

#head(ddiff_77)
plot_77 <- ggplot(ddiff_77, aes(x= pos, y = minallele_freq, colour = Treatment))
plot2_77 <- plot_77+geom_point() + ggtitle(77)
print(plot2_77)

#Combine plots together to show in one window:
multiplot(ddiff_plot2, plot2_77, plot2_38, plot2_115,cols=2)


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

multiplot(plot2_sel, plot2_con, cols =2)
#Group by position== average per replicate?




###### Changes from previous generation
genDiff <- dat %>%
  group_by(pos) %>%
  summarise(ref=ref, minallele=minallele, majallele=majallele, ConR1_115=ConR1_115-ConR1_77, ConR2_115=ConR2_115-ConR2_77, SelR1_115=SelR1_115-SelR1_77, SelR2_115=SelR2_115-SelR2_77, ConR1_38=ConR1_38-AncR1_0, ConR2_38=ConR2_38-AncR1_0, SelR1_38=SelR1_38-AncR1_0, SelR2_38=SelR2_38-AncR1_0, ConR1_77=ConR1_77-ConR1_38, ConR2_77=ConR2_77-ConR2_38, SelR1_77=SelR1_77-SelR1_38, SelR2_77=SelR2_77-SelR2_38, AncR1_0=AncR1_0-AncR1_0)

genDiff_long <- gather(genDiff, population_diff, minallele_freq, ConR1_115:AncR1_0)

#Removes 0 values 
genDiff_long <- genDiff_long[ -which(genDiff_long$minallele_freq==0),]

#head(genDiff_long)
gendiff_plot <- ggplot(genDiff_long, aes(x= pos, y = minallele_freq, colour = population_diff))
gendiff_plot2 <- gendiff_plot+geom_point() + ggtitle("Generation Differences; 115-77, 77-38, 38-base")
print(gendiff_plot2)

#gendiff_long2
gendiff_long2 <- genDiff_long %>% 
  separate(population_diff, c('Treatment', "Generation"), "_")


#Needed to order 
gendiff_long2$Generation <- as.numeric(gendiff_long2$Generation)
gendiff_long2 <- gendiff_long2[order(gendiff_long2$Generation),]
gendiff_long2$Generation <- as.factor(gendiff_long2$Generation)



gendiff_plot55 <- ggplot(gendiff_long2, aes(x=Generation, y=minallele_freq, colour=Treatment))
gendiff_plot56 <- gendiff_plot55 + geom_jitter() + ggtitle('Generation Spread') 
#+ guides(colour=T)
print(gendiff_plot56)

gendiff_plot44 <- ggplot(gendiff_long2, aes(x=Treatment, y=minallele_freq, colour=Generation))
gendiff_plot45 <- gendiff_plot44 + geom_jitter() + ggtitle('Generation Spread')
#+ guides(colour=T)
print(gendiff_plot45)

#gendiff_long2$Generation <- as.numeric(gendiff_long2$Generation)
#gendiff_plot66 <- ggplot(gendiff_long2, aes(x=Generation, y=minallele_freq, colour=pos))
#gendiff_plot65 <- gendiff_plot66 + geom_point() + ggtitle('Generation') + geom_
#+ guides(colour=T)
#print(gendiff_plot65)
dat




