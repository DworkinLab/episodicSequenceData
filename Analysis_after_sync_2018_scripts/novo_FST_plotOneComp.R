#Create a Plot of a single .csv file (created from https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Fst_Split_Comparisons.R)

# Plotting: One .csv file

require(data.table)
require(tidyverse)
setwd("/Users/paulknoops/Bioinformatics/Analysis_EpisodicSequenceData_2018/FST")

ddat2 <- fread('fst_1:2.csv')
tttle <- 'PLOT TITLE'

ddat2$num <- 1:length(ddat2$meanFst)

xc_mean <- mean(ddat2$meanFst)
xs_median <- median(ddat2$meanFst)

# Keep chrososomes of interest (in list below) if needed: 

#Obviously Unhash if needed:
# ddat2 <- filter(ddat2, chr %in% c('X','2L','2R','3L','3R','4')

g <- nrow(ddat2[which(ddat2$chr=='2L'),])
h <- nrow(ddat2[which(ddat2$chr=='2R'),])
i <- nrow(ddat2[which(ddat2$chr=='3L'),])
j <- nrow(ddat2[which(ddat2$chr=='3R'),])
k <- nrow(ddat2[which(ddat2$chr=='4'),])
l <- nrow(ddat2[which(ddat2$chr=='X'),])

#To change the order for X to be first:
# need to figure out the order the chromosomes are put in the data frame to give each the correct number in the sequence

#Example: I want X first but the rest the same, this will have X have numbers 1:l (last line) and then start with 2L (first line)

#2L-2R-3L-3R-4-X
ddat2$number <- c((l+1):(l+g), 
                  (l+g+1):(l+g+h), 
                  (l+g+h+1):(l+g+h+i),
                  (l+g+h+i+1):(l+g+h+i+j),
                  (l+g+h+i+j+1):(l+g+h+i+j+k), 
                  (1:l))

# Example: different order but same effect, X is second but want it first, so second line becomes 1:l. 2L follows ((l+1):(l+g)) and so on.
#3R-X-2L-3L-2R-4
#ddat2$number <- c((l+g+h+i+1):(l+g+h+i+j),
#                  (1:l),
#                  (l+1):(l+g), 
#                  (l+g+h+1):(l+g+h+i),
#                  (l+g+1):(l+g+h),
#                  (l+g+h+i+j+1):(l+g+h+i+j+k))

#Basically: each chromosome will correspond to one split, just need to move it around based initial order (assumeing the order you want is X-2L-2R-3L-3R-4)

ggxcv <-  ggplot(data = ddat2, aes(x=num, y=meanFst, color=chr))
ggxcv2 <- ggxcv + 
  geom_point(size=0.5, show.legend = F, alpha = 0.75) + 
  theme(panel.background = element_blank()) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
  geom_hline(yintercept = xc_mean) + 
  xlab("Chromosome") +
  ggtitle(tttle) +
  scale_x_discrete(limits=c(l/2, l+(g/2), (l+g+(h/2)), (l+g+h+(i/2)), (l+g+h+i+(j/2)), (l+g+h+i+j+(k/2))), labels = c("X","2L", "2R", '3L', '3R', "4")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))

print(ggxcv2)
