#MGD3 Plotting - Bowtie
bowAnc <- read.table("MGD3_bowtie.pi")
colnames(bowAnc) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')

bowAnc$chr <- as.character(bowAnc$chr)
bowAnc2 <- bowAnc[-which(bowAnc$chr=="YHet"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="2RHet"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="2LHet"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="3LHet"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="3RHet"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="U"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="XHet"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="dmel_mitochondrion_genome"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$chr=="Uextra"),]
bowAnc2 <- bowAnc2[-which(bowAnc2$Pi=="na"),]

BowAnc3R <- bowAnc2[which(bowAnc2$chr=="3R"),]
BowAnc3R$number <- 8565:11253
BowAnc3L <- bowAnc2[which(bowAnc2$chr=="3L"),]
BowAnc3L$number <- 6253:8564
BowAnc2R <- bowAnc2[which(bowAnc2$chr=="2R"),]
BowAnc2R$number <- 4259:6252
BowAnc2L <- bowAnc2[which(bowAnc2$chr=="2L"),]
BowAnc2L$number <- 2092:4258
BowAncX <- bowAnc2[which(bowAnc2$chr=="X"),]
BowAncX$number <- 1:2091
BowAnc4 <- bowAnc2[which(bowAnc2$chr=="4"),]
BowAnc4$number <- 11254:11372

BowAnc3 <- rbind(BowAncX, BowAnc2L, BowAnc2R, BowAnc3L, BowAnc3R, BowAnc4)
head(BowAnc3)

#Anc3$number <- 1:11419
require(ggplot2)

BowAnc3$Pi=as.numeric(levels(BowAnc3$Pi))[BowAnc3$Pi]

BowGG <- ggplot(BowAnc3, aes(x = number, y= Pi, colour = chr)) 

Bowgg2 <- BowGG + 
  geom_point(size=0.3, show.legend = F) +
  scale_y_continuous(limits=c(0, 0.02), breaks=seq(0, 0.02, 0.005)) + 
  xlab("") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))

Bowgg2
