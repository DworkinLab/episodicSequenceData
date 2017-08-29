#MGD3 Plotting - BWA
Anc <- read.table("MGD3.pi")
colnames(Anc) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')

Anc$chr <- as.character(Anc$chr)
Anc2 <- Anc[-which(Anc$chr=="YHet"),]
Anc2 <- Anc2[-which(Anc2$chr=="2RHet"),]
Anc2 <- Anc2[-which(Anc2$chr=="2LHet"),]
Anc2 <- Anc2[-which(Anc2$chr=="3LHet"),]
Anc2 <- Anc2[-which(Anc2$chr=="3RHet"),]
Anc2 <- Anc2[-which(Anc2$chr=="U"),]
Anc2 <- Anc2[-which(Anc2$chr=="XHet"),]
Anc2 <- Anc2[-which(Anc2$chr=="dmel_mitochondrion_genome"),]
Anc2 <- Anc2[-which(Anc2$chr=="Uextra"),]
Anc2 <- Anc2[-which(Anc2$Pi=="na"),]

Anc3R <- Anc2[which(Anc2$chr=="3R"),]
Anc3R$number <- 8606:11297
Anc3L <- Anc2[which(Anc2$chr=="3L"),]
Anc3L$number <- 6281:8605
Anc2R <- Anc2[which(Anc2$chr=="2R"),]
Anc2R$number <- 4274:6280
Anc2L <- Anc2[which(Anc2$chr=="2L"),]
Anc2L$number <- 2099:4273
AncX <- Anc2[which(Anc2$chr=="X"),]
AncX$number <- 1:2098
Anc4 <- Anc2[which(Anc2$chr=="4"),]
Anc4$number <- 11298:11419

Anc3 <- rbind(AncX, Anc2L, Anc2R, Anc3L, Anc3R, Anc4)
head(Anc3)

#Anc3$number <- 1:11419
require(ggplot2)

Anc3$Pi=as.numeric(levels(Anc3$Pi))[Anc3$Pi]

gganc <- ggplot(Anc3, aes(x = number, y= Pi, colour = chr)) 

gX <- gganc + 
  geom_point(size=0.3, show.legend = F) +
  scale_y_continuous(limits=c(0, 0.02), breaks=seq(0, 0.02, 0.005)) + 
  xlab("") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))

gX
