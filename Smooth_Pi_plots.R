#Pi Plots: Smooth Function with Overlay plots of different generations:

#Function to read all and overlay plot
setwd("~/Bioinformatics/episodic_practice/Bowtie_Pi")
MyPi <- list.files(pattern="bowtie.pi")

for (file in MyPi){
  require(ggplot2)
  x2 <- gsub("\\_.*","",file)
  
  #Read in the data:
  Datt2 <- read.table(file)
  colnames(Datt2) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  Datt2$chr <- as.character(Datt2$chr)
  
  #Remove "na" pi values
  Datt2 <- Datt2[-which(Datt2$Pi=="na"),]
  
  #Need the numbers for chromosomes for labelling and colours:
  DattX <- Datt2[which(Datt2$chr=="X"),]
  a <- dim(DattX)[1]
  DattX$number <- 1:a
  
  Datt2L <- Datt2[which(Datt2$chr=="2L"),]
  b <- dim(Datt2L)[1]
  Datt2L$number <- (a+1):(a+b)
  
  Datt2R <- Datt2[which(Datt2$chr=="2R"),]
  c <- dim(Datt2R)[1]
  Datt2R$number <- (a+b+1):(a+b+c)
  
  Datt3L <- Datt2[which(Datt2$chr=="3L"),]
  d <- dim(Datt3L)[1]
  Datt3L$number <- (a+b+c+1):(a+b+c+d)
  
  Datt3R <- Datt2[which(Datt2$chr=="3R"),]
  e <- dim(Datt3R)[1]
  Datt3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  Datt4 <- Datt2[which(Datt2$chr=="4"),]
  f <- dim(Datt4)[1]
  Datt4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  #Full data frame of necessary chromosomes
  DattFull <- rbind(DattX, Datt2L, Datt2R, Datt3L, Datt3R, Datt4)
  
  #Pi as numeric
  DattFull$Pi=as.numeric(levels(DattFull$Pi))[DattFull$Pi]
  DattFull$Seq <- x2
  assign(paste("Datt", x2, sep="_"),DattFull)
}

ggplot(Datt_MGD3, aes(x = number, y= Pi, colour = chr)) + 
  geom_smooth(method = "loess", show.legend = F) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))  +  
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
  geom_smooth(method = "loess", data = Datt_F115ConR1, linetype = "dashed", aes(x = number, y= Pi, colour = chr), show_guide = F) +
  geom_smooth(method = "loess", data = Datt_F115ConR2, linetype = "dotted", aes(x = number, y= Pi, colour = chr), show_guide = F) +
  geom_smooth(method = "loess", data = Datt_F115SelR1, linetype = "longdash",  aes(x = number, y= Pi, colour = chr), show_guide = F) +
  geom_smooth(method = "loess", data = Datt_F115SelR2, linetype = "dotdash",  aes(x = number, y= Pi, colour = chr), show_guide = F)


#keep scale_x_discrete: Not colour (for SEQ?)


#DATTTT <- rbind(Datt_F115ConR1, Datt_F115ConR2, Datt_F115SelR1, Datt_F115SelR2, Datt_F38ConR1, Datt_F38ConR2, Datt_F38SelR1, Datt_F38SelR2, Datt_F77ConR1, Datt_F77ConR2, Datt_F77SelR1, Datt_F77SelR2, Datt_MGD3)  

#DATTTT <- rbind(Datt_F115ConR1, Datt_MGD3)  


#ggplot(DATTTT, aes(x = number, y= Pi, colour=Seq)) + 
#  geom_smooth(method = "loess", show.legend = F) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
#  xlab("") +
#  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
#  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15))

ggplot(Datt_MGD3, aes(x = number, y= Pi, linetype = chr)) + 
  geom_smooth(method = "loess", size=1.25) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#E69F00", 'wheat3', 'lemonchiffon4', "#56B4E9"))  +  
  # theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), legend.position="none") +
  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
  geom_smooth(method = "loess", data = Datt_F115ConR1, aes(x = number, y= Pi)) +
  geom_smooth(method = "loess", data = Datt_F115ConR2, aes(x = number, y= Pi)) +
  geom_smooth(method = "loess", data = Datt_F115SelR1,  aes(x = number, y= Pi)) +
  geom_smooth(method = "loess", data = Datt_F115SelR2,  aes(x = number, y= Pi))




DATTTT_115 <- rbind(Datt_F115ConR1, Datt_F115ConR2, Datt_F115SelR1, Datt_F115SelR2, Datt_MGD3)
DATTTT_77<- rbind(Datt_F77ConR1, Datt_F77ConR2, Datt_F77SelR1, Datt_F77SelR2, Datt_MGD3)
DATTTT_38 <- rbind(Datt_F38ConR1, Datt_F38ConR2, Datt_F38SelR1, Datt_F38SelR2, Datt_MGD3)

ggplot(DATTTT_38, aes(x = number, y= Pi, linetype = chr, colour=Seq)) +
  geom_smooth(method = "loess", size=1.25) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  scale_linetype_manual(values=c("solid", "solid", "dotted", "dotted", "solid", "dotted"), guide = "none") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#E69F00", 'darkorange1', 'cornflowerblue', "#56B4E9", 'grey0'))  + 
  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), legend.text=element_text(size=7.5))

ggplot(DATTTT_77, aes(x = number, y= Pi, linetype = chr, colour=Seq)) +
  scale_linetype_manual(values=c("solid", "solid", "dotted", "dotted", "solid", "dotted"), guide = "none") +
  geom_smooth(method = "loess", size=1.25) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#E69F00", 'darkorange1', 'cornflowerblue', "#56B4E9", 'grey0'))  + 
  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), legend.text=element_text(size=7.5))

ggplot(DATTTT_115, aes(x = number, y= Pi, linetype = chr, colour=Seq)) +
  scale_linetype_manual(values=c("solid", "solid", "dotted", "dotted", "solid", "dotted"), guide = "none") +
  geom_smooth(method = "loess", size=1.25) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#E69F00", 'darkorange1', 'cornflowerblue', "#56B4E9", 'grey0'))  + 
  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), legend.text=element_text(size=7.5))



