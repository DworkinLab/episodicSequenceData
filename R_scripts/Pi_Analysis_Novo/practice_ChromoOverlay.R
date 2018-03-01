
#If you want to make a plot like this (at a generation)
DATTTT_115_2R <- DATTTT_115[which(DATTTT_115$chr=='2R'),]
head(DATTTT_115_2R)

gg_115_2R <- ggplot(DATTTT_115_2R, aes(x = window, y= Pi, linetype = Seq)) +
  scale_linetype_manual(values=c("solid", "dotted", "dashed", "longdash", "twodash")) +
  geom_smooth(method = "loess", size=1.25, color="#56B4E9") +   
  scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), legend.text=element_text(size=7.5))

print(gg_115_2R)






#### FORGET ALL THIS: EASIER METHOD (BELOW IS A POSSIBLE METHOD IF THIS DOES NOT WORK I AM KEEPING AS A BACKUP)

## Either do the same start (create many files with all chromosome) and subset for each chromosome: Will need to to the multiple times for each chromosome:
Datt_F115ConR1_3L <- Datt_F115ConR1[which(Datt_F115ConR1$chr=="3L"),]

####

# OR....... easier option: change the loop:

####
setwd("~/Bioinformatics/R-projects_git/episodicSequenceData/R_scripts/Pi_Analysis_Novo")
require('tidyverse')
MyPi <- list.files(pattern="novo.pi")

for (file in MyPi){
  require(ggplot2)
  x2 <- gsub("\\_.*","",file)
  
  #Read in the data:
  Datt2 <- read.table(file)
  colnames(Datt2) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  Datt2$chr <- as.character(Datt2$chr)
  #Remove "na" pi values
  Datt2 <- Datt2[-which(Datt2$Pi=="na"),]
  
  #Change Hash marks as desired:
  #DattX <- Datt2[which(Datt2$chr=="X"),]
  #DattFull <- DattX
  
  Datt2L <- Datt2[which(Datt2$chr=="2L"),]
  DattFull <- Datt2L
  
  #Datt2R <- Datt2[which(Datt2$chr=="2R"),]
  #DattFull <- Datt2R
  
  #Datt3L <- Datt2[which(Datt2$chr=="3L"),]
  #DattFull <- Datt3L
  
  #Datt3R <- Datt2[which(Datt2$chr=="3R"),]
  #DattFull <- Datt3R
  
  #Datt4 <- Datt2[which(Datt2$chr=="4"),]
  #DattFull <- Datt4

    #Pi as numeric
  DattFull$Pi=as.numeric(levels(DattFull$Pi))[DattFull$Pi]
  DattFull$Seq <- x2
  
  #name chrmosome to be more careful:
  #assign(paste("DattX", x2, sep="_"),DattFull)
  assign(paste("Datt2L", x2, sep="_"),DattFull)
  #assign(paste("Datt2R", x2, sep="_"),DattFull)
  #assign(paste("Datt3L", x2, sep="_"),DattFull)
  #assign(paste("Datt3R", x2, sep="_"),DattFull)
  #assign(paste("Datt4", x2, sep="_"),DattFull)
}

#should have many Datt2L_etc.. files: with only chromosome X

#create ggplots similarly: Not sure which you use, but the basic things to change should be the same:

  #change x = to window

  # REMOVE: #scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +

  # change scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))+
      #Just one input colour (colour that matches that chromosome)

  #(and obviously data input/output names)

head(Datt2L_MGD3)
datt2LGGPLOT <- ggplot(Datt2L_MGD3, aes(x = window, y= Pi, colour = chr)) + 
  geom_smooth(method = "loess", show.legend = F) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  #scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  #scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))  +  
  scale_colour_manual(values=c("#56B4E9"))  +
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
  geom_smooth(method = "loess", data = Datt2L_F115ConR1, linetype = "dashed", aes(x = window, y= Pi, colour = chr), show_guide = F) +
  geom_smooth(method = "loess", data = Datt2L_F115ConR2, linetype = "dotted", aes(x = window, y= Pi, colour = chr), show_guide = F) +
  geom_smooth(method = "loess", data = Datt2L_F115SelR1, linetype = "longdash",  aes(x = window, y= Pi, colour = chr), show_guide = F) +
  geom_smooth(method = "loess", data = Datt2L_F115SelR2, linetype = "dotdash",  aes(x = window, y= Pi, colour = chr), show_guide = F)

print(datt2LGGPLOT)






##
setwd("~/Bioinformatics/R-projects_git/episodicSequenceData/R_scripts/Pi_Analysis_Novo")
require('tidyverse')

#Can change the files you want to look at (change pattern of directories and move things)
MyPi <- list.files(pattern="novo.pi")

#change chromosomes of interest (can just be one if you like)
#mychromo <- c('X', '2L', '2R','3L','3R','4')
mychromo <- c('2L')

for (file in MyPi){
  require(ggplot2)
  x2 <- gsub("\\_.*","",file)
  Datt2 <- read.table(file)
  colnames(Datt2) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  Datt2$chr <- as.character(Datt2$chr)
  Datt2 <- Datt2[-which(Datt2$Pi=="na"),]
  for (chromo in mychromo){
    DattFull <- Datt2[which(Datt2$chr==chromo),]
  DattFull$Pi=as.numeric(levels(DattFull$Pi))[DattFull$Pi]
  DattFull$Seq <- x2

  assign(paste("Datt", chromo, x2, sep="_"),DattFull)
  }
}

Dddat_115_2L <- rbind(Datt_2L_F115ConR1, Datt_2L_F115ConR2, Datt_2L_F115SelR1, Datt_2L_F115SelR2)
Dddat_77_2L <-rbind(Datt_2L_F77ConR1, Datt_2L_F77ConR2, Datt_2L_F77SelR1, Datt_2L_F77SelR2)
Dddat_38_2L <- rbind(Datt_2L_F38ConR1, Datt_2L_F38ConR2, Datt_2L_F38SelR1, Datt_2L_F38SelR2)

pplot <- ggplot(Dddat_115_2L, aes(x = window, y= Pi, linetype = Seq, colour= Seq)) +
  #geom_point() +
  scale_linetype_manual(values=c("solid", "dotted", "dashed", "longdash", "twodash")) +
  geom_smooth(method = "loess", size=1.25, color="#56B4E9") +   
  scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), legend.text=element_text(size=7.5)) 
pplot


Datt_SelR1 <- rbind(Datt_2L_F38SelR1, Datt_2L_F77SelR1, Datt_2L_F115SelR1)

pplot <- ggplot(Datt_SelR1, aes(x = window, y= Pi, colour= Seq)) +
  #geom_point() +
  scale_color_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'wheat3' )) +
  geom_smooth(method = "loess", size=1) +   
  scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), legend.text=element_text(size=7.5)) 
pplot




