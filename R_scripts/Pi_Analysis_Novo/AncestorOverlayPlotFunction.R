#AncestorOverlayPlotFunction <- function(x, x2, y, y2, z, z2) {
  require(ggplot2)
  setwd("~/Bioinformatics/episodic_practice/MGD3")
  x <- "MGD3_SO_CAGATC_novo.pi"
  y <- "MGD3_bowtie.pi"
  z <- "MGD3.pi"
  x2 <- "novo"
  y2 <- "bowtie"
  z2 <- "bwa"
  #Read in the data:
  
  Datt_X <- read.table(x)
  colnames(Datt_X) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  #Remove unnecessary regions: Not necessary based on later steps
  Datt_X$chr <- as.character(Datt_X$chr)
  
  #Remove "na" pi values
  Datt_X <- Datt_X[-which(Datt_X$Pi=="na"),]
  
  #Need the numbers for chromosomes for labelling and colours:
  DattX <- Datt_X[which(Datt_X$chr=="X"),]
  a <- dim(DattX)[1]
  DattX$number <- 1:a
  
  Datt2L <- Datt_X[which(Datt_X$chr=="2L"),]
  b <- dim(Datt2L)[1]
  Datt2L$number <- (a+1):(a+b)
  
  Datt2R <- Datt_X[which(Datt_X$chr=="2R"),]
  c <- dim(Datt2R)[1]
  Datt2R$number <- (a+b+1):(a+b+c)
  
  Datt3L <- Datt_X[which(Datt_X$chr=="3L"),]
  d <- dim(Datt3L)[1]
  Datt3L$number <- (a+b+c+1):(a+b+c+d)
  
  Datt3R <- Datt_X[which(Datt_X$chr=="3R"),]
  e <- dim(Datt3R)[1]
  Datt3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  Datt4 <- Datt_X[which(Datt_X$chr=="4"),]
  f <- dim(Datt4)[1]
  Datt4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  #Full data frame of necessary chromosomes
  DattFull_X <- rbind(DattX, Datt2L, Datt2R, Datt3L, Datt3R, Datt4)
  
  
  
  #Pi as numeric
  DattFull_X$Pi=as.numeric(levels(DattFull_X$Pi))[DattFull_X$Pi]
  
  DattFull_X$Seq <- x2
  assign(paste("Datt", x2, sep="_"),DattFull_X)

  
  ####
  
  Datt_Y <- read.table(y)
  colnames(Datt_Y) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  #Remove unnecessary regions: Not necessary based on later steps
  Datt_Y$chr <- as.character(Datt_Y$chr)
  
  #Remove "na" pi values
  Datt_Y <- Datt_Y[-which(Datt_Y$Pi=="na"),]
  
  #Need the numbers for chromosomes for labelling and colours:
  DattX <- Datt_Y[which(Datt_Y$chr=="X"),]
  a <- dim(DattX)[1]
  DattX$number <- 1:a
  
  Datt2L <- Datt_Y[which(Datt_Y$chr=="2L"),]
  b <- dim(Datt2L)[1]
  Datt2L$number <- (a+1):(a+b)
  
  Datt2R <- Datt_Y[which(Datt_Y$chr=="2R"),]
  c <- dim(Datt2R)[1]
  Datt2R$number <- (a+b+1):(a+b+c)
  
  Datt3L <- Datt_Y[which(Datt_Y$chr=="3L"),]
  d <- dim(Datt3L)[1]
  Datt3L$number <- (a+b+c+1):(a+b+c+d)
  
  Datt3R <- Datt_Y[which(Datt_Y$chr=="3R"),]
  e <- dim(Datt3R)[1]
  Datt3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  Datt4 <- Datt_Y[which(Datt_Y$chr=="4"),]
  f <- dim(Datt4)[1]
  Datt4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  #Full data frame of necessary chromosomes
  DattFull_Y <- rbind(DattX, Datt2L, Datt2R, Datt3L, Datt3R, Datt4)
  
  
  
  #Pi as numeric
  DattFull_Y$Pi=as.numeric(levels(DattFull_Y$Pi))[DattFull_Y$Pi]
  
  DattFull_Y$Seq <- y2
  assign(paste("Datt", y2, sep="_"),DattFull_Y)
 
   ####
  Datt_Z <- read.table(z)
  colnames(Datt_Z) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  #Remove unnecessary regions: Not necessary based on later steps
  Datt_Z$chr <- as.character(Datt_Z$chr)
  
  #Remove "na" pi values
  Datt_Z <- Datt_Z[-which(Datt_Z$Pi=="na"),]
  
  #Need the numbers for chromosomes for labelling and colours:
  DattX <- Datt_Z[which(Datt_Z$chr=="X"),]
  a <- dim(DattX)[1]
  DattX$number <- 1:a
  
  Datt2L <- Datt_Z[which(Datt_Z$chr=="2L"),]
  b <- dim(Datt2L)[1]
  Datt2L$number <- (a+1):(a+b)
  
  Datt2R <- Datt_Z[which(Datt_Z$chr=="2R"),]
  c <- dim(Datt2R)[1]
  Datt2R$number <- (a+b+1):(a+b+c)
  
  Datt3L <- Datt_Z[which(Datt_Z$chr=="3L"),]
  d <- dim(Datt3L)[1]
  Datt3L$number <- (a+b+c+1):(a+b+c+d)
  
  Datt3R <- Datt_Z[which(Datt_Z$chr=="3R"),]
  e <- dim(Datt3R)[1]
  Datt3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  Datt4 <- Datt_Z[which(Datt_Z$chr=="4"),]
  f <- dim(Datt4)[1]
  Datt4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  #Full data frame of necessary chromosomes
  DattFull_Z <- rbind(DattX, Datt2L, Datt2R, Datt3L, Datt3R, Datt4)
  
  
  
  #Pi as numeric
  DattFull_Z$Pi=as.numeric(levels(DattFull_Z$Pi))[DattFull_Z$Pi]
  
  DattFull_Z$Seq <- z2
  assign(paste("Datt", z2, sep="_"),DattFull_Z)
  
  ##PLOTS:
  XGG <- ggplot(DattFull_X, aes(x = number, y= Pi, colour = chr)) 
  YGG <- ggplot(DattFull_Y, aes(x = number, y= Pi, colour = chr))
  ZGG <- ggplot(DattFull_Z, aes(x = number, y= Pi, colour = chr))
  
  X_plot <- XGG + geom_smooth(method = "loess", show.legend = F) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) +
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))
  
  print(X_plot)
  
  
  Y_plot <- YGG + geom_smooth(method = "loess", show.legend = F, linetype = "dotted") +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))
  print(Y_plot)
  
  Z_plot <- ZGG + geom_smooth(method = "loess", show.legend = F, linetype = "dashed") +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))
  print(Z_plot)
  
  ggplot(DattFull_X, aes(x = number, y= Pi, colour = chr)) + 
    geom_smooth(method = "loess", show.legend = T, linetype = "dotted") +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))  +  
    theme(text = element_text(size=20), 
          axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
    geom_smooth(method = "loess", linetype="dashed", data = DattFull_Y, aes(x = number, y= Pi, colour = chr), show_guide = TRUE) +
    geom_smooth(method = "loess", data = DattFull_Z, aes(x = number, y= Pi, colour = chr), show_guide = TRUE) 
    

head(Datt_bwa)
head(Datt_bowtie)
head(Datt_novo)
DATTT <- rbind(Datt_bowtie, Datt_bwa, Datt_novo)

ggplot(DATTT, aes(x = number, y= Pi, colour = chr)) + 
  geom_smooth(method = "loess", show.legend = T, aes(linetype=Seq)) +   scale_y_continuous(limits=c(0, 0.009), breaks=seq(0, 0.009, 0.001)) + 
  xlab("") +
  scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))  +  
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15))
