#Episodic data analysis:

source("episodic_packages.R")


#episodic_long <- read.csv("X", h=TRUE)

  

### BWA
  
#episodic_data_3R.csv
#episodic_data_2R.csv
#episodic_data_3L.csv
#episodic_data_2L.csv
#episodic_data_4.csv
#episodic_data_X.csv
  
#### Bowtie
#episodic_data_bowtie_3R.csv
#episodic_data_bowtie_2R.csv
#episodic_data_bowtie_3L.csv
#episodic_data_bowtie_2L.csv
#episodic_data_bowtie_4.csv
#episodic_data_bowtie_X.csv

#The Data: in long format, each position with Treatment, Cage and Generation, along with the Major and Mnor allele counts correponding to the ancestral major/minor allele

#head(episodic_long)

#The full model:

#Call each position
position <- unique(episodic_long$pos)
no.pos <- length(position)

#Remove N/A's -- possibly not needed so hashed out.
#episodic_long <- na.omit(episodic_long)

#Make list to store model
modlist_2 <- as.list(1:no.pos)
#Each model of a position, named for each mod will be position
names(modlist_2) <- position

#Empty Data Frame to store all coeffecients of model
coeffs_df <- data.frame(NULL)


#Run the  model for each position
for(i in position){
  #Don't need to print this out necessarily, just nice to track
  print(paste("Running entity:", i, "which is", which(position==i), "out of", no.pos))
  
  #Temporary data frame for only the one position
  tmp2 <- episodic_long[episodic_long$pos == i,]
  
  #The model: major vs. minor counts by Treatment, Generation and Treatment:Generation
  modlist_2[[i]] <- 
    glm(cbind(Major_count, Minor_count) ~ Treatment*Generation,
        data = tmp2, family = "binomial")
  
  #Turn this model into a data frame of coefficients
  x <- as.data.frame(summary(modlist_2[[i]] )$coefficients, row.names = FALSE)
  
  #Name the columns for the model; i.e the levels of Effects
  x$Effects <- c("Intercept", "TreatmentSel", "Generation", "TreatmentsSel:Generation")
  
  #Name the position of this model results with i
  x$position <- i
  
  #Add to data frame (total for the whole data set == coeffs_df + the newly made X)
  coeffs_df <- rbind(coeffs_df, x)
  
  #Remove i for safety and it starts over
  rm(i)
}


#head(coeffs_df)

#Make the p-values into -log10 p values
coeffs_df$log_p <- -log10(coeffs_df$`Pr(>|z|)`)

#Save data frame???


#gg1 <- ggplot(data = coeffs_df, aes(x=position, y=log_p, colour=Effects))
#gg1 + geom_point(size = 1) + geom_hline(yintercept = 1.5)
