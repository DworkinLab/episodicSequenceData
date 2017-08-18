#Episodic data analysis:

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

episodic_long <- read.csv("episodic_data_bowtie_4.csv", h=TRUE)
#The .csv files will be in an odd format: may need to loop like with Sync_to_counts_smallerSyncs.R

#The Data: in long format, each position with Treatment, Cage and Generation, along with the Major and Mnor allele counts correponding to the ancestral major/minor allele

#The full model:

#Call each position
position <- unique(episodic_long$pos)
no.pos <- length(position)

#Remove N/A's -- possibly not needed so hashed out.
episodic_long <- na.omit(episodic_long)

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
  x <- as.data.frame(summary(modlist_2[[i]] )$coefficients)
  
  #Name the position of this model results with i
  x$position <- i
  
  #Add to data frame (total for the whole data set == coeffs_df + the newly made X)
  coeffs_df <- rbind(coeffs_df, x)
  
  #Remove i for safety and it starts over
  rm(i)
}

#Change column names to workable
colnames(coeffs_df) <- c("Estimate", "Standard_error", "z-value", "p-value", "position")

coeffs_df$Effects<-rownames(coeffs_df)


coeffs_df$Effects_2 <- ifelse(grepl("TreatmentSel:Generation",coeffs_df$Effects),'T_Sel:Gen', ifelse(grepl("Intercept",coeffs_df$Effects),'Int', coeffs_df$Effects ))

coeffs_df$Effects_2 <- ifelse(grepl("TreatmentSel",coeffs_df$Effects_2),'T_Sel', ifelse(grepl("Generation",coeffs_df$Effects_2),'Gen', coeffs_df$Effects_2))


rownames(coeffs_df) <- c()

#Make the p-values into -log10 p values
coeffs_df$log_p <- -log10(coeffs_df$`p-value`)

coeffs_df <- subset(coeffs_df, select = -c(Effects))

coeffs_df$Effects <- ifelse(coeffs_df$Effects_2=='T_Sel', 'TreatmentSel', ifelse(coeffs_df$Effects_2=='Gen', 'Generation', ifelse(coeffs_df$Effects_2=='Int', 'Intercept', 'TreatmentSel:Generation')))

coeffs_df <- subset(coeffs_df, select = -c(Effects_2))


coeffs_df <- coeffs_df[-which(coeffs_df$log_p==0),]

#coeffs_df <- coeffs_df[grep("Treatment", coeffs_df$Effects), ]
#coeffs_df <- coeffs_df[which(coeffs_df$Effects=="TreatmentSel"),]


#Save as Data frame
#Change based on the data input!

write.csv(coeffs_df, file="episodic_bowtie_4_coeffs.csv", row.names = FALSE)

