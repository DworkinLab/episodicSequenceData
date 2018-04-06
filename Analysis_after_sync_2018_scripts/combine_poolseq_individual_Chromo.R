#setwd('/home/paul/episodicData/novoalign/novo_mpileup/splitsync_dir/novo_episodic_2L_Con_Split')
setwd('/home/paul/episodicData/novoalign/novo_mpileup/splitsync_dir/novo_episodic_2L_Sel_Split')

mycsvs <- list.files(pattern='.csv')
X <- NULL
for (file in mycsvs){
  X2 <- read.csv(file, h=T)
  X <- rbind(X, X2)
}
X$chr <- '2L'
#write.csv(X, file='/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_2L_Con.csv', row.names = FALSE)
write.csv(X, file='/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_2L_Sel.csv', row.names = FALSE)
