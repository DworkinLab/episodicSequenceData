library(GWASTools)
setwd("/Users/paulknoops/Sequence_analysis_2016")
episodic_target40_13_24 <- read.table("episodic_data_target40_1-3,2-4.cmh.gwas", h=T)
manhattanPlot(episodic_target40_13_24$P, episodic_target40_13_24$CHR, main = "Episodic_target40", ylim = c(0, (log10(length(episodic_target40_13_24$P)) +35)))
summary(episodic_target40_13_24)
episodic_target40_13_24
?manhattanPlot
#snpCorrelationPlot(episodic_target40_13_24$P, episodic_target40_13_24$CHR, ylim = c(0, 5.0))
