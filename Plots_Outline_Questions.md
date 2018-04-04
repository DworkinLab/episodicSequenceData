# File with the outline of plots and questions about these first draft of plots

## Pi: Ancestral Pi for Novoalign:

### Outline

The ancestral nucleotide diversity:

![Ancestral Pi Plot for Novoalign](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Ancestral_Pi.png)

### Questions

  1. Necessary for all populations? Have all populations (and for bowtie and novoalign mappers)

  2. Average Pi for all mappers?? 

  3. Overlay for changes in diversity over time?


## Fst Plots:

### Outline

Average pairwise Fst between control and selection replicates

**Generation 38:**
![meanFst for F38](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F38_meanFstPlot.png)

**Generation 77:**
![meanFst for F77](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F77_meanFstPlot.png)

**Generation 115:** 
![meanFst for F115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_meanFstPlot.png)


### Questions

  1. Can these values be adjusted? p.adjust for FST values??
   
  2. Downscaling (available for all generations): Necessary and proposed methods
  
    -- previous ideas were (Fst_C:C + Fst_S:S)/2 for scaling
  ___________________________________________________
**meanFst: Selection vs. Control: Generation115**
![meanFst for F115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_meanFstPlot.png)
**meanFst: Control vs. Control: Generation115**
![Fst_Con:Con_115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_Control:Control_FST.png)
**meanFst: Selection vs. Selection: Generation115**
![Fst_Sel:Sel_115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_Selection:Selection_FST.png)
____________________________________________________
  
  3. Cut off for positions? 
  
    -- Currently keeping anything with an Fst value for Con:Sel_115 comparison: any way to filter more deeply for peeks 


