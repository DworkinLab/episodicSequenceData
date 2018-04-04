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
  
  Ex. FDR adjust Fst:
  ![FDR_Adjust_FST](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Fst_FDRAdjust_Sel:Con_115.png)
   
  2. Downscaling (available for all generations): Necessary? and methods?
  
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

## Model Outputs

### Outline

**None corrected P values: TxG -log10(meanP-value)**
![FullGenomeTxGPlot](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/CHROMO_meanP.png)

**FDR Corrected P-values: TxG -log10(meanP-value)**
![FDRcorrection](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Fdr_adjustP.png)

### Questions: 
   1. is Bonferroni a better visualization for the paper (much less going on)
   
**TxG: -log10(meanP) with Bonferroni Correction for multiple comparisons**
![BonferroniCorrection_2200x1100](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Bonferroni_p.adjust_TxG.png)

  -- Can create a plot one the regular plot (first one) with coloured sig. values: Would not look good with FDR:
  
  ![Coloured Sig](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fdr_bonf_adjustP_sigColoured.png)
  
  
  2. The size of plots
  
**Bonferroni: smaller size (comparison for plot sizes)**
![Bonferonni_1100x550](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Bonf_AdjustP.png)






