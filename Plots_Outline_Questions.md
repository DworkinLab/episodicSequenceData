# File with the outline of plots and questions about these first draft of plots

## Pi: Ancestral Pi for Novoalign:

### Outline

The ancestral nucleotide diversity:

![Ancestral Pi Plot for Novoalign](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Ancestral_Pi.png)

### Questions

  1. Necessary for all populations?
  
    -- Have all populations (and for bowtie and novoalign mappers)

  2. Average Pi for all mappers??
  
    -- Calculate bwa Pi and average between three? or show one (/2) mappers as a represenation?

  3. Overlay for changes in diversity over time?
   
    -- Do we want overlay plots with ~splines showing the change in diversity from Ancestor --> 115?


## Fst Plots:

### Outline

Average pairwise Fst between control and selection replicates

  -- Average b/w mappers and replicates

**Generation 38:**
![meanFst for F38](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F38_meanFstPlot.png)

**Generation 77:**
![meanFst for F77](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F77_meanFstPlot.png)

**Generation 115:** 
![meanFst for F115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_meanFstPlot.png)


### Questions

  **1. Can these values be adjusted? p.adjust for FST values??**
  
    -- Does it make sense to use a FDR adjustment on these values?
  
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
    
    -- Should I keep the top 50%? the top 10% Fst values?

## Model Outputs

### Outline

Plots for original values and FDR adjusted

**None corrected P values: TxG -log10(meanP-value)**
![FullGenomeTxGPlot](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/CHROMO_meanP.png)

**FDR Corrected P-values: TxG -log10(meanP-value)**
![FDRcorrection](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Fdr_adjustP.png)

### Questions: 

   1. is Bonferroni a better visualization for the paper (much less going on)
   
**TxG: -log10(meanP) with Bonferroni Correction for multiple comparisons**
![BonferroniCorrection_2200x1100](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Bonferroni_p.adjust_TxG.png)

 Advantage with this: Can create a plot on with valued of the regular plot (first one) with coloured sig. values: Would not look good with FDR:
  
  ![Coloured Sig](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/fdr_bonf_adjustP_sigColoured.png)
  


## Poolseq outputs:

### Outline

Output from PoolSeq package: the significant selection coeffients that were significant for Predation lines and not for controls

Ongoing with the slow pace of Poolseq: 3L and 3R almost completed

This is the average b/w two mappers (bwa and novoalign), keeping the least significant pvalue.

![poolseq_2L](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Chromo_2_selcoef.png)

### Questions

1. Plot like above (and all chromo eventually)? 

2. Any cut off for selection coefficients or just any significant selection coefficients unique to predator lines?


## Trajectories and positions:
 
### Outline
 
 Filtered positions for:
 
 -- pvalues <0.05 after FDR from model output 
 
 -- similar positions in the poolseq significant selection coefficients
 
 -- Then found any overlapping windows with these positions with Fst values != 0. 
 
  Ended up with ~400 positions for both 2L and 2R each
 
 Trajectories are the mean absolute difference the treatments had from the ancestor
 
 ![Trajectory_2ndChromo](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Trajectory_diff_2Chromo.png)

### Questions:

1. Plots of individual positions?

2. The location of these positions on one of the above plots 

  -- larger and coloured positions on the model output for example: 2L with FDR:
  
  ![]()
  
  

