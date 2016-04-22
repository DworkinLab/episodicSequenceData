Step 1) raw_data file made and all raw reads into file
Step 2) md5sum for each different md5.txt
```
[paul@info114 raw_data]$ md5sum -c md5_F115_F38.txt 
F115ConR1_TAGCTT_L001_R1_001.fastq.gz: OK
F115ConR1_TAGCTT_L001_R2_001.fastq.gz: OK
F115ConR1_TAGCTT_L002_R1_001.fastq.gz: OK
F115ConR1_TAGCTT_L002_R2_001.fastq.gz: OK
F115ConR2_GGCTAC_L001_R1_001.fastq.gz: OK
F115ConR2_GGCTAC_L001_R2_001.fastq.gz: OK
F115ConR2_GGCTAC_L002_R1_001.fastq.gz: OK
F115ConR2_GGCTAC_L002_R2_001.fastq.gz: OK
F115SelR1_GTTTCG_L001_R1_001.fastq.gz: OK
F115SelR1_GTTTCG_L001_R2_001.fastq.gz: OK
F115SelR1_GTTTCG_L002_R1_001.fastq.gz: OK
F115SelR1_GTTTCG_L002_R2_001.fastq.gz: OK
F115SelR2_GTGGCC_L001_R1_001.fastq.gz: OK
F115SelR2_GTGGCC_L001_R2_001.fastq.gz: OK
F115SelR2_GTGGCC_L002_R1_001.fastq.gz: OK
F115SelR2_GTGGCC_L002_R2_001.fastq.gz: OK
F38ConR1_ATCACG_L001_R1_001.fastq.gz: OK
F38ConR1_ATCACG_L001_R2_001.fastq.gz: OK
F38ConR1_ATCACG_L002_R1_001.fastq.gz: OK
F38ConR1_ATCACG_L002_R2_001.fastq.gz: OK
F38ConR2_TTAGGC_L001_R1_001.fastq.gz: OK
F38ConR2_TTAGGC_L001_R2_001.fastq.gz: OK
F38ConR2_TTAGGC_L002_R1_001.fastq.gz: OK
F38ConR2_TTAGGC_L002_R2_001.fastq.gz: OK
F38SelR1_ACTTGA_L001_R1_001.fastq.gz: OK
F38SelR1_ACTTGA_L001_R2_001.fastq.gz: OK
F38SelR1_ACTTGA_L002_R1_001.fastq.gz: OK
F38SelR1_ACTTGA_L002_R2_001.fastq.gz: OK
F38SelR2_GATCAG_L001_R1_001.fastq.gz: OK
F38SelR2_GATCAG_L001_R2_001.fastq.gz: OK
F38SelR2_GATCAG_L002_R1_001.fastq.gz: OK
F38SelR2_GATCAG_L002_R2_001.fastq.gz: OK
```

```
[paul@info114 raw_data]$ md5sum -c md5_F77.txt 
Con_R1_F77_ATGTCA_L003_R1_001.fastq.gz: FAILED
Con_R1_F77_ATGTCA_L003_R2_001.fastq.gz: FAILED
Con_R1_F77_ATGTCA_L004_R1_001.fastq.gz: FAILED
Con_R1_F77_ATGTCA_L004_R2_001.fastq.gz: FAILED
Con_R2_F77_ATTCCT_L003_R1_001.fastq.gz: FAILED
Con_R2_F77_ATTCCT_L003_R2_001.fastq.gz: FAILED
Con_R2_F77_ATTCCT_L004_R1_001.fastq.gz: FAILED
Con_R2_F77_ATTCCT_L004_R2_001.fastq.gz: FAILED
Sel_R1_F77_TTAGGC_L003_R1_001.fastq.gz: FAILED
Sel_R1_F77_TTAGGC_L003_R2_001.fastq.gz: FAILED
Sel_R1_F77_TTAGGC_L004_R1_001.fastq.gz: FAILED
Sel_R1_F77_TTAGGC_L004_R2_001.fastq.gz: FAILED
Sel_R2_F77_GATCAG_L003_R1_001.fastq.gz: OK
Sel_R2_F77_GATCAG_L003_R2_001.fastq.gz: OK
Sel_R2_F77_GATCAG_L004_R1_001.fastq.gz: OK
Sel_R2_F77_GATCAG_L004_R2_001.fastq.gz: OK
md5sum: WARNING: 12 of 16 computed checksums did NOT match
```

```
MGD also failed (long with extra files in md5.txt)
```

Test in scratch folder (* on the head)

Step 3)



