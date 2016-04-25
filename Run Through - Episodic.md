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


Appears F77 files had some errors, will run through with them but re do later with a new transfer from Ian


Step 3) Fastqc quality control
```
mkdir /home/paul/episodicData/fastqc
```
```
fastqc -o /home/paul/episodicData/fastqc /home/paul/episodicData/raw_data/*.fastq.gz
```

Step 4) Renaming
For some of the steps, need a specific naming style, matching F38SelR2_GATCAG_L001_R2_001.fastq.gz
Change names using mv and copy and paste from below necessary files (do from raw_dir)
ex. from raw_dir
```
mv Con_R1_F77_ATGTCA_L003_R1_001.fastq.gz F77ConR1_ATGTCA_L001_R1_001.fastq.gz
```

```
Con_R1_F77_ATGTCA_L003_R1_001.fastq.gz F77ConR1_ATGTCA_L001_R1_001.fastq.gz
Con_R1_F77_ATGTCA_L003_R2_001.fastq.gz F77ConR1_ATGTCA_L001_R2_001.fastq.gz
Con_R1_F77_ATGTCA_L004_R1_001.fastq.gz F77ConR1_ATGTCA_L002_R1_001.fastq.gz
Con_R1_F77_ATGTCA_L004_R2_001.fastq.gz F77ConR1_ATGTCA_L002_R2_001.fastq.gz
Con_R2_F77_ATTCCT_L003_R1_001.fastq.gz F77ConR2_ATTCCT_L001_R1_001.fastq.gz
Con_R2_F77_ATTCCT_L003_R2_001.fastq.gz F77ConR2_ATTCCT_L001_R2_001.fastq.gz
Con_R2_F77_ATTCCT_L004_R1_001.fastq.gz F77ConR2_ATTCCT_L002_R1_001.fastq.gz
Con_R2_F77_ATTCCT_L004_R2_001.fastq.gz F77ConR2_ATTCCT_L002_R2_001.fastq.gz
F115ConR1_TAGCTT_L001_R1_001.fastq.gz OK
F115ConR1_TAGCTT_L001_R2_001.fastq.gz OK
F115ConR1_TAGCTT_L002_R1_001.fastq.gz OK
F115ConR1_TAGCTT_L002_R2_001.fastq.gz OK
F115ConR2_GGCTAC_L001_R1_001.fastq.gz OK
F115ConR2_GGCTAC_L001_R2_001.fastq.gz OK
F115ConR2_GGCTAC_L002_R1_001.fastq.gz OK
F115ConR2_GGCTAC_L002_R2_001.fastq.gz OK
F115SelR1_GTTTCG_L001_R1_001.fastq.gz OK
F115SelR1_GTTTCG_L001_R2_001.fastq.gz OK
F115SelR1_GTTTCG_L002_R1_001.fastq.gz OK
F115SelR1_GTTTCG_L002_R2_001.fastq.gz OK
F115SelR2_GTGGCC_L001_R1_001.fastq.gz OK
F115SelR2_GTGGCC_L001_R2_001.fastq.gz OK
F115SelR2_GTGGCC_L002_R1_001.fastq.gz OK
F115SelR2_GTGGCC_L002_R2_001.fastq.gz OK
F38ConR1_ATCACG_L001_R1_001.fastq.gz OK
F38ConR1_ATCACG_L001_R2_001.fastq.gz OK
F38ConR1_ATCACG_L002_R1_001.fastq.gz OK
F38ConR1_ATCACG_L002_R2_001.fastq.gz OK
F38ConR2_TTAGGC_L001_R1_001.fastq.gz OK
F38ConR2_TTAGGC_L001_R2_001.fastq.gz OK
F38ConR2_TTAGGC_L002_R1_001.fastq.gz OK
F38ConR2_TTAGGC_L002_R2_001.fastq.gz OK
F38SelR1_ACTTGA_L001_R1_001.fastq.gz OK
F38SelR1_ACTTGA_L001_R2_001.fastq.gz OK
F38SelR1_ACTTGA_L002_R1_001.fastq.gz OK
F38SelR1_ACTTGA_L002_R2_001.fastq.gz OK
F38SelR2_GATCAG_L001_R1_001.fastq.gz OK
F38SelR2_GATCAG_L001_R2_001.fastq.gz OK
F38SelR2_GATCAG_L002_R1_001.fastq.gz OK
F38SelR2_GATCAG_L002_R2_001.fastq.gz OK
MGD2_SO_CAGATC_L005_R1_001.fastq.gz MGD2_SO_CAGATC_L001_R1_001.fastq.gz
MGD2_SO_CAGATC_L005_R2_001.fastq.gz MGD2_SO_CAGATC_L001_R2_001.fastq.gz
MGD2_SO_CAGATC_L006_R1_001.fastq.gz MGD2_SO_CAGATC_L002_R1_001.fastq.gz
MGD2_SO_CAGATC_L006_R2_001.fastq.gz MGD2_SO_CAGATC_L002_R2_001.fastq.gz
MGD_SO_CAGATC_L005_R1_001.fastq.gz MGD_SO_CAGATC_L001_R1_001.fastq.gz
MGD_SO_CAGATC_L005_R2_001.fastq.gz MGD_SO_CAGATC_L001_R2_001.fastq.gz
MGD_SO_CAGATC_L006_R1_001.fastq.gz MGD_SO_CAGATC_L002_R1_001.fastq.gz
MGD_SO_CAGATC_L006_R2_001.fastq.gz MGD_SO_CAGATC_L002_R2_001.fastq.gz
Sel_R1_F77_TTAGGC_L003_R1_001.fastq.gz F77SelR1_TTAGGC_L001_R1_001.fastq.gz
Sel_R1_F77_TTAGGC_L003_R2_001.fastq.gz F77SelR1_TTAGGC_L001_R2_001.fastq.gz
Sel_R1_F77_TTAGGC_L004_R1_001.fastq.gz F77SelR1_TTAGGC_L002_R1_001.fastq.gz
Sel_R1_F77_TTAGGC_L004_R2_001.fastq.gz F77SelR1_TTAGGC_L002_R2_001.fastq.gz
Sel_R2_F77_GATCAG_L003_R1_001.fastq.gz F77SelR2_GATCAG_L001_R1_001.fastq.gz
Sel_R2_F77_GATCAG_L003_R2_001.fastq.gz F77SelR2_GATCAG_L001_R2_001.fastq.gz
Sel_R2_F77_GATCAG_L004_R1_001.fastq.gz F77SelR2_GATCAG_L002_R1_001.fastq.gz
Sel_R2_F77_GATCAG_L004_R2_001.fastq.gz F77SelR2_GATCAG_L002_R2_001.fastq.gz
```
Step 5) mkdir / redefine Def_dir



