## Analyzing Sequence Data from Pooled Illumina Sequences

The structure of this data set gives an overview of the process to work with the raw outputs from Illumina sequencing of Pooled sequence data (in .fastq.gz format). The Data set I am using is from .......

DNA extraction of populations occurred at MSU by Dr. Michael DeNieu and Mauricio Losilla using the Zymo DNA extractor for insects with population pool sizes of 60 flies. Samples were submitted to RTSF (Research Technology Support Facility) Genomics Core (MSU), where next-generation sequencing libraries were prepared (Illumina TruSeq Nano DNA Library prep kit) and samples were loaded onto two lanes of an Illumina HiSeq 2500 rapid flow cell (version 1). Sequencing was done with Rapid SBS reagents in a 2x150bp paired end format, with the bases called using Illumina Real Time Analysis RTA v1.18.61, and converted to FastQ format. 

A total of 13 sequences were run through Illumina sequencing, the shared ancestor and the 4 treatments/replicates at generation 38, generation 77 and generation 115.


### Initial Analysis: Quality Control

The first step is to insepect the raw sequence reads available using md5sum and Fastqc





