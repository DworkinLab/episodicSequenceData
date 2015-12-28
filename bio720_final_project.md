#Paul Knoops
##Bio 720 Final Project
## Determining patterns in alleleic variation for experimentally evolved Drosophila melanogaster under predation threat
_____________________________________________________________________________
## Introduction

###The Question?

One of the major forces in evolution is natural selection, with beneficial traits for survival persisting and adapting species based on the selective pressures they face. However, natural selection is functionally not well understood, as studies into speciation and evolution can generally only focus on the final product of selection and only infer natural histories and related phylogenies to study how organisms evolve due to selection. Using a combination of whole genome sequencing and laboratory evolved population under specific selective pressure, we can now study organisms adaptation to new environments and the underlying genetic evolution (Tobler *et al* 2014, Turner *et al* 2012). These Evolve and Resequence (E&R) studies using sequenced genomes of pooled populations (Pool-Seq), genome-wide allele frequencies between selection regiems can be found, and by comparing different generations, patterns of allelic variation and trajectories of allele frequencies (Orozco-terWengel *et al.* 2012). By looking at 4 generations of pooled *Drosophila melanogaster* sequences, raised under different predatory selection, we can hope to better understand how natural selection effects genetic adaptations from a polymorphic base population.

By utilizing experimental evolution (EE), observations can be made on functional evolution of species within a controlled environment. Using laboratory experimentally evolved populations, evolutionary histories can be studied from a base population evolving due to varying selective pressures manipulated by the experimenters. Experimental evolution usually falls under two catagoiries based on the base population. These categories involve either using a base population of genetically invariable population allow adaptation through the accumulation of beneficial mutations for survival or by starting with a polymorphic population with high variability, which may already have beneficial alleles for selection to act upon (Schlotterer *et al.* 2014). Generally, using invariable populations are limited to microorganisms and larger, polymorphic eukaryotes are studied as a variable base population if they have short generation time and are easy to work with in the lab. This second approach, in connection with the resequencing methods, can allow researchers to see how alleles present in a population can potentially evolve under specific selective regiems, and how natural selection effects alleles to increase or decrease in frequency. Allowing natural selection to act on these populations allows for independent evolution, not restricted by researcher bias for particular behaviours or traits, allowing for any number of possible outcomes, and not restricted to a narrow scope of experimentation. Experimental evolution has been used for many years to discover the evolutionary outcomes for differences in morphology and behaviours, but now with the advancement of next-generation sequencing for a more affordable cost, the evolution of genes and alleles can now be studied. Using bioinformatic tools, the evolutionary dynamics of selected alleles can be observed for sequence data from selected populations

In order to elucidate population-wide trends for evolution of selected groups, a large number of individual sequencing and comparisons must be done. However, this is not a feasible method with current practices in genome sequencing, and an alternative method can be used. The most feasible option is to sequence multiple indviduals (a pool from the population) together in a method known as Pool-Seq. This Pool-Seq data is able to show variation between populations and differences in genome allele frequency estimates (reviewed in Rellstab et al., 2013; Schlötterer et al., 2014). The mixture of E&R with Pool-Seq analysis is relatively new compared to other E&R methods, and has few programs to analyze Pool-Seq data, but researchers working with this data have produced two notable software packages; PoPoolation and PoPoolation2 (Kofler *et al.* 2011a, Kofler *et al.* 2011b). PoPoolation was released as a method for analysis of population genetic inferences of a single population of Pool-Seq data, with PoPoolation2 developed as a comparative software for populations. When evaluating populations evolved under different selective pressures, PoPoolation2 is useful to compare different selection regiems and also different generational time points of the selected individuals. 

Natural selection allows beneficial traits for survival to evolve and persist in adapting populations based on the different selective pressures they face (Darwin 1859). Using experimental evolution and Pool-Seq analysis, natural selection can be studied for populations evolving populations at a genetic level. One strong selective pressure on evolving species is predation, where traits favourable for predator avoidance and escape are selected for and survivors pass genes on to further generations (Lima and Dill 1990). The evolutionary outcomes of predator selection can be seen within many species through prey species crypsis, mimicry, startle behaviours, or physiological defences. However, how selection has evolved populations to display complex behaviours and morphologies was only be inferred based on natural histories, but can not be functionally studied within the laboratory (Lima and Bednekoff 1999). Using *Drosophila melanogaster* sequence data for populations evolved under variable predation threat, genomic evolution of anti-predatory behaviour will be studied. 

Experimental evolution studies on polymorphic multicellular organisms require both short generation time and can be easily studied in a laboratory setting, along with a strong response to the selective regiem that will be used (Schlotterer *et al.* 2014). *Drosophila melanogaster (D. mel)* is a perfect fit for these studies, as it is one of the most commonly studied organisms in the lab, with well understood and widely studied genomics, behaviour, and development (Hall 1994, Sokolowski 2001). Populations of *D. mel* within the lab have been under predatory selection for 122 generations, with exposure to predators for 24 hours once per generation. The "episodic" populations experience a high mortality for a brief exposure to predators, and are expected to have adapted to this pressure. Only those that are able to survive the high predation threat are able to lay eggs and pass on their genes to the next generation. Natural selection should favour those that are able to survive, and over many generations, those survivors should have beneficial alleles increase in frequency, and alleles detrimental to survival should decrease. These populations either experience predation from the ambush hunting Chinese praying mantids (*Tenodera aridifolia sinensis*) for the 24 hours (predation treatment) or do not have the threat of predatory mortality at all (control treatment). Both treatments have 2 replicates, which can allow for independent evolution and different evolutionary trajectories. One advantage for experimental evolution with natural selection is that evolution is independent, where these populations can evolve along any trajectories and not on any linear, predetermined path. Although we can make predictions on what genes or alleles may be increasing or decreasing, the randomness that is evolution makes it difficult to know what to expect from these populations, and replicates can allow for these different adaptive paths to evolve (DeNieu *et al.* 2014). 

By comparing the two treatments, observations on differences in gene evolution due to selection can be studied. By following the methods outlined for PoPoolation2, along with papers using these pipelines (Kofler *et al* 2011b), comparisons between 4 generations of the evolved populations can be studied. Using this information, we can better understand not only how predation may shape prey species, but how natural selection acts on adapting genomes and their allele freguencies. 



##Methods/Pipeline
###1. Sequencing and populations

DNA extraction was performed Dr. Michael DeNieu and Dr. Sudarshan Chari at Michigan State University using the Zymo DNA extration for insects, and prepared NGS libraries with Illumina TruSeq Nano DNA Library preparation kit. Sequencing was done on two lanes of Illumina HiSeq 2500 Rapid Run flow call (v1.) at two times (2012 and 2015). In addition, generation 0 (base population) was sequenced twice, with the second providing greater coverage and higher accuracy for the sequence, which will be used together for most accurate full sequences. Many FastQ files of raw reads were made available and analyzed from the Golding server at McMaster University. 


###2. Quality Control and Fastqc

The first step is to check if the data uploaded correctly to my home directory using md5sum. A matched set using md5sum on md5.txt from raw data directory (with the raw sequence files and the md5.txt file in it) will output either "FAILED" or "OK" to know if the raw reads match the line in the md5.txt file, to signify a correct or incorrect transfer.
```
md5sum - c md5.txt
```
Flags;
"-c" = report if checksums match contents of files (OK).

Next check with Fastqc (Babraham Bioinformatics, 2015) for quality of the reads. Mistakes can occur with sequencing and using Fastqc allows one to view the quality of the reads that the sequencer has found. The process will output two files (*fastqc.html and *fastqc.zip). The *fastqc.html will be loaded to local machine and opened in web browser to view files. Moving to local machine also shown below while on the local machine.

```
fastqc -o /home/paul/episodicData/fastqcOutputs /home/paul/episodicData/rawData*.fastq.gz
```
Flags
"-o" sends all output files to output directory.
	
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/fastqcOutputs/*_fastqc.html /Users/paulknoops/episodicWork/data/fastqcOutputs
```

###3. Quality Trimming files (Trimmomatic)

Trimmomatic (Bolger *et al.* 2014) will be used to remove (trim) low quality reads or technical sequences (adapters) from the raw fastq files. Trimmomatic version 0.33 is used with the following outlines:

```
* IlluminaClip = adapter removal (TruSeq3-PE:2:30:10)
* LEADING & TRAILING = 3; removal at start end end if below quality
* MINLEN = minimum length of 36
* MAXINFO = adaptive quality (balance b/w length and quality) = 0.5. 0.5 is the intermediate between strict (0.8) and relaxed (0.2)
```

These will trim the data sequence to a higher quality of confidence that these are accurate sequences to be correctly mapped to the reference genome (step 4.)

```
#!/bin/bash

# make variable for trimmomatic program location
trim_dir=/usr/local/trimmomatic

# make input directory for raw reads
dir=/home/paul/episodicData/rawData/

# make output directory from trimmomatic outputs
out_dir=/home/paul/episodicData/trimmomaticOutputs

# make path to adapter sequences (to be used with ILLUMINACLIP)
adapter_path=/usr/local/trimmomatic/adapters

files=(${dir}*_R1_001.fastq.gz)
# echo ${files[@]}

for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`
java -jar ${trim_dir}/trimmomatic-0.33.jar PE -phred33 -trimlog ${out_dir}/trimlog.txt ${dir}${base}_R1_001.fastq.gz ${dir}${base}_R2_001.fastq.gz ${out_dir}/${base}_R1_PE_phred33.fastq.gz ${out_dir}/${base}_R1_SE_phred33.fastq.gz ${out_dir}/${base}_R2_PE_phred33.fastq.gz ${out_dir}/${base}_R2_SE_phred33.fastq.gz ILLUMINACLIP:${adapter_path}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done
```

Fastqc can also be done for the trimmed files to compare, but is not necessary

```
fastqc -o /home/paul/episodicData/fastqcTrimOutputs /home/paul/episodicData/trimmomaticOutputs/*.fastq.gz 
```

	
###4. Bring in reference sequence 

Using Drosophila (version r5.57.fasta.gz) uploaded from flybase.com using curl from index sequence directory.

```
index_dir=/home/paul/episodicData/indexSequence
cd ${index_dir}
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz
```
Set bwa index genome 

```
bwa index dmel-all-chromosome-r5.57.fasta.gz
```	

###5. BWA Mapping
Aligning the reads to the known reference genome will take all the trimmed raw reads and piece them together using BWA (Burrows-Wheeler Alignment, Li and Durbin 2010). The aligned output files are multiple SAM files for each treatment. While reviewing the outputs, raw reads from the same pooled populations are separated based on lanes, and have not been fully aligned. To account for this, the two lanes will be merged (after converted sam to bam files; see step 7.).
	
```
#!/bin/bash

# Log onto remote server
# make BWA directory path
bwa_dir=/usr/local/bwa/0.7.8

cd ${bwa_dir}

# make variable for working directory for
dir=/home/paul/episodicData/trimmomaticOutputs

#make variable for reference genome
ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz

# make variable for output directory
sam_dir=/home/paul/episodicData/mappedSequence/SAM_files

#make an array for each file in the directory "dir" that ends in _R1_PE_phred33.fastq.gz
files=(${dir}/*_R1_PE_phred33.fastq.gz)

#Check with echo
#echo ${files[1]}
#echo ${files[@]}


#Use "for loop" to map reads with the same "basename" to ref_genome
#Two flags for bwa mem
# -t = number of processors
# -M	Mark shorter split hits as secondary (for Picard compatibility(see step 8).
# Do I need qstat

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE_phred33.fastq.gz`
bwa mem -t 8 -M ${ref_genome} ${dir}/${base}_R1_PE_phred33.fastq.gz ${dir}/${base}_R2_PE_phred33.fastq.gz > ${sam_dir}/${base}_aligned_pe.SAM
done
```



###6. SAM to BAM
Sequence Alignment/Mapped files (sam) (Li *et al.* 2009) need to be converted to the easier to use bam files (Binary Alignment/Mapped). bam format will allow for merging files and creating the mpileup file needed to ran PoPoolation2 program and scripts. This is done with samtools view (Li *et al.* 2009) with the -q flag, allowing some quality control removing improperly mapped reads below a score of 20 (as recommended by PoPoolation2 creators). The other flags indicate sam (-s) or bam (-b) format. This is piped to samtools sort, ordering the file by the leftmost coordinate for ease of use further down the pipeline. The final output will be a bam file matching that of the sam file.

```
#! /bin/bash
sam_dir=/home/paul/episodicData/mappedSequence/SAM_files/
bam_dir=/home/paul/episodicData/mappedSequence/BAM_files/
files=(${sam_dir}*.SAM)
echo ${files[@]}
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .SAM`
samtools view -b -S -q 20 ${sam_dir}${base}.SAM | samtools sort - ${bam_dir}${base}
done
```

###7. Merge different lanes but same pooled sequence populations: samtools merge

Based on the mapping parameters, the output data from one pooled population (treatment/generation) are still separated by the lanes of Illumina. 

```
Con_R1_F77_ATGTCA_L001_aligned_pe.bam  F38ConR2_TTAGGC_L001_aligned_pe.bam
Con_R1_F77_ATGTCA_L002_aligned_pe.bam  F38ConR2_TTAGGC_L002_aligned_pe.bam
Con_R2_F77_ATTCCT_L001_aligned_pe.bam  F38SelR1_ACTTGA_L001_aligned_pe.bam
Con_R2_F77_ATTCCT_L002_aligned_pe.bam  F38SelR1_ACTTGA_L002_aligned_pe.bam
F115ConR1_TAGCTT_L001_aligned_pe.bam   F38SelR2_GATCAG_L001_aligned_pe.bam
F115ConR1_TAGCTT_L002_aligned_pe.bam   F38SelR2_GATCAG_L002_aligned_pe.bam
F115ConR2_GGCTAC_L001_aligned_pe.bam   MGD2_SO_CAGATC_L005_aligned_pe.bam
F115ConR2_GGCTAC_L002_aligned_pe.bam   MGD2_SO_CAGATC_L006_aligned_pe.bam
F115SelR1_GTTTCG_L001_aligned_pe.bam   MGD_SO_CAGATC_L005_aligned_pe.bam
F115SelR1_GTTTCG_L002_aligned_pe.bam   MGD_SO_CAGATC_L006_aligned_pe.bam
F115SelR2_GTGGCC_L001_aligned_pe.bam   Sel_R1_F77_TTAGGC_L003_aligned_pe.bam
F115SelR2_GTGGCC_L002_aligned_pe.bam   Sel_R1_F77_TTAGGC_L004_aligned_pe.bam
F38ConR1_ATCACG_L001_aligned_pe.bam    Sel_R2_F77_GATCAG_L003_aligned_pe.bam
F38ConR1_ATCACG_L002_aligned_pe.bam    Sel_R2_F77_GATCAG_L004_aligned_pe.bam
```

In order to join these two into one singe read, samtools merge (Li *et al.* 2009) is used. The basic outline of samtools merge is "samtools merge --output --input1.bam --input2.bam .... inputN.bam" The input will be based on lane, however, lanes vary and are not constant (L001/L002). The options to solve this could be to run it multiple times (changing input each time), or change the names of files to match (i.e L003/L004 becomes L001/L002). The easiest method to simply change lane names to match L001/L002 using move function (mv X_L003_aligned_pe.bam X_L001_aligned_pe.bam etc.) than perform script below is used.

```
#!/bin/bash
#Only works if all lanes are L001/L002

bam_dir=/home/paul/episodicData/mappedSequence/BAM_files
merged_dir=/home/paul/episodicData/mappedSequence/merged_bam_files
files=(${bam_dir}/*_L001_aligned_pe.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.bam`
samtools merge ${merged_dir}/${base}_merged_aligned_pe.bam ${bam_dir}/${base}_L001_aligned_pe.bam ${bam_dir}/${base}_L002_aligned_pe.bam
done
```

*Alternatively, in order to merge lanes earlier, include second BWA mem (after initial mapping) for L001 and L002*

```#!/bin/bash

bwa_dir=/usr/local/bwa/0.7.8
cd ${bwa_dir}
sam_dir=/home/paul/episodicData/mappedSequence/SAM_files
ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
out_dir=/home/paul/episodicData/mappedSequence/merge_SAM
files=(${sam_dir}/*_L001_aligned_pe.SAM)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.SAM`
bwa mem -t 8 -M ${ref_genome} ${sam_dir}/${base}_L001_aligned_pe.SAM ${sam_dir}/${base}_L002_aligned_pe.SAM > ${out_dir}/${base}_aligned_pe.SAM
done
```

*Side Note: PoPoolation2 manual and examples did not have the merged lanes to deal with, so this changed the pipeline followed slightly (sam-bam can be done later). Also, steps 8,9,10 follow that of the PoPoolation pipeline (the predecessor to PoPoolation2, looking at single populations of Pool-Seq data), and not used with PoPoolation2 manual and outlines. However, based on information (and code) obtained from one of the authors of PoPoolation2 (Dr. Schlotterer), the additional steps are recommended for some PoPoolation2 pipelines. Below I am including these steps, which serve as a precaution to ensure proper mapping quality and duplicates are removed.*


###8. Sort with Picard
In order to complete the next step (mark and remove duplicates) the file must be sorted by Picard (http://broadinstitute.github.io/picard/index.html) (as it will not accept samtools sorting). In order to remove any duplicated areas within the reads, Picard is used (step 9). This needs to be sorted properly for Picard to read, and the function SortSam from Picard Tools sorts the files properly. The flags include -Xmx2g, which allocated Java 2 Gb of memory, and SO, which is the sort order (in this case based on coordinate). Silencing VALIDATION_SRINGENCY stops Picard from reporting every issue that would ultimately be displayed that would not aid the final results for this analysis.

```
#! /bin/bash

merged_dir=/home/paul/episodicData/mappedSequence/merged_bam_files
sort_dir=/home/paul/episodicData/mappedSequence/sort_bam_files
files=(${merged_dir}/*)
pic=/usr/local/picard-tools-1.131/picard.jar
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
java -Xmx2g -jar ${pic} SortSam I= ${merged_dir}/${base}.bam O= ${sort_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate
done
```
*For output of sorted Picard; the control and selection treatments for generation 77 had a file size of zero. Appears the files were not sorted correctly and gave no output. Nothing was sorted and the problems with this is being investigated, however, I continued with the methods in order to get results, but will redo analysis with F77 once the issues have been resolved (would add the results to the mpileup file and sync file later down the pipeline).*


###9. Remove duplicates
As said above, duplicates should be removed from the sequences. Duplicates can come from two sources; during sequencing isolated clusters are identified as two when they are a single cluster (optical duplicates), or duplicated during PCR. Duplicates are removed with Picard Tools, using the files generated by sorting with Picard (Step 8). The flag "-M" creates an output file of statistics of duplicates found, and marking remove duplicates as true will get rid of any found duplicated regions.

```
#! /bin/bash

sort_dir=/home/paul/episodicData/mappedSequence/sort_bam_files
rmd_dir=/home/paul/episodicData/mappedSequence/rmd_bam_files
files=(${sort_dir}/*)
pic=/usr/local/picard-tools-1.131/picard.jar
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${sort_dir}/*.sort.bam O= ${rmd_dir}/${base}.rmd.sort.bam M= ${rmd_dir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```

MGD (Gen 0) did not output a file, so ran line below just for the one file.

```
java -Xmx2g -jar /usr/local/picard-tools-1.131/picard.jar MarkDuplicates I= /home/paul/episodicData/mappedSequence/sort_bam_files/MGD_2_SO_CAGATC_merged_aligned_pe.sort.bam O= /home/paul/episodicData/mappedSequence/rmd_bam_files/MGD_2_SO_CAGATC_merged_aligned_pe.rmd.sort.bam M= /home/paul/episodicData/mappedSequence/rmd_bam_files/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
```

But output for this file said:

```
Exception in thread "main" htsjdk.samtools.SAMException: Value was put into PairInfoMap more than once.
```

There may have been an error merging the two Generation 0 files so I ran SortSam for the two original Gen0 files merged, which have different degrees of coverage and were sequenced at different times. This was to see if errors merging files had any possible effect, and to have some files that would still work for analysis.

Picard Sort:

```
#! /bin/bash

merged_dir=/home/paul/episodicData/mappedSequence/merged_bam_files/seperate_merged_MGD
sort_dir=/home/paul/episodicData/mappedSequence/sort_bam_files/seperate_merged_MGD_sort
files=(${merged_dir}/*)
pic=/usr/local/picard-tools-1.131/picard.jar
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
java -Xmx2g -jar ${pic} SortSam I= ${merged_dir}/${base}.bam O= ${sort_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate
done
```

Remove duplicates:

```
#! /bin/bash

sort_dir=/home/paul/episodicData/mappedSequence/sort_bam_files/seperate_merged_MGD_sort
rmd_dir=/home/paul/episodicData/mappedSequence/rmd_bam_files/seperate_merged_MGD_rmd
files=(${sort_dir}/*)
pic=/usr/local/picard-tools-1.131/picard.jar
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${sort_dir}/${base}.sort.bam O= ${rmd_dir}/${base}.rmd.sort.bam M= ${rmd_dir}/dupstatGen0.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```

As these two work not merged, I will continue analysis with the two files, which should be have enough coverage separately, with merging used to create a more confident generation 0 sequence. 
*The merged file error will be investigated and a merged file will be done for reanalysis (with Gen77)*

###10. Remove low quality mapping
This step may be ambiguos, where the flag "-q" was already placed when converting sam to bam. However, other flags can be added (and can be moved above in future work) to remove any ambiguously mapped reads or any unmapped reads. The flag "-F 0x0004" will remove reads that are not mapped, and the final output will be the final bam files to be used to create the mpileup files and sync files.
	
```
#! /bin/bash

rmd_dir=/home/paul/episodicData/mappedSequence/rmd_bam_files
final_bam=/home/paul/episodicData/mappedSequence/final_bam_files
files=(${rmd_dir}/*.rmd.sort.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .rmd.sort.bam`
samtools view -q 20 -F 0x0004 -b ${rmd_dir}/${base}.rmd.sort.bam > ${final_bam}/${base}.final.bam
done
```


###11. Create mpileup file with samtools
In order to run PoPoolation2, a mpileup (multiple pileup) file format must be used. The file has information from each sample, including chromosome name and position, coverage, reference base, the coverage, and base mapping / quality numbers. The "read base column" holds information on if the region is a match, mismatch, indel or low quality to signify variance to the reference base. This can allow for variant calling and locating variants for evolved populations. For this analysis, the mpileup is used for creating a sync folder for PoPoolation2. The flags for samtools mpileup are; "-B" = disable BAQ (base alignment quality) computation (which stops probabilistic realignment for potential misaligned SNPs, but some false SNPs should not affect the results for the evolved pooled populations), "-f" = path to reference sequence,  "-6" =illumina v.1.3+ (Illumina 1.9 for this data, signifies encoding type).

```
#! /bin/bash

ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
final_bam = /home/paul/episodicData/mappedSequence/final_bam_files
mpileup_dir=/home/paul/episodicData/mappedSequence/
samtools mpileup -6 -B -f ${ref_genome} ${final_bam}/*.bam > ${mpileup_dir}episodicData.mpileup
```

###12. Convert to sync file 
For PoPoolation2 to run, the mpileup file must by synced with PoPoolation2 (Kofler *et al.* 2011b) to allow for quicker analysis of the file. The sync file will only need to be parsed once during further analysis to make the process faster. This is done with the PoPoolation2 script "mpileup2sync.jar" with defined input, output, quality and number of processors (threads). Fastq type is needed, based on the encoding of sequences (Illumina 1.9 for these sequences). The output file (.sync) has columns for each chromosomal position, with the reference character to compare to columns for each population (showing allele frequencies in format A:T:C:G:N:del). This will be the input for PoPoolation2 scripts below.

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

java -jar /usr/local/popoolation/mpileup2sync.jar --input ${map_dir}/episodicData.mpileup --output ${map_dir}/episodicData.sync --fastq-type illumina --min-qual 20 --threads 8
```

Upon inspection of the files, no allele frequencies are shown (A:T:C:G:N:del = 0:0:0:0:0:0), indicating some errors previously in the pipeline. Without data to find values; none of the tests outlined below can be analyzed.


### Alternative methods; testing where the problem may be

A thought is that using Picard (which is not included in some PoPoolation2 tutorials) may be interfering with the end results, so I will re-run the analysis from step 7, skipping strait to step 11 and create an mpileup file for them. 
*This is being done to hopefully have some initial results, with less stringency on the sequences to analyze, and to hopefully discover where the errors may be present*

From merged bam files, mpileup and sync are done:

```
#! /bin/bash

ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
merged_bam=/home/paul/episodicData/mappedSequence/merged_bam_files
mpileup_dir=/home/paul/episodicData/mappedSequence/
samtools mpileup -6 -B -f ${ref_genome} ${merged_bam}/*.bam > ${mpileup_dir}episodicData_noPicard.mpileup
```

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

java -jar /usr/local/popoolation/mpileup2sync.jar --input ${map_dir}/episodicData_noPicard.mpileup --output ${map_dir}/episodicData_noPicard.sync --fastq-type illumina --min-qual 20 --threads 8
```

Once again, no output allele frequencies (all 0's)

_______________________
Using non-merged files (test if merge is mistake)
Adding -Q = skip bases with base quality smaller than the given value (for 0 does not skip any low quality)

```
#! /bin/bash

ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
BAM_bam=/home/paul/episodicData/mappedSequence/BAM_files
mpileup_dir=/home/paul/episodicData/mappedSequence/
samtools mpileup -6 -B -Q 0 -f ${ref_genome} ${BAM_bam}/*.bam > ${mpileup_dir}episodicData_nomerge.mpileup
```
This step had an output (checked with less ${mpileup_dir}episodicData_nomerge.mpileup) of file containing many symbols for each sample combined. Output looks like this:

```
YHet	3	G	4	....	*+!$	4	....	**+*	1	.	)	2	..	()	2	..	$%	..	*%	2	..	'%	2	..	%%	3	...	)*#	5	.....	++!+%	5	.....	#&&$#	2	..	'%	0	*	*	0	*	*	0	*	*	5	....,	*%&&%	1	.	#	2	..	#%	1	.	%	4	....	+#**	2	..	(*	1	.	+	2	..	(*	1	..	*+	0	*	*	5	....,	+*++%
YHet	4	G	4	....	*+!$	4	....	+#)*	1	.	&	2	..	))	2	..	&#	..	*#	2	..	!%	2	..	*%	3	...	%*#	5	.....	*)!*%	5	.....	##&$%	2	..	'%	0	*	*	0	*	*	0	*	*	5	....,	*%(&%	1	.	%	2	..	%#	1	.	%	4	....	)#+*	2	..	(+	1	.	+	2	..	(+	1	..	++	1	^].	$	5	....,	++++%
YHet	5	T	4	....	!#!%	4	....	!%&#	1	.	!	2	..	!)	2	..	($	..	!%	2	..	"$	2	..	($	3	...	'(%	5	.....	!)#%$	5	.....	!&#!%	2	..	%!	0	*	*	0	*	*	0	*	*	5	....,	*%(%%	1	.	%	2	..	!$	1	.	"	4	....	!$'!	2	..	!'	1	.	'	2	..	!'	1	..	'(	1	.	$	5	....,	''&(%
YHet	6	C	4	....	')!%	4	....	!)$(	1	.	!	2	..	%'	2	..	""	..	%%	2	..	%$	2	..	)&	3	...	'*$	5	.....	$(()$	5	.....	"&!"$	2	..	)$	0	*	*	0	*	*	0	*	*	5	....,	+%)(%	1	.	%	2	..	!%	1	.	$	4	....	&$&'	2	..	$)	1	.	)	2	..	$)	1	..	()	1	.	$	5	....,	)())%
YHet	7	A	4	....	)*!'	4	....	$(&)	1	.	!	2	..	''	2	..	!#	..	'%	2	..	'#	2	..	*#	3	...	&&#	5	.....	$))*#	5	.....	!&'!%	2	..	)#	0	*	*	0	*	*	0	*	*	5	....,	*%*)#	1	.	%	2	..	!%	1	.	%	4	....	'#))	2	..	%+	1	.	+	2	..	%+	1	..	)*	1	.	'	5	....,	*(**#
YHet	8	C	4	....	**''	4	....	(*)+	1	.	$	2	..	)'	2	..	'$	..	)%	2	..	*$	2	..	)$	3	...	$)"	5	.....	))*(%	5	.....	!&!!%	2	..	)$	0	*	*	0	*	*	0	*	*	5	....,	*%**%	1	.	%	2	..	$%	1	.	$	4	....	'$**	2	..	)+	1	.	+	2	..	)+	1	..	*+	1	.	'	5	....,	*)+*%
YHet	9	G	4	....	&+''	4	....	&)'+	1	.	%	2	..	!'	2	..	)%	..	(%	2	..	*%	2	..	*!	3	...	#*#	5	.....	!'*)%	5	.....	!$(!%	2	..	)#	0	*	*	0	*	*	0	*	*	5	....,	*%)#&	1	.	%	2	..	%%	1	.	%	4	....	'!**	2	..	*+	1	.	+	2	..	*+	1	..	**	1	.	'	5	....,	++++%
YHet	10	T	4	....	%(')	4	....	()&*	1	.	&	2	..	!&	2	..	*%	..	)%	2	..	*%	2	..	!!	3	...	!)!	5	.....	''$!%	5	.....	!#'!%	2	..	)#	0	*	*	0	*	*	0	*	*	5	....,	+%'(&	1	.	!	2	..	!#	1	.	%	4	....	!#**	2	..	!)	1	.	)	2	..	!)	1	..	%%	1	.	'	5	....,	')**%
```

And this matches the PoPoolation sample output for mpileup

```
YHet    4067    N       9       ttttTtttt       aaab_Za_b       2       tt      `b      2       tt      \b
YHet    4068    N       9       c$cccCcccc      a_a_a]a_`       2       cc      ab      2       cc      `b
YHet    4069    N       8       aaaAaaaa        ]aaa_^__        2       aa      ab      2       aa      \b
YHet    4070    N       8       a$a$aAaaaa      \a]^YX_a        2       a$a     ab      2       aa      Wb
```

And to re-run the sync file

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

java -ea -Xmx7g -jar /usr/local/popoolation/mpileup2sync.jar --input ${map_dir}/episodicData_nomerge.mpileup --output ${map_dir}/episodicData_nomerge.sync --fastq-type illumina --min-qual 20 --threads 8
```

However, this step only produced lines of null

```
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
null
```

*Was not difference with Illumina or Sanger as fastq-type (ran also and gave only outputs of 0:0:0:0:0:0)*

_______________________
Merge issue

Looks like may be a merge issue, merging with alternative methods (bwa mem with lanes) shown in step 7, than rerun mpileup and sync.
*Changed for BAM files (removed SAM files to save space)*

```
#! /bin/bash

bwa_dir=/usr/local/bwa/0.7.8
cd ${bwa_dir}
bam_dir=/home/paul/episodicData/mappedSequence/BAM_files
ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
out_dir=/home/paul/episodicData/mappedSequence/merge_BAM
files=(${bam_dir}/*_L001_aligned_pe.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.bam`
bwa mem -t 8 -M ${ref_genome} ${bam_dir}/${base}_L001_aligned_pe.bam ${bam_dir}/${base}_L002_aligned_pe.bam > ${out_dir}/${base}_aligned_pe.bam
done
```

Ran very quick, may not have worked (and alternative done below - -6 flag change) so will not continue this for now.

```
samtools mpileup -6 -B -Q 0 -f /home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz /home/paul/episodicData/mappedSequence/merge_BAM/*.bam > /home/paul/episodicData/mappedSequence/merge_BAM/episodicData_bwaMerge.mpileup
```

```
java -ea -Xmx7g -jar /usr/local/popoolation/mpileup2sync.jar --input /home/paul/episodicData/mappedSequence/merge_BAM/episodicData_bwaMerge.mpileup --output /home/paul/episodicData/mappedSequence/merge_BAM/episodicData_bwaMerge.sync --fastq-type illumina --min-qual 20 --threads 8
```


________________________



Could be by having .sync
Ran java with .txt


```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

java -ea -Xmx7g -jar /usr/local/popoolation/mpileup2sync.jar --input ${map_dir}/episodicData_nomerge.mpileup --output ${map_dir}/episodicData_nomerge_sync.txt --fastq-type illumina --min-qual 20 --threads 8
```

Nope
___________________________

Is the issue the -6 flag (illumina) when is actually sanger?

```
#! /bin/bash

ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
BAM_bam=/home/paul/episodicData/mappedSequence/BAM_files
mpileup_dir=/home/paul/episodicData/mappedSequence/
samtools mpileup -B -Q 0 -f ${ref_genome} ${BAM_bam}/*.bam > ${mpileup_dir}episodicData_nomerge_Sanger.mpileup 
```

Output:

```
YHet    1       A       3       ^].^].^].       GJC     4       ^].^].^].^].    IGJI    1       ^].     ;       2       ^].^].  IJ
      2       ^].^].  GE      2       ^].^].  GD      2       ^].^O.  FD      2       ^].^].  JF      3       ^].^].^].       HJA     5       ^].^].^].^].^]. JJHJD   5       ^].^].^].^].^]. CHBFC   2       ^].^].  FF      0       *       *       0       *       *       0       *       *       5       ^].^].^].^].^], JDEGD   1       ^].     D       2       ^[.^].  DD      1       ^].     D       4
       ^].^].^].^].    IDHI    2       ^].^].  JJ      1       ^].     H       2       ^].^].  JJ      1       ^].     H       1       ^].     I       2       ^].^].  IJ      0       *       *       5       ^].^].^].^].^], IJJJD
YHet    2       G       4       ...^].  IJFB    4       ....    GHIJ    1       .       C       2       ..      IJ      2       ..
      FD      2       ..      GD      2       ..      FD      2       ..      JF      3       ...     GJA     5       .....   IJGID   5       .....   BAD@D   2       ..      FD      0       *       *       0       *       *       0       *       *       5       ....,   GDGED   1       .       B       2       ..      DD      1       .       D       4       ....    JDHH    2       ..      HJ      1
       .       H       2       ..      HJ      1       .       H       1       .       J       2       ..      JJ      0       *       *       5       ....,   JIJJD
YHet    3       G       4       ....    IJ1C    4       ....    IIJI    1       .       H       2       ..      GH      2       ..
      CD      2       ..      ID      2       ..      FD      2       ..      DD      3       ...     HIB     5       .....   JJ;JD   5       .....   BEECB   2       ..      FD      0       *       *       0       *       *       0       *       *       5       ....,   IDEED   1       .       B       2       ..      BD      1       .       D       4       ....    JBII    2       ..      GI      1
       .       J       2       ..      GI      1       .       J       1       .       J       2       ..      IJ      0       *       *       5       ....,   JIJJD
YHet    4       G       4       ....    IJ:C    4       ....    JBHI    1       .       E       2       ..      HH      2       ..
      EB      2       ..      IB      2       ..      ?D      2       ..      ID      3       ...     DIB     5       .....   IH:ID   5       .....   BBECD   2       ..      FD      0       *       *       0       *       *       0       *       *       5       ....,   IDGED   1       .       D       2       ..      DB      1       .       D       4       ....    HBJI    2       ..      GJ      1
       .       J       2       ..      GJ      1       .       J       1       .       J       2       ..      JJ      1       ^].     C       5       ....,   JJJJD
YHet    5       T       4       ....    ?B1D    4       ....    )DEB    1       .       =       2       ..      ?H      2       ..
      GC      2       ..      ?D      2       ..      AC      2       ..      GC      3       ...     FGD     5       .....   :HBDC   5       .....   5EB@D   2       ..      D>      0       *       *       0       *       *       0       *       *       5       ....,   IDGDD   1       .       D       2       ..      <C      1       .       A       4       ....    ?CF?    2       ..      =F      1
       .       F       2       ..      =F      1       .       F       1       .       D       2       ..      FG      1       .       C       5       ....,   FFEGD

```

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

java -ea -Xmx7g -jar /usr/local/popoolation/mpileup2sync.jar --input ${map_dir}/episodicData_nomerge_Sanger.mpileup --output ${map_dir}/episodicData_nomerge_Sanger.sync --fastq-type sanger --min-qual 20 --threads 2
```

Output:

```
YHet    1       A       3:0:0:0:0:0     4:0:0:0:0:0     1:0:0:0:0:0     2:0:0:0:0:0     2:0:0:0:0:0     2:0:0:0:0:0     2:0:0:0:0:0
     2:0:0:0:0:0     3:0:0:0:0:0     5:0:0:0:0:0     5:0:0:0:0:0     2:0:0:0:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:0:0:0     5:0:0:0:0:0     1:0:0:0:0:0     2:0:0:0:0:0     1:0:0:0:0:0     4:0:0:0:0:0     2:0:0:0:0:0     1:0:0:0:0:0     2:0:0:0:0:0     1:0:0:0:0:0     1:0:0:0:0:0     2:0:0:0:0:0     0:0:0:0:0:0     5:0:0:0:0:0
YHet    2       G       0:0:0:4:0:0     0:0:0:4:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:2:0:0     0:0:0:2:0:0     0:0:0:2:0:0
     0:0:0:2:0:0     0:0:0:3:0:0     0:0:0:5:0:0     0:0:0:5:0:0     0:0:0:2:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:5:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:4:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:0:0:0     0:0:0:5:0:0
YHet    3       G       0:0:0:3:0:0     0:0:0:4:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:2:0:0     0:0:0:2:0:0     0:0:0:2:0:0
     0:0:0:2:0:0     0:0:0:3:0:0     0:0:0:5:0:0     0:0:0:5:0:0     0:0:0:2:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:5:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:4:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:0:0:0     0:0:0:5:0:0
YHet    4       G       0:0:0:4:0:0     0:0:0:4:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:2:0:0     0:0:0:2:0:0     0:0:0:2:0:0
     0:0:0:2:0:0     0:0:0:3:0:0     0:0:0:5:0:0     0:0:0:5:0:0     0:0:0:2:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:5:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:4:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:1:0:0     0:0:0:2:0:0     0:0:0:1:0:0     0:0:0:5:0:0
YHet    5       T       0:3:0:0:0:0     0:3:0:0:0:0     0:1:0:0:0:0     0:2:0:0:0:0     0:2:0:0:0:0     0:2:0:0:0:0     0:2:0:0:0:0
     0:2:0:0:0:0     0:3:0:0:0:0     0:5:0:0:0:0     0:5:0:0:0:0     0:2:0:0:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:0:0:0:0:0     0:5:0:0:0:0     0:1:0:0:0:0     0:2:0:0:0:0     0:1:0:0:0:0     0:4:0:0:0:0     0:2:0:0:0:0     0:1:0:0:0:0     0:2:0:0:0:0     0:1:0:0:0:0     0:1:0:0:0:0     0:2:0:0:0:0     0:1:0:0:0:0     0:5:0:0:0:0
```

__________________
###*Will rerun final bam (with picard sort etc. for mpileip etc.) and all below (Step 13 and onward; Fst and CMH etc.) will be run with sync files = /home/paul/episodicData/mappedSequence/episodicData_Sanger.sync found below*

```
#! /bin/bash

ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
final_bam=/home/paul/episodicData/mappedSequence/final_bam_files
mpileup_dir=/home/paul/episodicData/mappedSequence/
samtools mpileup -B -Q 0 -f ${ref_genome} ${final_bam}/*.bam > ${mpileup_dir}episodicData_Sanger.mpileup
```

Output: Worked, so issue may have just been the -6
```
YHet    1       A       4       ^].^].^].^].    GEGD    4       ^].^O.^].^].    FDJF    8       ^].^].^].^].^].^].^].^].        HJAJJHJD        7       ^].^].^].^].^].^].^].   CHBFCFF 0       *       *       5       ^].^].^].^].^], JDEGD   3       ^].^[.^].       DDD     5       ^].^].^].^].^]. DIDHI   3       ^].^].^].       JJH     3       ^].^].^].       JJH
YHet    2       G       4       ....    FDGD    4       ....    FDJF    8       ........        GJAIJGID        7       ....... BAD@DFD 0       *       *       5       ....,   GDGED   3       ...     BDD     5       .....   DJDHH   3       ...     HJH     3       ...     HJH
YHet    3       G       4       ....    CDID    4       ....    FDDD    8       ........        HIBJJ;JD        7       ....... BEECBFD 0       *       *       5       ....,   IDEED   3       ...     BBD     5       .....   DJBII   3       ...     GIJ     3       ...     GIJ
YHet    4       G       4       ....    EBIB    4       ....    ?DID    8       ........        DIBIH:ID        7       ....... BBECDFD 0       *       *       5       ....,   IDGED   3       ...     DDB     5       .....   DHBJI   3       ...     GJJ     3       ...     GJJ
YHet    5       T       4       ....    GC?D    4       ....    ACGC    8       ........        FGD:HBDC        7       ....... 5EB@DD> 0       *       *       5       ....,   IDGDD   3       ...     D<C     5       .....   A?CF?   3       ...     =FF     3       ...     =FF
```

sync:

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

java -ea -Xmx7g -jar /usr/local/popoolation/mpileup2sync.jar --input ${map_dir}/episodicData_Sanger.mpileup --output ${map_dir}/episodicData_Sanger.sync --fastq-type sanger --min-qual 20 --threads 2
```

Output:

```
YHet    1       A       4:0:0:0:0:0     4:0:0:0:0:0     8:0:0:0:0:0     7:0:0:0:0:0     0:0:0:0:0:0     5:0:0:0:0:0     3:0:0:0:0:0     5:0:0:0:0:0     3:0:0:0:0:0     3:0:0:0:0:0
YHet    2       G       0:0:0:4:0:0     0:0:0:4:0:0     0:0:0:8:0:0     0:0:0:7:0:0     0:0:0:0:0:0     0:0:0:5:0:0     0:0:0:3:0:0     0:0:0:5:0:0     0:0:0:3:0:0     0:0:0:3:0:0
YHet    3       G       0:0:0:4:0:0     0:0:0:4:0:0     0:0:0:8:0:0     0:0:0:7:0:0     0:0:0:0:0:0     0:0:0:5:0:0     0:0:0:3:0:0     0:0:0:5:0:0     0:0:0:3:0:0     0:0:0:3:0:0
YHet    4       G       0:0:0:4:0:0     0:0:0:4:0:0     0:0:0:8:0:0     0:0:0:7:0:0     0:0:0:0:0:0     0:0:0:5:0:0     0:0:0:3:0:0     0:0:0:5:0:0     0:0:0:3:0:0     0:0:0:3:0:0
YHet    5       T       0:4:0:0:0:0     0:4:0:0:0:0     0:8:0:0:0:0     0:7:0:0:0:0     0:0:0:0:0:0     0:5:0:0:0:0     0:3:0:0:0:0     0:5:0:0:0:0     0:3:0:0:0:0     0:3:0:0:0:0
```


*NOTE* The output below are only based on numbers (i.e compare 1 to 2 etc.) and the numbers are the input order from the final bam files directory

```
1 = F115ConR1_TAGCTT_merged_aligned_pe.final.bam
2 = F115ConR2_GGCTAC_merged_aligned_pe.final.bam
3 = F115SelR1_GTTTCG_merged_aligned_pe.final.bam
4 = F115SelR2_GTGGCC_merged_aligned_pe.final.bam
5 = F38ConR1_ATCACG_merged_aligned_pe.final.bam
6 = F38ConR2_TTAGGC_merged_aligned_pe.final.bam
7 = F38SelR1_ACTTGA_merged_aligned_pe.final.bam
8 = F38SelR2_GATCAG_merged_aligned_pe.final.bam
9 = MGD2_SO_CAGATC_merged_aligned_pe.final.bam
10 = MGD_SO_CAGATC_merged_aligned_pe.final.bam
```

		
###13. Can now run scripts from popoolation2 directory to find / visualize Fst, CMH Tests, and Fisher's Exact Test
____________________________________________________________________

###Fst
Using PoPoolation2 script *fst-test.pl* all pairwise comparisons can be tested for a window size of 500 based of the sequence. Fst values show the genetic differentiation between populations. Each comparison (i.e F0 vs. Control Replicate 1 F115 or Selection Replicate 2 F38 vs. Control Replicate 1 F115) will be calculated, with the expectation that the later generations for selected populations will show the greatest differentiation at certain sites when compared to all others, but the intermediate generations (F38 and F77) show increases as well at the same sites. This can give approximate locations of areas of selection on the genome and if the common trend is seen as increasing differentiation (i.e adaptation) with increasing generations, those locations are strong candidates for genes driven by natural selection. Another possible outcome is to see a quick increase in frequency of an allele, but a plateau in later generations (i.e increase up to F77, than flat line after to F115) as seen by Orozco-terWengel (2012), which may be due to alleles reaching the optimum frequency. Depending on how alleles respond can tell us more on the function and process of selection on genes

The following parameters are used for this script;

```
Using input     /home/paul/episodicData/mappedSequence/episodicData_Sanger.sync
Using output    /home/paul/episodicData/mappedSequence/fst_Sanger.txt
Using min-count 3
Using min-coverage      10
Using max-coverage      250
Using min-covered-fraction      1.0
Using pool-size 60
Using window-size       500
Using step-size 500
Using asympt-unbiased   0
Using suppress na       1
Using test      0
Using help      0
```


```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

perl /usr/local/popoolation/fst-sliding.pl --window-size 500 --step-size 500 --suppress-noninformative --input ${map_dir}/episodicData_Sanger.sync --min-covered-fraction 1.0 --min-coverage 10 --max-coverage 250 --min-count 3 --output ${map_dir}/fst_Sanger.txt --pool-size 60
```

In order to view these results graphically, the Fst.txt file (with PoPoolation2 script *pwc2igv.pl*) will be loaded to IGV (Integrative Genomics Viewer).

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

perl /usr/local/popoolation/export/pwc2igv.pl --input ${map_dir}/fst_Sanger.txt --output ${map_dir}/fst_Sanger.igv 
```

This output file is moved to my local machine 

```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mappedSequence/fst_Sanger.igv /Users/paulknoops/episodicWork
```

and IGV is opened with command:

```
java -Xmx750m -jar igv.jar
```
 
and the file is opened in GUI format to open file and commpared using reference *D. mel* sequence r5.57

The output image for full pairwise comparisons;


*could not find way to fit all in yet*




###The CMH test (Cochran-Mantel-Haenszel test) 
This test can test the statistical significance between groups, depending on the input. The CMH test only tests significance of allele frequency changes between generations, with each populations only in the input once (may need to run multiple times dependent on the desired data). This test can identify SNPs with allele frequency changes among different time points. This is run with PoPoolation2 script *cmh-test.pl*, and diplayed to IGV with *cmh2gwas.pl*

--population 1-2,3-4; this will change based on what is being tested, but each number can only be included once. The input will depend on the .sync file based on population placement, and can change based on the input order. In this case, 1-2 compared the control treatments for F115 (R1 and R2) and 3-4 does the same for selection treatment, and continues like this (5-6, 7-8), 9-10 would only compare the two same unmerged reference sequences. 

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence
sync_file=/home/paul/episodicData/mappedSequence/episodicData_Sanger.sync

perl /usr/local/popoolation/cmh-test.pl --min-count 3 --min-coverage 10 --max-coverage 250 --population 1-2,3-4,5-6,7-8 --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/cmhtest_Sanger.txt

perl /usr/local/popoolation/export/cmh2gwas.pl --input ${map_dir}/cmhtest_Sanger.txt --output ${map_dir}/cmh_Sanger.gwas --min-pvalue 1.0e-20
```

Move to local machine (scp) and open IGV GUI format

```
java -Xmx2g -jar /Users/paulknoops/episodicWork/IGV_2.3.67/igv.jar
```

Other comparisons that can be made, and for this analysis, I will just compare generation 0 with generation 115 (both Selection and Control). Comparisons between replicates are completed, and below are inputs comparing gen 0 (using MGD2_SO..., which was seqeunced with more coverage) and generation 115. When comparing all different comparisons, a script will be set up to analyze all necessary input options.

Comparison between Control and Selection (R1 with R1)
```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence
sync_file=/home/paul/episodicData/mappedSequence/episodicData_Sanger.sync

perl /usr/local/popoolation/cmh-test.pl --min-count 3 --min-coverage 10 --max-coverage 250 --population 1-3,2-4 --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/cmhtest_115_R1.txt

perl /usr/local/popoolation/export/cmh2gwas.pl --input ${map_dir}/cmhtest_115_R1.txt --output ${map_dir}/cmh_115_R1.gwas --min-pvalue 1.0e-20
```

Comparison between Control and Selection (R1 with R2)
```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence
sync_file=/home/paul/episodicData/mappedSequence/episodicData_Sanger.sync

perl /usr/local/popoolation/cmh-test.pl --min-count 3 --min-coverage 10 --max-coverage 250 --population 1-4,2-3 --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/cmhtest_115_R1vsR2.txt

perl /usr/local/popoolation/export/cmh2gwas.pl --input ${map_dir}/cmhtest_115_R1vsR2.txt --output ${map_dir}/cmh_115_R1vsR2.gwas --min-pvalue 1.0e-20
```

Comparsion between Gen0 and each 115 treatment

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence
sync_file=/home/paul/episodicData/mappedSequence/episodicData_Sanger.sync

#F115ConR1 vs. MGD2

perl /usr/local/popoolation/cmh-test.pl --min-count 3 --min-coverage 10 --max-coverage 250 --population 1-9 --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/cmhtest_19.txt

perl /usr/local/popoolation/export/cmh2gwas.pl --input ${map_dir}/cmhtest_19.txt --output ${map_dir}/cmh_19.gwas --min-pvalue 1.0e-20

#F115ConR2 vs. MGD2

perl /usr/local/popoolation/cmh-test.pl --min-count 3 --min-coverage 10 --max-coverage 250 --population 2-9 --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/cmhtest_29.txt

perl /usr/local/popoolation/export/cmh2gwas.pl --input ${map_dir}/cmhtest_29.txt --output ${map_dir}/cmh_29.gwas --min-pvalue 1.0e-20

#F115SelR1 vs. MGD2

perl /usr/local/popoolation/cmh-test.pl --min-count 3 --min-coverage 10 --max-coverage 250 --population 3-9 --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/cmhtest_39.txt

perl /usr/local/popoolation/export/cmh2gwas.pl --input ${map_dir}/cmhtest_39.txt --output ${map_dir}/cmh_39.gwas --min-pvalue 1.0e-20

#F115SelR2 vs. MGD2

perl /usr/local/popoolation/cmh-test.pl --min-count 3 --min-coverage 10 --max-coverage 250 --population 4-9 --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/cmhtest_49.txt

perl /usr/local/popoolation/export/cmh2gwas.pl --input ${map_dir}/cmhtest_49.txt --output ${map_dir}/cmh_49.gwas --min-pvalue 1.0e-20
```



			
###Fisher's Exact Test 
This test gives the statistical significance of allele frequency differences using script *fisher-test.pl* and viewed using *pwc2igv.pl* and IGV.

```
#! /bin/bash

map_dir=/home/paul/episodicData/mappedSequence

perl /usr/local/popoolation/fisher-test.pl --input ${map_dir}/episodicData_Sanger.sync --output ${map_dir}/episodicData_Sanger.fet --min-count 3 --min-coverage 10 --max-coverage 250 --suppress-noninformative

perl /usr/local/popoolation/export/pwc2igv.pl --input ${map_dir}/episodicData_Sanger.fet --output ${map_dir}/episodicData_Sanger.fet.igv

# Load to IGV java -Xmx2g -jar /usr/local/igv/IGV_2.1.21/igv.jar
```

Output when run (and even when using help command for Popoolation perl fisher-test.pl --help) appears to be a script issue, and fishers will not be included yet.

```
Can't locate Text/NSP/Measures/2D/Fisher/twotailed.pm in @INC (@INC contains: /usr/local/popoolation /usr/local/popoolation/Modules /usr/local/ensembl/current/bioperl-1.2.3 /usr/local/ensembl/current/ensembl/modules /usr/local/ensembl/current/ensembl-compara/modules /usr/local/ensembl/current/ensembl-variation/modules /usr/local/ensembl/current/ensembl-funcgen/modules /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5 .) at /usr/local/popoolation/Modules/FET.pm line 10.
BEGIN failed--compilation aborted at /usr/local/popoolation/Modules/FET.pm line 10.
Compilation failed in require at fisher-test.pl line 9.
BEGIN failed--compilation aborted at fisher-test.pl line 9.
```






## Work Cited

^"Babraham Bioinformatics - FastQC A Quality Control Tool for High Throughput Sequence Data." Babraham Bioinformatics - FastQC A Quality Control Tool for High Throughput Sequence Data. Web. 18 Dec. 2015.

^Bolger, A.M., Lohse, M., & Usadel, B. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. *Bioinformatics*, btu170.

^Darwin, C. 1859. *The origin of species*. John Murray, London.

^DeNieu, M., Pitchers, W. & Dworkin, I. 2014. Adaptation to a novel predator in Drosophila melanogaser: How well are we able to predict evolutionary responses? In revisions, doi: http://dx.doi.org/10.1101/005322. 

^Hall, J.C 1994. The mating of a fly. *Science*. 264:1702–1714.

^Kofler R, Orozco-terWengel P, De Maio N, Pandey R.V, Nolte V, Futschik A, Kosiol C, Schlotterer C. 2011a. PoPoolation: a toolbox for popula- tion genetic analysis of next generation sequencing data from pooled individuals. *PLoS One* 6:e15925.

^Kofler R, Pandey R.V, Schlotterer C. 2011b. PoPoolation2: identifying dif- ferentiation between populations using sequencing of pooled DNA samples (Pool-Seq). *Bioinformatics* 27:3435–3436.

^Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*. 25:2078–2079.

^Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler Transform. *Bioinformatics*, Epub. [PMID: 20080505]

^Lima, S.L. & Bednekoff, P.A. 1999. Temporal variation in danger drives antipredator behaviour: The predation risk allocation hypothesis. *The American Naturalist*. 153:649–659

^Lima, S.L. & Dill, L.M. 1990. Behavioural decisions made under the risk of predation: a review and prospectus. *Canadian Journal of Zoology*. 68, 619–640.

^Orozco-terWengel P, Kapun M, Nolte V, Kofler R, Flatt T, Schlotterer C. 2012. Adaptation of Drosophila to a novel laboratory environment reveals temporally heterogeneous trajectories of selected alleles. *Molecular Ecology*. 21: 4931–4941.

^Schlotterer C, Kofler R, Versace E, Tobler R, Franssen S.U. 2014. Combining experimental evolution with next-generation sequencing: a powerful tool to study adaptation from standing genetic variation. *Heredity* doi:10.1038/hdy.2014.86.

^Sokolowski, M.B. 2001. *Drosophila*: Genetics meets behaviour. *Nature Reviews Genetics*. 2:879–890.

^Tobler R, Franssen S.U, Kofler R, Orozco-terWengel P, Nolte V, Hermisson J, Schlotterer C. 2014. Massive habitat-specific genomic response in D. melanogaster populations during experimental evolution in hot and cold environments. *Molecular Biology and Evolution*. 31:364–375.

^Turner T.L., Miller P.M. 2012. Investigating natural variation in Drosophila courtship song by the evolve and resequence approach. *Genetics* 191:633–642.
