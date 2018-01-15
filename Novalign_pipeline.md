# Mapping with Novoalign and the subsequent steps to clean data:
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________

## Running Novoalign mapper

Starting with sequence files that have already been inspected with md5sum and fastqc and have been trimmed using trimmomatic (See other files for running trimmomatic)

Novoalign link tutorial: http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/basic-short-read-mapping/

### Need to create directory for project and to house mapping outputs

```
# Project Directory
mkdir novoalign 
```

```
# mapping outputs
mkdir novo_dir
```

### Find path to call Novoalign variable

```
#Variable for novoalign
novoalign=/usr/local/novoalign
```


### Make a scripts directory to house scripts be ran 

```
mkdir novo_scripts
```

### Novoindex reference

Need index dir for novoindex files

```
mkdir novo_index
```

The reference genome needs to be indexed for novoalign mapping (with novoindex)

Script: novo_index.sh

```
#! /bin/bash

#Create variable for location of reference genome (fasta vs. fasta.gz?)
ref_genome=/home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta

#Variable for project
project_dir=/home/paul/episodicData/novoalign

#Variable for novoalign
novoalign=/usr/local/novoalign

#Variable for output directory
novo_index=${project_dir}/novo_index

#Index the reference with novoindex

${novoalign}/novoindex ${novo_index}/dmel-all-chromosome-r5.57_2.nix  ${ref_genome}

```


### Unzip Files

Note: Compressed read files are not supported in unlicensed versions.

 - The unlicensed version of Novoalign (used here) does not support the zipped files, so need to unzip trimmomatic outputs

```
#From trim_dir

gunzip *.gz
```

### Novoalign 

Flags: 

- d -- Full pathname of indexed reference sequence from novoindex

- f -- Files containing the read sequences to be aligned  

- o -- Specifies output report format and options (SAM)  

- i ###,## -- Sets fragment orientation and approximate fragment length for proper pairs.
    ex. -i 250 50  Defaults to paired end Illumina or Mate Pair ABI with 250bp insert and 50bp standard deviation (possible check below)
     
     - 400, 100 found based on initial mapping with novoalign first run through
     
     - using 500, 150 (below)

### Checking insert Size

Using output from Picard Sort (if available) from Bowtie2 or BWA mem previously (before removing duplicates) use Picard CollectInsertSizeMetrics.jar to get summary statistics on file for insert size and other information

If other mappers not available, map using defaults or extreme insert / SD values (i.e. 0, 500) and then run this script on .bam files (needs to be sorted with Picard)

```
java -jar /usr/local/picard-tools-1.131/picard.jar CollectInsertSizeMetrics \
    I=/home/paul/episodicData/novoalign/novo_rmd/F115ConR1_TAGCTT_novo_merge_novo_rmd.bam \
    O=/home/paul/episodicData/novoalign/novo_rmd/insert_size_metrics.txt \
    H=/home/paul/episodicData/novoalign/novo_rmd/insert_size_histogram.pdf
```

Example output:

```
MEDIAN_INSERT_SIZE|MEDIAN_ABSOLUTE_DEVIATION|MIN_INSERT_SIZE|MAX_INSERT_SIZE|MEAN_INSERT_SIZE|STANDARD_DEVIATION|READ_PAIRS|PAIR_ORIENTATION|WIDTH_OF_10_PERCENT|WIDTH_OF_20_PERCENT|WIDTH_OF_30_PERCENT|WIDTH_OF_40_PERCENT|WIDTH_OF_50_PERCENT|WIDTH_OF_60_PERCENT|WIDTH_OF_70_PERCENT|WIDTH_OF_80_PERCENT|WIDTH_OF_90_PERCENT|WIDTH_OF_99_PERCENT|SAMPLE|LIBRARY|READ_GROUP

542|94|30|28389293|542.112611|156.897954|25505882|FR|35|69|107|145|189|241|305|391|533|893  
     
For generation 115 alone: mean = 542, SD = 156
```

### Running Novoalign (File By File)

This process will run one file at a time

The script: novo_map.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=/home/paul/episodicData/trim_dir

#Path to output directory for mapped files
novo_dir=${project_dir}/novo_dir

files=(${trim_dir}/*_R1_PE.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq`

${novoalign}/novoalign -d ${novo_index} \
    -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq \ 
    -i 500,150 -o SAM > ${novo_dir}/${base}_novo.sam

done
```

This takes a long time, as the unlicensed version can only uses 1 thread (100% computer)

*** From novoalign reference manual: -c 99 Sets the number of threads to be used. On licensed versions it defaults 
to the number of CPUs as reported by sysinfo(). On free version the option is disabled ***

### Running Novoalign (Running In Parallel)

A solution to run each file seperatly in a simple splitting method

*** alternative that may be an option: add the & after the code in the for Loop (before the done) and it should push that "for" to the background and run the next file in sequence *** 


__1) Make the script to make multiple scripts__

Make dir for all output scripts: 

```
mkdir split_mappingScripts
```

Script to create many scripts

- each different file normally looped through one by one above are put into a seperate script

Script: novo_map_scriptMaker.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=/home/paul/episodicData/trim_dir

#Path to output directory
novo_dir=${project_dir}/novo_dir

files=(${trim_dir}/*_R1_PE.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq`
echo "${novoalign}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 500,150 -o SAM > ${novo_dir}/${base}_novo.sam" > ./split_mappingScripts/${base}.sh

done
```

__2) Create script to call all and run in parallel (use "&" which puts job in background then multiple can run at a time)__

This creates a file that has all the scripts made in step 1) in a list with ''&'' at the end to run in parrallel

One method to run only half is make files/basename based on lane (i.e 001 or 002)

Script: novo_createParallel_scipt.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Variable for each script location
map_scripts=${project_dir}/novo_scripts/split_mappingScripts

#Variable for location of whole script (novo_scripts)
scripts=${project_dir}/novo_scripts


files=(${map_scripts}/*.sh)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sh`
echo "${map_scripts}/${base}.sh &" >> ${scripts}/novo_parallel_map.sh

done
```

__3) Change permissions and run novo_parallel_map.sh__

If needed based on the computer space available, change input parametes to run subsets on different days (one option above)

Run on screen

Screen can be named with -S (ex. screen -S IDENTIFIERTITLE)

Can save all outputs of screen using script (ex. script LOGTITLE.log) and finish script with "exit"

```
novo_parallel_map.sh
```

This should map each file seperate in unison


### Save space again with trimmed files: Rezip

Rezip files in trim_dir (saves space)

From the trim_dir:
```
gzip *.fastq
```

### Left with mapped sequences with Novoalign: can continue with cleaning the sequence data to final .bam files
____________________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________

## Cleaning the aligned Data

__Steps used for all mapping outputs: changing parameters for input/output directories__

### Change SAM files to BAM files: novo_samTobam.sh

-- Saves space (BAM files are binary compressed versions of SAM files)

-- need bam directory for .bam files (mkdir novo_bam)

Flags:

- b -- output is .bam

- S -- input is .sam

- q 20 -- quality mapping score of 20 (standard throughout all experiments)
     

Script: novo_samTobam.sh

```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_dir=${project_dir}/novo_dir

#Path to output directory
novo_bam=${project_dir}/novo_bam

files=(${novo_dir}/*.sam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sam`
samtools view -b -S -q 20 ${novo_dir}/${base}.sam | samtools sort -o ${novo_bam}/${base}.bam
done
```

### Merge 
```
mkdir novo_merge
```

no flags, just merging the two lanes of the illumina sequencing run

The script: novo_merge.sh
```    
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_merge=${project_dir}/novo_merge


files=(${novo_bam}/*_L001_novo.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_novo.bam`
samtools merge ${novo_merge}/${base}_novo_merge.bam ${novo_bam}/${base}_L001_novo.bam ${novo_bam}/${base}_L002_novo.bam
done
```

### Picard Sort 

Need to sord with Picard to mark and remove duplicates (needs to be sorted with Picard in order for downstream analysis)

Need a directory for outputs, and a temporary directory for space allocation
```
mkdir novo_pic

mkdir novo_tmp
```

Flags:

 - Xmx2g -- 2 Gb of memory allocated
 
 - Djava.io.tmpdir=${tmp} -- using my temporary directory due to errors in space allocation to avoid errors while running (not necessary but useful)

 - I -- input
 
 - O -- output

 - VALIDATION_STRINGENCY=SILENT -- stops Picard from reporting every issue that would ultimately be displayed

 - SO=coordinate -- sort order based on coordinate
  
 
Script: novo_picard_sort.sh 
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_merge=${project_dir}/novo_merge

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to output directory
novo_pic=${project_dir}/novo_pic

#Path to tmp
novo_tmp=${project_dir}/novo_tmp


files=(${novo_merge}/*.bam)
for file in ${files[@]}
do
name=${file}

base=`basename ${name} .bam`
java -Xmx2g -Djava.io.tmpdir=${novo_tmp} -jar ${pic} SortSam \
I= ${novo_merge}/${base}.bam \
O= ${novo_pic}/${base}_novo_sort.bam \
VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${novo_tmp}

done
```

### Remove Duplicates
```
mkdir novo_rmd
```
Flags: Similar to above

- Xmx2g -- ""
        
- MarkDuplicates -- ""
        
- I -- ""
        
- O -- ""
        
- M -- creates an output file of statistics of duplicates found

- VALIDATION_STRINGENCY=SILENT -- ""

- REMOVE_DUPLICATES= true -- get rid of any found duplicated regions


Script: novo_rmd.sh
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_pic=${project_dir}/novo_pic

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to output directory
novo_rmd=${project_dir}/novo_rmd

files=(${novo_pic}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _novo_sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${novo_pic}/${base}_novo_sort.bam O= ${novo_rmd}/${base}_novo_rmd.bam M= ${novo_rmd}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```

### More QC and make final bam files

```
mkdir novo_final
```

Flags:

- q 20 -- ""

- F 0x0004 -- remove any unmapped reads (hexidecimal value for unmapped = 0x0004)

- b -- ""

Script: novo_final_qc.sh
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_rmd=${project_dir}/novo_rmd

#Path to output directory
novo_final=${project_dir}/novo_final


files=(${novo_rmd}/*_novo_rmd.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _novo_rmd.bam`
samtools view -q 20 -F 0x0004 -b ${novo_rmd}/${base}_novo_rmd.bam > ${novo_final}/${base}_novo_final.bam
done
```

### Merge the 2 ancestors

Need to merge the base generation additionaly (two sequence runs for ancestor need to merge: MGD2 and MGD)

Script: 
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to final directory
novo_final=${project_dir}/novo_final

samtools merge ${novo_final}/MGD3_SO_CAGATC_novo_merge_novo_final.bam \
${novo_final}/MGD2_SO_CAGATC_novo_merge_novo_final.bam \
${novo_final}/MGD_SO_CAGATC_novo_merge_novo_final.bam

mkdir ${novo_final}/Anc_unmerged
mv ${novo_final}/ MGD2_SO_CAGATC_novo_merge_novo_final.bam ${novo_final}/Anc_unmerged
mv ${novo_final}/ MGD_SO_CAGATC_novo_merge_novo_final.bam ${novo_final}/Anc_unmerged
```

### Indel Realigner with GATK


__1) Need an unzipped copy of the reference genome__

Made a second copy and unzipping the second copy
```
gunzip dmel-all-chromosome-r5.57_2.fasta.gz
```

__2) make a gatk directory__
```
mkdir novo_GATK
```

__3) Set up reference genome__

a) Create .dict file for reference: This creates a dictionary file for the ref genome with a header but no sam records (the header is only sequence records)

script: novo_dict_index.sh
```
#! /bin/bash

pic=/usr/local/picard-tools-1.131/picard.jar
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

java -jar ${pic} CreateSequenceDictionary R=${ref_genome} O=${index_dir}/dmel-all-chromosome-r5.57_2.dict
```

b) Create a .fai to reference genome: 
```
samtools faidx dmel-all-chromosome-r5.57_2.fasta
```

__4) Need Read Group headers__

For GATK to run, read group headers are needed to read for Indel Realignment, __BUT__ the details are not necessary and don't need to be accurate, just need to be there.

Flags:

- RGID --Read Group Identifier; for Illumina, are composed using the flowcell + lane name and number [using Lanes L001_L002 for now]

- RGLB -- DNA Preperation Library Identifier [library1 as place holder]

- RGPL -- platform/technology used to produce the read [Illumina]

- RGPU -- Platform Unit; details on the sequencing unit (i.e run barcode) [None, used for practice]

- RGSM -- Sample [Using the basename which is each unique sequence]

Script: novo_readgroups.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to .bam files
novo_final=${project_dir}/novo_final

files=(${novo_final}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${pic} AddOrReplaceReadGroups I=${novo_final}/${base}.bam O=${novo_final}/${base}_RG.bam RGID=L001_L002 RGLB=library1 RGPL=illumina RGPU=None RGSM=${base}

done
```

__5) Index the read group Bam files__

Need to have the .bam files indexed prior to the indel realignment

Script: novo_indexBam.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .bam files
novo_final=${project_dir}/novo_final

files=(${novo_final}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`
samtools index ${novo_final}/${base}_RG.bam
done
```

__6) Run GATK indel realigner__

GATK indel realigner takes two steps, 1) target the indels to be raligned (.intervals file) and 2) realign the indels (.realigned files)

Script: novo_gatk.sh
```
#!/bin/bash

#Variable for project name (file name)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_final=${project_dir}/novo_final

#Path to output directory
novo_GATK=${project_dir}/novo_GATK

#Variable for reference genome (non-zipped)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar

files=(${novo_final}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

java -Xmx8g -jar ${gatk} -I ${novo_final}/${base}_RG.bam -R ${ref_genome} -T RealignerTargetCreator -o ${novo_GATK}/${base}.intervals

java -Xmx8g -jar ${gatk} -I ${novo_final}/${base}_RG.bam -R ${ref_genome} -T IndelRealigner -targetIntervals ${novo_GATK}/${base}.intervals -o ${novo_GATK}/${base}_realigned.bam

done
```

## Analysis of the variation

### Tajima's Pi

1) Create Directories

Will need two new directories, one for pileup files, and one for the output Pi files: Be sure they are in project dir

```
mkdir novo_pileup
mkdir novo_pi
```

2) Create pileup files of every .bam file

Each generated .bam file needs to be in mpileup format for the use in Popoolation Scripts:

Flags:

- B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment

- Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)

- f -- path to reference sequence
        
Script: novo_Pi_pileups.sh
```
#! /bin/bash

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_GATK

# Path to output novoalign pileup files
output=${project_dir}/novo_pileup

index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

files=(${input}/*_merge_novo_final_realigned.bam)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _merge_novo_final_realigned.bam`

samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_merge_novo_final_realigned.bam > ${output}/${base}.pileup

done
```

3) Run script to calcualte Tajima's Pi

Using the Variance-sliding.pl script from Popoolation1

Flags:

- input -- input pileup file

- output -- output file with Tajima's Pi calculated

- measure [pi] -- Options include Tajima's Pi or Wattersons Theta or Tajima's D along chromosomes using a sliding window approach

- window-size [10000] -- size of the sliding window 

- step-size [10000] -- how far to move along with chromosome (if step size smaller, windows will overlap)
 
- min-count [2] -- minimum allele count 

- min-coverage [4] -- minimum coverage (not important if subsampling done..)

- max-coverage [400] --maximum coverage

- min-qual [20] -- minimum base quality (already filtered for 20 multiple times)

- pool-size [120] -- number of chromosomes (So double the number of individuals per pool)

- fastq-type [sanger] -- depending on the encoding of the fastq files

- min-covered-fraction [0.5] -- minimum percentage of sites having sufficient coverage in the given window -- 0.5 from example


Script: novo_tajima_pi.sh

```
#! /bin/bash

# Path to PoPoolation1 (Currently in Paul's Home directory)
popoolation=/home/paul/popoolation_1.2.2

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_pileup

# Path to output Tajima Pi files
output=${project_dir}/novo_pi

files=(${input}/*.pileup)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} .pileup`

perl ${popoolation}/Variance-sliding.pl \
	--input ${input}/${base}.pileup \
	--output ${output}/${base}.pi \
	--measure pi \
	--window-size 10000 \
	--step-size 10000 \
	--min-count 2 \
	--min-coverage 4 \
	--max-coverage 400 \
	--min-qual 20 \
	--pool-size 120 \
	--fastq-type sanger \
	--snp-output ${output}/${base}.snps \
	--min-covered-fraction 0.5
done
```
Need to check the files were created properly (check that there are no empty files).


4) Bring to local machine
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_pi/*.pi /Users/paulknoops/Bioinformatics/R-projects_git/episodicSequenceData/R_scripts/Pi_Analysis_Novo
```

5) Function to create plots of tajima Pi data

R Function to run and create plots

See: R_scripts/Pi_Analysis_Novo/Pi_PlotFunction.R
```
Pi_PlotFunction <- function(x, y) {
  require(ggplot2)
  x2 <- gsub("\\_.*","",x)
  y2 <- y
  
  #bowtie.pi
  #novo.pi
  #Read in the data:
  Datt <- read.table(x)
  colnames(Datt) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  #Remove unnecessary regions: Not necessary based on later steps
  Datt$chr <- as.character(Datt$chr)
  Datt2 <- Datt

  #Remove "na" pi values
  Datt2 <- Datt2[-which(Datt2$Pi=="na"),]
  
  #Need the numbers for chromosomes for labelling and colours:
  DattX <- Datt2[which(Datt2$chr=="X"),]
  a <- dim(DattX)[1]
  DattX$number <- 1:a
  
  Datt2L <- Datt2[which(Datt2$chr=="2L"),]
  b <- dim(Datt2L)[1]
  Datt2L$number <- (a+1):(a+b)
  
  Datt2R <- Datt2[which(Datt2$chr=="2R"),]
  c <- dim(Datt2R)[1]
  Datt2R$number <- (a+b+1):(a+b+c)
  
  Datt3L <- Datt2[which(Datt2$chr=="3L"),]
  d <- dim(Datt3L)[1]
  Datt3L$number <- (a+b+c+1):(a+b+c+d)
  
  Datt3R <- Datt2[which(Datt2$chr=="3R"),]
  e <- dim(Datt3R)[1]
  Datt3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  Datt4 <- Datt2[which(Datt2$chr=="4"),]
  f <- dim(Datt4)[1]
  Datt4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  #Full data frame of necessary chromosomes
  DattFull <- rbind(DattX, Datt2L, Datt2R, Datt3L, Datt3R, Datt4)
  
  #Pi as numeric
  DattFull$Pi=as.numeric(levels(DattFull$Pi))[DattFull$Pi]
  
  #Title:
  z2 <- paste(x2, y2, sep="_")
  
  # The plots: 
  Pi_plot <- ggplot(DattFull, aes(x = number, y= Pi, colour = chr)) 
  
  Pi_plot_2 <- Pi_plot + 
    geom_point(size=0.3, show.legend = F) +
    scale_y_continuous(limits=c(0, 0.02), breaks=seq(0, 0.02, 0.005)) + 
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    theme(text = element_text(size=20), 
          axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) +
    ggtitle(z2)
  
  return(Pi_plot_2)
}
```

6) Run the function for each .pi file

```
Pi_PlotFunction('F115ConR1_TAGCTT_novo.pi', "Novoalign")
```














_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
## Combining files together for analysis


### Create mpileup

```
mkdir novo_mpileup
```

Flags:

- B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment

- Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)

- f -- path to reference sequence
    

Script: novo_mpileup.sh
```
#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_final=${project_dir}/novo_GATK

#Path to output directory
novo_mpileup=${project_dir}/novo_mpileup

#Variable for reference genome (non-indexed for novoAlign)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

samtools mpileup -B -Q 0 -f ${ref_genome} ${novo_final}/*.bam > ${novo_mpileup}/${project_name}.mpileup
```

### Create .sync file 

--use mpileup dir

Flags:

- Xmx7g -- ""

- input -- ""

- output -- ""

- fastq-type -- needed for base encoding

- min-qual 20 -- already set to 20 before, but since custome script, use 20 as safer assumption

- threads 2 -- ""

Script: 
```
#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input/output directory
novo_mpileup=${project_dir}/novo_mpileup

#Path and variable for script from PoPoolation to create .sync files
sync=/usr/local/popoolation/mpileup2sync.jar

java -ea -Xmx7g -jar ${sync} --input ${novo_mpileup}/${project_name}.mpileup --output ${novo_mpileup}/${project_name}.sync --fastq-type sanger --min-qual 20 --threads 2
```

## Running Custom R script: Logistic Regression

### Split .sync file into chromosomes (easier to work with)
```
#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
novo_mpileup=${project_dir}/novo_mpileup

grep -v 'Het' ${novo_mpileup}/${project_name}.sync > ${novo_mpileup}/${project_name}_less_het.sync

wait

grep -v 'U' ${novo_mpileup}/${project_name}_less_het.sync > ${novo_mpileup}/${project_name}_removed_U_Het.sync

wait

grep -v 'dmel_mitochondrion_genome' ${novo_mpileup}/${project_name}_removed_U_Het.sync > ${novo_mpileup}/${project_name}_main.sync

wait

rm -f ${novo_mpileup}/${project_name}_less_het.sync

rm -f ${novo_mpileup}/${project_name}_removed_U_Het.sync

grep '3R' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_3R.sync &
grep '2R' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_2R.sync &
grep '3L' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_3L.sync &
grep '2L' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_2L.sync &
grep '^4' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_4.sync &
grep 'X' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_X.sync 
```
Sitting with a different .sync file for each chromosome, can now split more, but since it all goes back together (and is annoying to do all the steps etc.) will make it all into one LONG script... maybe set up as a function?

### LONG SCRIPT: (Change input and output, add Rscripts, etc.)

Some assumptions are made in R scripts that may not work with this data: suggestion is to attempt a trail run with R scripts using 4th chrmosome (the smallest) to make sure all scripts work

```
#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup

#Output dir:
mkdir ${project_dir}/ChromoSubsets
subsets=${project_dir}/ChromoSubsets

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts


# The seperated .sync files
sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 

#Removing 4 b/c complete and short Run on its own (Does not need to be split)

for file in ${sync[@]}
do
name=${file}
base=`basename ${name} .sync`

mkdir ${subsets}/${base}_dir
basedir=${subsets}/${base}_dir

length=($(wc -l ${SyncFiles}/${base}.sync))
echo ${length}

#Split length into 11 segements (11th == length)
cut=$((${length}/11))
cut_2=$((${cut}*2))
cut_3=$((${cut}*3))
cut_4=$((${cut}*4))
cut_5=$((${cut}*5))
cut_6=$((${cut}*6))
cut_7=$((${cut}*7))
cut_8=$((${cut}*8))
cut_9=$((${cut}*9))
cut_10=$((${cut}*10))

###

sed -n " 1, ${cut} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_1.sync
sed -n " $((${cut} + 1)), ${cut_2} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_2.sync
sed -n " $((${cut_2} + 1)), ${cut_3} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_3.sync
sed -n " $((${cut_3} + 1)), ${cut_4} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_4.sync
sed -n " $((${cut_4} + 1)), ${cut_5} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_5.sync
sed -n " $((${cut_5} + 1)), ${cut_6} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_6.sync
sed -n " $((${cut_6} + 1)), ${cut_7} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_7.sync
sed -n " $((${cut_7} + 1)), ${cut_8} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_8.sync
sed -n " $((${cut_8} + 1)), ${cut_9} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_9.sync
sed -n " $((${cut_9} + 1)), ${cut_10} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_10.sync
sed -n " $((${cut_10} + 1)), ${length} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_11.sync
done


mkdir ${subsets}/novo_episodic_4_dir
cp ${SyncFiles}/novo_episodic_4.sync ${subsets}/novo_episodic_4_dir

echo 'Done Splitting Files'

# Should now have 11 different .sync files to work with


###################

# Running Rscript below to go through each subset file created and convert to counts (see Sync_to_counts.R)

echo 'Running R script Sync To Counts'

Rscript ${Rscripts}/Sync_to_counts.R ${subsets}

echo ' Done R script Sync to Counts'

# can remove the sync files? (need a loop to enter each ${subsets} dir to remove *.sync

sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 
sync[5]=${SyncFiles}/novo_episodic_4.sync 

for file in ${sync[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	basedir=${subsets}/${base}_dir
	rm -f ${basedir}/*.sync
done

echo 'Removed Sync Files'

########################

# Run Model On each seperate file (creating a coeffs file)

echo 'Rscript Model'

Rscript ${Rscripts}/Counts_to_model.R ${subsets}

echo 'Done Model'

# Left with the ${subsets} directory full of .csv files and .coeffs.csv files
# Note: this step takes the longest time, possibly search for method to parallelize this script (per chromosome, per file etc..)

#########################

#Combine Coeffs file into one large chromosome file

#mkdir in main project location for combined coeffs 

# Takes into assumption the files are named "episodic_data_2L_11.sync.csv.coeffs.csv" to create a final .csv named "episodic_data_2L_chromo.csv" (removes all past the last _)

mkdir ${project_dir}/novo_coeffs
coeff_dir=${project_dir}/novo_coeffs

Rscript ${Rscripts}/Combine_chromo.R ${subsets} ${coeff_dir}

########################

#Remove all .csv intermediates from random files:

#for file in ${sync[@]}
#	do
#	name=${file}
#	base=`basename ${name} .sync`
#	basedir=${subsets}/${base}_dir
#	rm -f ${basedir}/*.csv
#	rmdir ${base}_dir
#done

#rmdir ${subsets}

#left with the final output files (coeffs) for novoalign files

echo 'DONE'
```


# WORKING SO FAR

Create R-script: need to call variable for output directory above (i.e. the outputs OR subsets above).


### Script: Sync_to_counts.R
To run on own:
```
Rscript Sync_to_counts.R '/home/paul/episodicData/novoalign/novo_mpileup'
```
```
# For loop sync to counts

## need next line to call arguments:

args <- commandArgs(trailingOnly = TRUE)


## Convert a .sync file into long format, filter somewhat, and have only position, treatment, Cage, Generation and Maj/Min counts

## Packages source code: only need these two for this script (need to be this order)

require('tidyr')
require('dplyr')

#1) Need to change details as needed above and below string of #####

#2) Needs a .sync file made by popoolation2

#3) Need to change most importantly for analysis the read in and read out names 

# Read in Data: Big Data Sets

#pwd a direcotry containing only the directories of interest (made with other sed -n script)


mydirs <- list.dirs(path = args[1], recursive = FALSE)

#includes that actual dir.. not with recursive = FALSE

for (dir in mydirs){

    setwd(dir)
  
  mysyncs <- list.files(pattern=".sync")
  
  for (sync in mysyncs){
  
      episodic_counts <- read.table(sync)
    
    #adjust colnames
    print("data read in")
    name.Columns <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR2_115", "SelR1_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "SelR1_0")
    colnames(episodic_counts) <- name.Columns
    
    #Add "replicates" of ancestor -- all are equal
    episodic_counts$SelR2_0 <- episodic_counts$SelR1_0
    episodic_counts$ConR1_0 <- episodic_counts$SelR1_0
    episodic_counts$ConR2_0 <- episodic_counts$SelR1_0
    
    
    #Need the ancestor to stay (after making long) to call major/minor alleles later
    episodic_counts$Ancestor <- episodic_counts$SelR1_0
    
    # Make long by bring populations down
    print("making long")
    long_episodic <- gather(episodic_counts, Population, Allele_Freq , ConR1_115:ConR2_0, factor_key=TRUE)
    
    rm(episodic_counts)
    print("removed counts")
    
    #Error???
    
    # All geneneric below for sync files (only real issue through file is population naming convention)
    ###################################################
    
    #Seperate the allele counts into independent columns for each base
    print("splitting allele freq")
    Episodic_split_2 <- long_episodic %>% 
      separate(Allele_Freq, c("A","T","C","G","N","del"), ":")
    
    rm(long_episodic)
    
    print("removed long")
    
    #Seperate the ancestor to base later things on
    Episodic_split_2 <- Episodic_split_2 %>% 
      separate(Ancestor, c("A_0","T_0","C_0","G_0","N_0","del_0"), ":")
    
    # as.numeric to multiple columns:
    cols.num <- c("A_0", "T_0", "C_0", "G_0", "N_0", "del_0", "A", "T", "C", "G", "N", "del")
    
    #Seems to take a long time for this step?
    Episodic_split_2[cols.num] <- sapply(Episodic_split_2[cols.num],as.numeric)
    
    #Get the sum of all the rows (all the different bases) for each population position:
    
    print("getting row sums")
    Episodic_split_2$sum <- (rowSums(Episodic_split_2[,11:16]))
    
    #Ancestor Major_Allele and minor allele:
    
    # Major allele of ancestor == the maximum positional count
    Episodic_split_2$anc_max <- apply(Episodic_split_2[,4:9], 1, max)
    # Minor is the ancestor second highest count
    Episodic_split_2$anc_min <- apply(Episodic_split_2[,4:9], 1, 
                                      function(x)max(x[x!=max(x)]))
    
    #Major / Minor Base name: match the number of anc_max with the column to call the correct base:
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MajorAllele = ifelse(anc_max== Episodic_split_2[,4], "A", ifelse(anc_max== Episodic_split_2[,5],  "T", ifelse(anc_max== Episodic_split_2[,6],  "C",ifelse(anc_max== Episodic_split_2[,7],  "G", ifelse(anc_max== Episodic_split_2[,8],  "N", ifelse(anc_max== Episodic_split_2[,9],  "del", "N/A" ))))))})
    
    #Major Allele Count of evolved populations; match the Major allele with the count of certain columns for each population 
    
    Episodic_split_2 <- within(Episodic_split_2, {
      Maj_count = ifelse (MajorAllele == "A", Episodic_split_2[,11], ifelse (MajorAllele == "T", Episodic_split_2[,12], ifelse (MajorAllele == "C", Episodic_split_2[,13], ifelse (MajorAllele == "G", Episodic_split_2[,14], ifelse (MajorAllele == "N", Episodic_split_2[,15], ifelse (MajorAllele == "del", Episodic_split_2[,16], "N/A"))))))})
    
    
    # Same thing for minor allele: first ensure that if the sum of all counts == the Major coutn and the ancestor had no minor allele, their is no minor allele (N/A), then follow the same match of anc_min to a certain base
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MinorAllele = ifelse(Maj_count==Episodic_split_2[,17] & anc_min==0, "N/A", ifelse(anc_min== Episodic_split_2[,4], "A", ifelse(anc_min== Episodic_split_2[,5],  "T", ifelse(anc_min== Episodic_split_2[,6],  "C",ifelse(anc_min== Episodic_split_2[,7],  "G", ifelse(anc_min== Episodic_split_2[,8],  "N", ifelse(anc_min== Episodic_split_2[,9],  "del", "Z") ))))))})
    
    
    #Minor Allele Count of the ancestreal minor allele count
    Episodic_split_2 <- within(Episodic_split_2, {
      Min_count = ifelse (MinorAllele == "A", Episodic_split_2[,11], ifelse (MinorAllele == "T", Episodic_split_2[,12], ifelse (MinorAllele == "C", Episodic_split_2[,13], ifelse (MinorAllele == "G", Episodic_split_2[,14], ifelse (MinorAllele == "N", Episodic_split_2[,15],ifelse (MinorAllele == "del", Episodic_split_2[,16],"N/A"))))))})
    
    print("called major and minor alleles and counts")
    # To determine the minor allele base if not specified by the ancestor (new allele brough up etc.)
    
    #max for the population (could be the minor allele)
    Episodic_split_2$maj_all <- apply(Episodic_split_2[,11:16], 1, max)
    
    #alt== second highest count for populations
    Episodic_split_2$alt_allele <- apply(Episodic_split_2[,11:16], 1, 
                                         function(x)max(x[x!=max(x)]))
    
    print("define unknown alleles")
    Episodic_split_2 <- within(Episodic_split_2, {
      Min_count_2 = ifelse (Maj_count == sum, 0, ifelse(Maj_count==maj_all, alt_allele, maj_all))})
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MinorAllele_base = ifelse(Min_count_2==0, "N/A", ifelse(Min_count_2== Episodic_split_2[,11], "A", ifelse(Min_count_2== Episodic_split_2[,12],  "T", ifelse(Min_count_2== Episodic_split_2[,13],  "C",ifelse(Min_count_2== Episodic_split_2[,14],  "G", ifelse(Min_count_2== Episodic_split_2[,15],  "N", ifelse(Min_count_2== Episodic_split_2[,16],  "del", "Z") ))))))})
    
    # Remove unneeded columns (6,7,8,9,10,11,13,14,15)
    Episodic_split_2 <- subset(Episodic_split_2, select = -c(A_0,T_0,C_0,G_0,N_0,del_0,A,T,C,G,N,del,anc_max,anc_min, MinorAllele, Min_count, maj_all, alt_allele))
    
    print("removed unneeded columns")
    nam.col <- c("chr", "pos", "ref", "Population", "sum", "MajorAllele", "Major_count", "Minor_count", "MinorAllele")
    colnames(Episodic_split_2) <- nam.col
    
    
    #Remove unneccessary Columns (as needed)
    #Keep them all for now (except sum) as may be needed later
    #Episodic_split_2 <- subset( Episodic_split_2, select = -ref )
    #Episodic_split_2 <- subset( Episodic_split_2, select = -chr)
    #Episodic_split_2 <- subset( Episodic_split_2, select = -MajorAllele )
    #Episodic_split_2 <- subset( Episodic_split_2, select = -MinorAllele)
    Episodic_split_2<- subset( Episodic_split_2, select = -sum)
    
    
    ## Depends on the filter method:
    print("begin filtering")
    #Filter method: take the sum of each position, and must have at least 5 counts called (i.e over the 16 populations, the total of alleles called for the minor allele must be over 5)
    grp <- Episodic_split_2 %>%
      group_by(pos) %>%
      summarise(sum=sum(Minor_count))
    grp2 <- grp[which(grp$sum<=5),]
    Episodic_split_2 <- Episodic_split_2[!(Episodic_split_2$pos %in% grp2$pos),]
    
    
    #check that the number of obs for episodic_long2 == obs for those without 0's sum (*16 for number of "populations") (needed later as well == grp3)
    
    #grp3 <- grp[-which(grp$sum<=5),]
    rm(grp)
    rm(grp2)
    
    print("remove filter inermediates")
    #################################################
    #Should be all genetic above (from start specificed)
    
    ## Below depends on the population name layout etc. made above
    
    
    #Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)
    
    print("seperate population to Treatment, Generation and Cage")
    episodic_long <- Episodic_split_2 %>%
      separate(Population, c("Treatment", "Generation"), "_")
    
    rm(Episodic_split_2)
    
    episodic_long <- episodic_long %>%
      separate(Treatment, c("Treatment", "Cage"), "R")
    
    cols.num <- c("Cage", "Generation", "Major_count", "Minor_count")
    episodic_long[cols.num] <- sapply(episodic_long[cols.num],as.numeric) 
    
    print("Have final episodic long; now write a csv")
  
    #will need to rename .csv files  
    write.csv(episodic_long, file=paste(sync, ".csv", sep=""))
  
    
    print("wrote csv and now done this .sync file")
  }
}

```


### Script: Counts_to_model.R
To run on own:
```
Rscript Counts_to_model.R 'DIRECTORY'
```
```
#Episodic data analysis: loop .csv files to run model:

args <- commandArgs(trailingOnly = TRUE)

#change to directory holding all directories:

mydirs <- list.dirs(path = args[1], recursive = FALSE)

for (dir in mydirs){

  setwd(dir)
  
  mycsvs <- list.files(pattern=".csv")
  
  for (file in mycsvs){
    
    episodic_long <- read.csv(file, h=T)
    
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
      print(paste("Running entity:", i, "which is", which(position==i), "out of", no.pos, "file=", file))
      
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
      x$chr <- tmp2$chr[1]
      #Add to data frame (total for the whole data set == coeffs_df + the newly made X)
      
      coeffs_df <- rbind(coeffs_df, x)
      
      #Remove i for safety and it starts over
      
      rm(i)
    }
    
    #Change column names to workable
    
    colnames(coeffs_df) <- c("Estimate", "Standard_error", "z-value", "p-value", "position", "chr")
    
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
    
    write.csv(coeffs_df, file=paste(file,".coeffs.csv", sep=""))
        rm(coeffs_df)
    rm(tmp2)
    rm(x)
    rm(modlist_2)
    rm(episodic_long)
    rm(no.pos)
    rm(position)
  }
}

```

### Script: Combine_chromo.R:
To run on own:
```
Rscript Combine_chromo.R 'DIRECTORY' 'OutputDIRECTORY'
```
```

# Combine .coeffs.csv for two mappers
require(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Change to directory holding all directories:

mydirs <- list.dirs(path = args[1], recursive = FALSE)

for (dir in mydirs){

  setwd(dir)
  
  print("Read coeffs.csv files")

  mycsvs <- list.files(pattern=".coeffs.csv")

  Novoalign_Chromosome <- NULL
  
  for (file in mycsvs){
    print(file)
    coeffs1 <- read.csv(file, h=T)
    Novoalign_Chromosome  <- rbind(Novoalign_Chromosome , coeffs1)
    rm(coeffs1)
    J2 <- gsub("*_","", file)
    J3 <- gsub("\\..*","",J2)
}

x3 <- gsub("\\..*","",file)
J3 <- gsub('(.*)_\\w+', '\\1', x3)

X <- args[2]

write.csv(Novoalign_Chromosome , file=paste(X,"/",J3,"_chromo.csv", sep=""), row.names = FALSE)
rm(x3)
rm(J3)
rm(Novoalign_Chromosome)

}
```






______________________________________________

______________________________________________
Links and Notes

Carful with basenames (don't make the outputs novo_aligned_novo_mapped_novo_final.bam etc.
- samtools workflow with GATK:               http://www.htslib.org/workflow/
- CRISP:                     https://github.com/vibansal/crisp
