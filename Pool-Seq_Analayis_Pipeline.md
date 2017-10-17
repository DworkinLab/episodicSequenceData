## Analyzing Sequence Data from Pooled Illumina Sequences

The structure of this data set gives an overview of the process to work with the raw outputs from Illumina sequencing of Pooled sequence data (in .fastq.gz format). The Data set I am using is from .......

DNA extraction of populations occurred at MSU by Dr. Michael DeNieu and Mauricio Losilla using the Zymo DNA extractor for insects with population pool sizes of 60 flies. Samples were submitted to RTSF (Research Technology Support Facility) Genomics Core (MSU), where next-generation sequencing libraries were prepared (Illumina TruSeq Nano DNA Library prep kit) and samples were loaded onto two lanes of an Illumina HiSeq 2500 rapid flow cell (version 1). Sequencing was done with Rapid SBS reagents in a 2x150bp paired end format, with the bases called using Illumina Real Time Analysis RTA v1.18.61, and converted to FastQ format. 

A total of 13 sequences were run through Illumina sequencing, the shared ancestor and the 4 treatments/replicates at generation 38, generation 77 and generation 115.

### Step 1: Project Directory

Create a project directory (variable ${project_dir}) for all analysis to be completed within. Many sub directories will be used so proper naming to know each section is vital.

Within this directory, you can create a raw data directory (i.e. raw_dir) for a personal copy of the data or (if the data is to large), links to the locations of the raw data may be more convinient. For explination from here, the variable ${raw_dir} refers to the location of the raw files.


### Quality Control: md5sum

The first step is to check if the data uploaded correctly using md5sum. Using the md5.txt file for the sequence reads will output either "FAILED" or "OK" to know if the raw reads match the line in the md5.txt file, to signify a correct or incorrect transfer.

Flags; 
  
  *-c == report if checksums match contents of files (OK)*
  
```
md5sum - c md5.txt
```

### Quality Control: Fastqc

Using Fastqc from Babraham Bioinformatics (2015), the quality of the reads can be measures. 

Mistakes can occur with sequencing and using Fastqc allows one to view the quality of the reads that the sequencer has found.


Flags;

  *-o == sends all output files to output directory*
  
```
mkdir ${project_dir}/fastqcOutputs
fastqc -o ${project_dir}/fastqcOutputs ${project_dir}/${raw_dir}/*.fastq.gz
```

The process will output two files (fastqc.html and fastqc.zip). 

The fastqc.html can be loaded to local machine and opened in web browser to view files.

```
scp paul@info.mcmaster.ca:${project_dir}/fastqcOutputs/*_fastqc.html ${LOCAL_PROJECT_DIR}/fastqcOutputs
```

### Trimmomatic:
Trimmomatic

Need to make an output directory for trimmomatic outputs: ${trim_dir}

```
mkdir trim_dir
```

Flags;
  
  *-phred33 == may not need to be specified (actually important for later if 33 used)*
  
  *-trimlog == log of trim outputs*
  
  *-LEADING:3 & TRAILING:3 == removal at start end end if below quality*
  
  *-MINLEN:36 == minimum length of 36*
  
  *-MAXINFO:40:0.5 == adaptive quality (balance b/w length and quality)*

  *-IlluminaClip == adapter removal*
        - Need to find the path to the adapter (trimmomatic directory has adapters) and create a variable for the adapter for sequences 
   

```
#! /bin/bash
# Project directory variable:
project_dir=/home/paul/episodicData

# make variable for trimmomatic program location
trim=/usr/local/trimmomatic

# make input directory for raw reads
raw_dir=${project_dir}/raw_dir

# make output directory from trimmomatic outputs
trim_dir=/${project_dir}/trim_dir

# make path to adapter sequences (to be used with ILLUMINACLIP)
adapter=${trim}/adapters/TruSeq3-PE.fa:2:30:10

#list all files to be read (all raw data)
files=(${raw_dir}/*_R1_001.fastq.gz)

#For loop over every file
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`

java -jar ${trim} PE -phred33 \
  -trimlog ${trim_dir}/trimlog.txt \
  ${raw_dir}/${base}_R1_001.fastq.gz \
  ${raw_dir}/${base}_R2_001.fastq.gz \
  ${trim_dir}/${base}_R1_PE.fastq.gz \
  ${trim_dir}/${base}_R1_SE.fastq.gz \
  ${trim_dir}/${base}_R2_PE.fastq.gz \
  ${trim_dir}/${base}_R2_SE.fastq.gz \
  ILLUMINACLIP:${adapter} \
  LEADING:3 \
  TRAILING:3 \
  MAXINFO:40:0.5 \
  MINLEN:36
  
done
```


