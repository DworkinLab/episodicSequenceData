# Pipeline; Paul Knoops
## Episodic Data
### Generation 0, 38, 77, 115

### 1. Quality Control and Fastqc

First; check if the data uploaded correctly (matched set) 

Use md5sum on md5.txt. 
	- "-c" = report if checksums match contents of files (OK).
	
	- Must be in rawData directory (with the raw sequence files and the md5.txt file)

```
md5sum - c md5.txt
```

Next check with Fastqc for quality of the reads. 

Fastqc flag "-o" sends all output files to output directory. 

The process will output two files (*fastqc.html and *fastqc.zip). 

The *fastqc.html will be loaded to local machine and opened in web browser to view files.

Moving to local machine also shown below while on the local machine.

```
fastqc -o /home/paul/episodicData/fastqcOutputs /home/paul/episodicData/rawData*.fastq.gz
```

```
scp paul@info.mcmaster.ca:/home/paul/episodicData/fastqcOutputs/*_fastqc.html /Users/paulknoops/episodicWork/data/fastqcOutputs
```

### 2. Quality Trimming files (Trimmomatic)
	FLAGS/INPUT SETTINGS
	- java program
	- trimmomatic v.0.33
	- IlluminaClip = adapter removal (TruSeq3-PE:2:30:10)
	- LEADING & TRAILING = 3; removal at start end end if below quality
	- MINLEN = minimum length (36)
	- MAXINFO = adaptive quality (balance b/w length and quality) = 0.5 (middle ground = safe bet)
	- alternatively, edit input and run PoPoolation2 perl script for trimming
	
	-trimlog <trimlog>

*Trimming F0 and F77

```
#!/bin/bash

cd /usr/local/trimmomatic
dir=/home/paul/episodicData/gen0_gen77_start_to_mapped/Rawdata_2/
files=(${dir}*_R1_001.fastq.gz)
echo ${files[@]}
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`
java -jar trimmomatic-0.33.jar PE -phred33 -trimlog /home/paul/episodicData/gen0_gen77_start_to_mapped/trim2/trimlog_gen0_gen77.txt ${dir}${base}_R1_001.fastq.gz ${dir}${base}_R2_001.fastq.gz /home/paul/episodicData/gen0_gen77_start_to_mapped/trim2/${base}_R1_PE_phred33.fastq.gz /home/paul/episodicData/gen0_gen77_start_to_mapped/trim2/${base}_R1_SE_phred33.fastq.gz /home/paul/episodicData/gen0_gen77_start_to_mapped/trim2/${base}_R2_PE_phred33.fastq.gz /home/paul/episodicData/gen0_gen77_start_to_mapped/trim2/${base}_R2_SE_phred33.fastq.gz ILLUMINACLIP:/usr/local/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done
```
	
### 3. Bring in reference sequence (version r5.57.fasta.gz)
	- curl -O (URL)
	- bwa index ${REFERENCE}

```
# Download index genome for D. mel
# Make outdirectory

index_dir=/home/paul/episodicData/indexSequence
cd ${index_dir}
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

# Set bwa index genome 

bwa index dmel-all-chromosome-r5.57.fasta.gz
```

### 4. BWA Mapping
	- Output = SAM files 
		* Two files (different lanes) for same sample
			-combine now or downstream
			- comine with bwa (or add to current script)
	- use basename:
	- bwa mem
		- t - processors (8)
		- M - Mark shorter split hits as secondary

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
# -M	Mark shorter split hits as secondary (for Picard compatibility).
# Do I need qstat

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE_phred33.fastq.gz`
bwa mem -t 8 -M ${ref_genome} ${dir}/${base}_R1_PE_phred33.fastq.gz ${dir}/${base}_R2_PE_phred33.fastq.gz > ${sam_dir}/${base}_aligned_pe.SAM
done
```

May want to combine lanes

- Can combine L001 with L002 with same start (R1)

- need to change to run those in L003 and L004 and L005 and L006

- find way to add to first script

- example below, but better to samtools merge (shown after SAM to BAM)

```
#!/bin/bash

bwa_dir=/usr/local/bwa/0.7.8
cd ${bwa_dir}
mapped_dir=/home/paul/episodicData/mappedSequence
ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
out_dir=/home/paul/episodicData/mappedSequence
files=(${mapped_dir}/*_L001_aligned_pe.SAM)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.SAM`
bwa mem -t 8 -M ${ref_genome} ${mapped_dir}/${base}_L001_aligned_pe.SAM ${mapped_dir}/${base}_L002_aligned_pe.SAM > ${out_dir}/${base}_aligned_pe.SAM
done
```

### 5. SAM to BAM

- samtools view;  FLAGS
	-Sb = sam to bam
	-q = quality (will do 20 as recommended by PoPoolation2)
	
- pipe to {sametools sort}

```
#! /bin/bash
dir = /home/paul/episodicData/mappedSequence
files = (${dir}/*.SAM)
echo ${files[@]}
for file in ${files[@]}
do
name=${file}
base=`basename ${name}.SAM`
sametools view -b -S -q 20 ${dir}${base}.SAM| sametools sort - ${dir}${base}
done
```

Test with L005

```
#! /bin/bash
dir=/home/paul/episodicData/mappedSequence/
files=(${dir}*_L005_aligned_pe.SAM)
echo ${files[@]}
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L005_aligned_pe.SAM`
samtools view -b -S -q 20 ${dir}${base}_L005_aligned_pe.SAM | samtools sort - ${dir}${base}_L005_aligned_pe
done > L005_SAM_TO_BAM.txt
```

* be sure to change output to signify lane (merge after) -- did above (lane specified)

doing each seperate (based on lane) to make sure get done before server shuts down

Different Outputs for organization: Final Script

```

#! /bin/bash
sam_dir=/home/paul/episodicData/mappedSequence/SAM_files/
bam_dir=/home/paul/episodicData/mappedSequence/BAM_files/
files=(${dir}*.SAM)
echo ${files[@]}
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .SAM`
samtools view -b -S -q 20 ${sam_dir}${base}.SAM | samtools sort - ${bam_dir}${base}
done
```


### 5.5. merge files:

samtools merge --output --input1.bam --input2.bam .... inputN.bam

should be based on lane (differently numbered)

- could change names (all L001 and L002) or set up for even and odd numbers)

```
#!/bin/bash

mapped_dir=/home/paul/episodicData/mappedSequence
files=(${mapped_dir}/*_L001_aligned_pe.SAM)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.SAM`
samtools merge ${mapped_dir}/${base}_merged_aligned_pe.SAM ${mapped_dir}/${base}_L001_aligned_pe.SAM ${mapped_dir}/${base}_L002_aligned_pe.SAM
```

- still need to change for L003/L004 and L005/L006
- tried running: needs to be in BAM format for samtools merge aparently.

FINAL; 

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

	 
### 6.* Sort with Picard?
```
#! /bin/bash

mapped_dir=/home/paul/episodicData/mappedSequence
files=(${mapped_dir}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name}.bam`
java -Xmx2g -jar ˜pic/SortSam.jar I= ${mapped_dir}/*.bam O= ${mapped_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate
done
```

Picard runs with Java

 Xmx2g give Java 2 Gb of memory

jar SortSam use the Java software SortSam

I= input

O= output

SO= sort order; sort by coordinate

VALIDATION STRINGENCY= Picard is like a Princess that is constantly complaining about every small deviation of our sam file from the most stringent requirements. I have never found a sam file satisfying all of Picards demands ⇒ ’shut up Picard’

### 7.* Remove duplicates
```
#! /bin/bash

mapped_dir=/home/paul/episodicData/mappedSequence
files=(${mapped_dir}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name}.sort.bam`
java -Xmx2g -jar ˜pic/MarkDuplicates.jar I= ${mapped_dir}/*.sort.bam O= ${mapped_dir}/${base}.rmd.sort.bam M= ${mapped_dir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```

- I= input file
- O= output file for reads
- M= output file of statistics (how many identified duplicates)
- REMOVE DUPLICATES= remove duplicates from the output file rather than just marking them (remember flag in sam-file 0x400)

### 8.* Remove low quality mapping

	- same as adding -q (step 5.) however, adds flags -f 0x0002 -F 0x0004 -F 0x0008 -b
	- only using -q (20), -F 0x0004, and -b (for bam files)
	- can be done with samtools view above 
	-ex: samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b pe/pe.rmd.sort.bam > pe/pe.qf.rmd.sort.bam
```
samtools view -q 20 -F 0x0004 -b pe/pe.rmd.sort.bam > pe/pe.qf.rmd.sort.bam
```

* 6,7,8 from PoPoolation; PoPoolation2 uses SAM to BAM with sort, than to step 9; including just as a precaution, may not be needed


### 9. Create mpileup file with samtools
	- samtools mpileup (output = .mpileup)
		-FLAGS
			-B = disable BAQ computation
		- BCF and VCF?

```
#! /bin/bash

ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
dir = /home/paul/episodicData/mappedSequence
files = (${dir}*.bam)
echo ${files[@]}
samtools mpileup -B ${ref_genome} ${files[@]} > episodicData.mpileup
```


### 10. Convert to sync file with popoolation2 (java)
	- with java and PoPoolation2 script (mpileup2sync.jar)
		- Flags:
			--input = input files (.mpileup)
			--output = output files (_java.sync)
			--fastq-type = sanger or illumina
			--min-qual = minimum quality score (20)
			--threads = 8
			
Ex.

java -jar ˜/programs/popoolation2/mpileup2sync.jar --input p1-2.mpileup --output p1-2.sync --fastqtype sanger --min-qual 20 --threads 2

```
java -jar ˜/programs/popoolation2/mpileup2sync.jar --input episodicData.mpileup --output episodicData.mpileup.sync --fastqtype illumina --min-qual 20 --threads 8
```
			
### 11. Can now run scripts from popoolation2 directory to find / visualize;
	- calculate allele frequency differences
	- Fst (differentiation b/w pops.) 
		-SCRIPT = fst-sliding
			--input
			--output
			--suppress-noninformative
			--min-count 3
			--min-coverage 10
			--max-coverage 250
			--min-covered-fraction 1
			--window-size = 1for SNP, 500 for sliding window
			--step-size = window size (500)
			--pool-size =
		-SCRIPT to load to IGV = pwc2igv.pl
			--input and --output flags
			-index bam files
		-has option for Fst for genes
		
	- The CMH test (Cochran-Mantel-Haenszel test) - constant allele frequency changes in several replicates
		- SCRIPT = cmh-test.pl
			--input =  (.sync)
			--output = (.cmh)
			--min-count 12 (6)
			--min-coverage 50 (4)
			--max-coverage 200 (120)
			--population 1-2,3-4 (each population tested once)
		- SCRIPT to display on IGV = cmh2gwas.pl
			-set min p-value (--min-pvalue 1.0e-20)
			
	- Fisher's Exact Test (significance of allele frequency differences)
		- SCRIPT = fisher-test.pl
			--input
			--output = (.fet)
			--min-count 6 
			--min-coverage 50 
			--max-coverage 200 
			--suppress-noninformative
		- SCRIPT to IGV = pwc2igv.pl
