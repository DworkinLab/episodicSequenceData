# Full run through of all sequence data - Trimming to mpileup file format
## Need a thourough log
### Change parameters here at top, rest should fall in place
Should run md5sum and fastqc seperatly (before running quality control)
  - set it up so user needs to create project name and create a raw directory (raw_dir) and this will automatically create other directories
  - need to move (mv) all raw files with md5sum files into {project_dir}/raw_dir
  - need known path to project name (i.e /home/paul/episodicData)
  - need to move paths to other directories on machine (i.e bwa or trim) at the top
  - ?? scripts: make a directory and move based on them

Set up/ edit this script
-- Rough locations needed:
```
trimmomatic = /usr/local/trimmomatic
location of trimmomatic on machine
adapt_path = /usr/local/trimmomatic/adapters
path to adapter sequences
?? need to change the apater type!!!
bwa_path = /usr/local/bwa/0.7.8

pic=/usr/local/picard-tools-1.131/picard.jar


merged= ${project_dir}/merged
sort_dir=/home/paul/episodicData/mappedSequence/sort_bam_files
rmd_dir=/home/paul/episodicData/mappedSequence/rmd_bam_files
final_bam=/home/paul/episodicData/mappedSequence/final_bam_files
mpileup_dir=/home/paul/episodicData/mappedSequence

```

###Change so all files like this:  _R1_001.fastq.gz (or R2_001...)




###Create all working Directories and bringing in reference sequence and indexing
Change trim to actual file / input (same for adapters)
```
#! /bin/bash

project_dir = /home/paul/episodicData
raw_dir = ${project_dir}/raw_dir

cd {project_dir}

mkdir ${project_dir}/trim_dir

mkdir ${project_dir}/index_dir
index_dir = ${project_dir}/index_dir

mkdir ${project_dir}/bwa_dir

mkdir ${project_dir}/sam_dir

mkdir ${project_dir}/bam_dir

mkdir ${project_dir}/merged

mkdir ${project_dir}/sort_dir

mkdir ${project_dir}/rmd_dir

mkdir ${project_dir}/final_bam

mkdir ${project_dir}/mpileup_dir

# Can change the index sequence here
cd ${index_dir}
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

bwa index dmel-all-chromosome-r5.57.fasta.gz

ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz

```


Defining all directories (copy to start of all scripts?)
```
project_dir = /home/paul/episodicData
raw_dir = ${project_dir}/raw_dir
trimmomatic = /usr/local/trimmomatic
adapt_path = /usr/local/trimmomatic/adapters
bwa_path = /usr/local/bwa/0.7.8
pic=/usr/local/picard-tools-1.131/picard.jar

trim_dir = ${project_dir}/trim_dir
index_dir = ${project_dir}/index_dir
ref_genome = ${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
bwa_dir = ${project_dir}/bwa_dir
sam_dir = ${project_dir}/sam_dir
bam_dir = ${project_dir}/bam_dir 
merged= ${project_dir}/merged
sort_dir= ${project_dir}/sort_dir
rmd_dir= ${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

```

###Scripts:

### Trimmomatic -- Check the trim log. adapter path
make $trimmomatic/trimmomatic-0.33.jar one imput (no need to change later)

```
#! /bin/bash

files=(${raw_dir}*_R1_001.fastq.gz)
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`
java -jar ${trimmomatic}/trimmomatic-0.33.jar PE -phred33 -trimlog ${trim_dir}/trimlog.txt ${raw_dir}${base}_R1_001.fastq.gz ${raw_dir}${base}_R2_001.fastq.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${out_dir}/${base}_R1_SE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz ${trim_dir}/${base}_R2_SE.fastq.gz ILLUMINACLIP:${adapt_path}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done
```

### BWA mapping
```
#!/bin/bash

cd ${bwa_path}

#make an array for each file in the directory "dir" that ends in _R1_PE_phred33.fastq.gz
files=(${trim_dir}/*_R1_PE.fastq.gz)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE_phred33.fastq.gz`
bwa mem -t 8 -M ${ref_genome} ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${sam_dir}/${base}_aligned_pe.SAM
done
```

### convert SAM to BAM
```
#! /bin/bash


files=(${sam_dir}*.SAM)
echo ${files[@]}
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .SAM`
samtools view -b -S -q 20 ${sam_dir}${base}.SAM | samtools sort - ${bam_dir}${base}
done
```
### Merge files
Only works if all lanes are L001/L002
Is alternative merge method (other, 1st script)

```
#!/bin/bash

files=(${bam_dir}/*_L001_aligned_pe.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.bam`
samtools merge ${merged}/${base}_merged_aligned_pe.bam ${bam_dir}/${base}_L001_aligned_pe.bam ${bam_dir}/${base}_L002_aligned_pe.bam
done
```

### sort with Picard
```
#! /bin/bash

files=(${merged_dir}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
java -Xmx2g -jar ${pic} SortSam I= ${merged_dir}/${base}.bam O= ${sort_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate
done
```

### Remove duplicates
```
#! /bin/bash

files=(${sort_dir}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${sort_dir}/*.sort.bam O= ${rmd_dir}/${base}.rmd.sort.bam M= ${rmd_dir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```

### Remove low quality reads
```
#! /bin/bash

files=(${rmd_dir}/*.rmd.sort.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .rmd.sort.bam`
samtools view -q 20 -F 0x0004 -b ${rmd_dir}/${base}.rmd.sort.bam > ${final_bam}/${base}.final.bam
done
```

### Create mpileup file format
check the flags (illumina or sanger) -6 removed as not what is needed
```
#! /bin/bash

ref_genome=/home/paul/episodicData/indexSequence/dmel-all-chromosome-r5.57.fasta.gz
final_bam=/home/paul/episodicData/mappedSequence/final_bam_files
mpileup_dir=/home/paul/episodicData/mappedSequence
samtools mpileup -B -Q 0 -f ${ref_genome} ${final_bam}/*.bam > ${mpileup_dir}/episodicData_Sanger.mpileup
```
