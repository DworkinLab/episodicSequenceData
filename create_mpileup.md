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
```

###Change so all files like this:  _R1_001.fastq.gz (or R2_001...)




###Create all working Directories and bringing in reference sequence and indexing

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
trim_dir = ${project_dir}/trim_dir
index_dir = ${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
bwa_dir = ${project_dir}/bwa_dir
sam_dir = ${project_dir}/sam_dir
bam_dir = ${project_dir}/bam_dir 
```

###Scripts:

### Trimmomatic -- Check the trim log. adapter path
```
#! /bin/bash

project_dir = /home/paul/episodicData
raw_dir = ${project_dir}/raw_dir
trimmomatic = /usr/local/trimmomatic
adapt_path = /usr/local/trimmomatic/adapters
bwa_path = /usr/local/bwa/0.7.8
trim_dir = ${project_dir}/trim_dir
index_dir = ${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
bwa_dir = ${project_dir}/bwa_dir
sam_dir = ${project_dir}/sam_dir
bam_dir = ${project_dir}/bam_dir 

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

project_dir = /home/paul/episodicData
raw_dir = ${project_dir}/raw_dir
trimmomatic = /usr/local/trimmomatic
adapt_path = /usr/local/trimmomatic/adapters
bwa_path = /usr/local/bwa/0.7.8
trim_dir = ${project_dir}/trim_dir
index_dir = ${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
bwa_dir = ${project_dir}/bwa_dir
sam_dir = ${project_dir}/sam_dir
bam_dir = ${project_dir}/bam_dir 

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


project_dir = /home/paul/episodicData
raw_dir = ${project_dir}/raw_dir
trimmomatic = /usr/local/trimmomatic
adapt_path = /usr/local/trimmomatic/adapters
bwa_path = /usr/local/bwa/0.7.8
trim_dir = ${project_dir}/trim_dir
index_dir = ${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
bwa_dir = ${project_dir}/bwa_dir
sam_dir = ${project_dir}/sam_dir
bam_dir = ${project_dir}/bam_dir 

files=(${sam_dir}*.SAM)
echo ${files[@]}
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .SAM`
samtools view -b -S -q 20 ${sam_dir}${base}.SAM | samtools sort - ${bam_dir}${base}
done
```






