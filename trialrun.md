## Trial

1) moved one set of Con_F77 into new project_dir = Trial

2) moved raw into new dir
 - mkdir raw_dir

3)  make script directory
  - mkdir scripts

4) nano script to make directories
  - change to new project dir
  - # hash out the curl functions (for the test)
 
```
#! /bin/bash
project_dir=/home/paul/Trial
raw_dir=${project_dir}/raw_dir
cd {project_dir}
mkdir ${project_dir}/trim_dir
mkdir ${project_dir}/index_dir
index_dir=${project_dir}/index_dir
mkdir ${project_dir}/bwa_dir
mkdir ${project_dir}/sam_dir
mkdir ${project_dir}/bam_dir
mkdir ${project_dir}/merged
mkdir ${project_dir}/sort_dir
mkdir ${project_dir}/rmd_dir
mkdir ${project_dir}/final_bam
mkdir ${project_dir}/mpileup_dir
# Can change the index sequence here
#cd ${index_dir}
#curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB201$
#bwa index dmel-all-chromosome-r5.57.fasta.gz
#ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
```

All directories here:
```
project_dir=/home/paul/Trial
raw_dir=${project_dir}/raw_dir
trimmomatic=/usr/local/trimmomatic
trim=${trimmomatic}/trimmomatic-0.33.jar
adapt_path=/usr/local/trimmomatic/adapters
adapter=${adapt_path}/TruSeq3-PE.fa:2:30:10
bwa_path=/usr/local/bwa/0.7.8
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir}/trim_dir
bwa_dir=${project_dir}/bwa_dir
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir
```


5) nano script indexing
- cp the hashed out bits
- now run with this
- ** interesting, worked even though ${index_dir} not specified
- worked : but all in home directory..... (take away interesting note)
- used mv dmel* /home/paul/Trial/index_dir to get them to correct place

```
#! /bin/bash

project_dir=/home/paul/Trial
index_dir=${project_dir}/index_dir
cd ${index_dir}
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014$
bwa index dmel-all-chromosome-r5.57.fasta.gz
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
```

6) script for trimming:
- errors: screen said no file; said raw_dir/raw_dir*R1_001 etc.
- change the spacing; be sure if short cut (i.e raw_dir) does not end with / and the line in code has the /
- Was just a missing /
- fixed and running (with a trim log try)
- trim log worked: unknown what it is

```
#! /bin/bash

project_dir=/home/paul/Trial
raw_dir=${project_dir}/raw_dir
trimmomatic=/usr/local/trimmomatic
trim=${trimmomatic}/trimmomatic-0.33.jar
adapt_path=/usr/local/trimmomatic/adapters
adapter=${adapt_path}/TruSeq3-PE.fa:2:30:10
bwa_path=/usr/local/bwa/0.7.8
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir}/trim_dir
bwa_dir=${project_dir}/bwa_dir
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${raw_dir}/*_R1_001.fastq.gz)
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`
java -jar ${trimmomatic}/trimmomatic-0.33.jar PE -phred33 -trimlog ${trim_dir}/trimlog.txt ${raw_dir}/${base}_R1_001.fastq.gz ${raw_dir}/${base}_R2_001.fastq.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R1_SE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz ${trim_dir}/${base}_R2_SE.fastq.gz ILLUMINACLIP:${adapt_path}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done
```

7) BWA mapping
- forgot script screen to make log
- need to change file base name to match trim file names
```
#!/bin/bash

project_dir=/home/paul/Trial
raw_dir=${project_dir}/raw_dir
trimmomatic=/usr/local/trimmomatic
trim=${trimmomatic}/trimmomatic-0.33.jar
adapt_path=/usr/local/trimmomatic/adapters
adapter=${adapt_path}/TruSeq3-PE.fa:2:30:10
bwa_path=/usr/local/bwa/0.7.8
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir}/trim_dir
bwa_dir=${project_dir}/bwa_dir
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

cd ${bwa_path}

#make an array for each file in the directory "dir" that ends in _R1_PE_phred33.fastq.gz
files=(${trim_dir}/*_R1_PE.fastq.gz)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
bwa mem -t 8 -M ${ref_genome} ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${sam_dir}/${base}_aligned_pe.SAM
done
```


