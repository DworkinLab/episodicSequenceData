# Note Book for the work I have been doing

### A lot of what is in here is practice work and very rough and messy (other files are clearer for step by step instructions)

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

Change names using mv and copy and paste from below necessary files

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
F38ConR2_TTAGGC_L001_R1_001.fastq.gzx OK
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

Make scripts directory
```
mkdir /home/paul/episodicData/scripts
```
Make script to make all directories and reference seuqence (nano and copy/paste below)
```
#! /bin/bash

project_dir=/home/paul/episodicData
raw_dir=${project_dir}/raw_dir

cd ${project_dir}

mkdir ${project_dir}/trim_dir

mkdir ${project_dir}/index_dir
index_dir=${project_dir}/index_dir

mkdir ${project_dir}/bwa_dir

mkdir ${project_dir}/sam_dir

mkdir ${project_dir}/bam_dir

mkdir ${project_dir}/merged

mkdir ${project_dir}/sort_dir

mkdir ${project_dir}/tmp

mkdir ${project_dir}/rmd_dir

mkdir ${project_dir}/final_bam

mkdir ${project_dir}/mpileup_dir

#Can change the index sequence here
cd ${index_dir}
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

bwa index dmel-all-chromosome-r5.57.fasta.gz

ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
```
named mkdir_and_ref
```
chmod +x mkdir_and_ref
```
made executable
```
screen
```
because reference sequence being brought in could take time
```
mkdir_and_ref
```
runs script

Now need to have all the directories defined to copy and paste into later scripts

Can edit sections here (ending at ref_genome) depending on adapters, reference genome used, project directory etc.
```
project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir
```

Step 6)
Trimmomatic

In scripts
```
nano trim_episodic
```
Copy and paste below

```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${raw_dir}/*_R1_001.fastq.gz)
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`
java -jar ${trim} PE -phred33 -trimlog ${trim_dir}/trimlog.txt ${raw_dir}/${base}_R1_001.fastq.gz ${raw_dir}/${base}_R2_001.fastq.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R1_SE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz ${trim_dir}/${base}_R2_SE.fastq.gz ILLUMINACLIP:${adapter} LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done
```
Change permissions
```
chmod +x trim_episodic
```

Run in screen

```
screen
```

Create a log of all outputs of trim screen
```
script trimscreen.log
```

Run script 
```
trim_episodic
```

Once complete: exit log
```
exit
```

Step 7) Create all scripts and change permissions while trimmomatic runs: Copy and paste each section as applicable
  - BWA mapping

```
nano map_episodic
```

```
#!/bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

cd ${bwa_path}
files=(${trim_dir}/*_R1_PE.fastq.gz)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
bwa mem -t 8 -M ${ref_genome} ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${sam_dir}/${base}_aligned_pe.SAM
done
```

```
chmod +x map_episodic
```

  - Sam to Bam
```
nano samTObam_episodic
```
```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${sam_dir}/*.SAM)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .SAM`
samtools view -b -S -q 20 ${sam_dir}/${base}.SAM | samtools sort - ${bam_dir}/${base}
done
```
```
chmod +x samTObam_episodic
```

  - Merge
```
nano merge_episodic
```

```
#!/bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${bam_dir}/*_L001_aligned_pe.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.bam`
samtools merge ${merged}/${base}_merged_aligned_pe.bam ${bam_dir}/${base}_L001_aligned_pe.bam ${bam_dir}/${base}_L002_aligned_pe.bam
done
```

```
chmod +x merge_episodic
```

  - Picard Sort
```
nano picard_episodic
```

```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${merged}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
java -Xmx2g -Djava.io.tmpdir=${tmp} -jar ${pic} SortSam I= ${merged}/${base}.bam O= ${sort_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${tmp}
done
```

```
chmod +x picard_episodic
```

  - Remove Duplicates

```
nano rmd_episodic
```

```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${sort_dir}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${sort_dir}/*.sort.bam O= ${rmd_dir}/${base}.rmd.sort.bam M= ${rmd_dir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```
Issue: base name in input?
```
files=(${sort_dir}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${sort_dir}/${base}.sort.bam O= ${rmd_dir}/${base}.rmd.sort.bam M= ${rmd_dir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```

```
chmod +x rmd_episodic
```

  - Remove low quality reads
```
nano quality_episodic
```

```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${rmd_dir}/*.rmd.sort.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .rmd.sort.bam`
samtools view -q 20 -F 0x0004 -b ${rmd_dir}/${base}.rmd.sort.bam > ${final_bam}/${base}.final.bam
done
```

```
chmod +x quality_episodic
```

  - Create mpileup
```
nano mpileup_episodic
```

```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

samtools mpileup -B -Q 0 -f ${ref_genome} ${final_bam}/*.bam > ${mpileup_dir}/${project_name}.mpileup
```

```
chmod +x mpileup_episodic
```

  - Create sync file
```
nano sync_episodic
```

```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

java -ea -Xmx7g -jar ${sync} --input ${mpileup_dir}/${project_name}.mpileup --output ${mpileup_dir}/${project_name}.sync --fastq-type sanger --min-qual 20 --threads 2
```

```
chmod +x sync_episodic
```

Step 8) Bwa Mapping
Running from scripts directory puts log's into that section.
```
screen -r
script bwaScreen.log
map_episodic
```
once done, exit to stop log
```
exit
```

Step 9) Sam - Bam 
```
screen -r
script samBamScreen.log
samTObam_episodic
```
Once Done
```
exit
```


Step 10) Merging
```
screen -r
script mergeScreen.log
merge_episodic
```

```
exit
```

Step 11) Picard Sort

```
screen -r
script picardScreen.log
picard_episodic
```
```
exit
```
Step 12) Remove Duplicates
```
screen -r
script rmdScreen.log
rmd_episodic
```

```
exit
```

Step 13) Remove low Quality Reads
```
screen -r
script qualityScreen.log
quality_episodic
```

```
exit
```

REALIZATION: The MGD are not merged and have two difference coverages for Gen 0
- Becuase I am lazy, will just use both as a comparision

Step 14) Create mpileup
```
screen -r
script mpileupScreen.log
mpileup_episodic

```

```
exit
```

Step 15) Create Sync

```
screen -r
script syncScreen.log
sync_episodic
```

```
exit
```

Franssen: indels / repeat masker?


Now have a sync file to be used for CMH test and Fst values.

Order? from final_bam
```
1. F115ConR1_TAGCTT_merged_aligned_pe.final.bam
2. F115ConR2_GGCTAC_merged_aligned_pe.final.bam
3. F115SelR1_GTTTCG_merged_aligned_pe.final.bam
4. F115SelR2_GTGGCC_merged_aligned_pe.final.bam
5. F38ConR1_ATCACG_merged_aligned_pe.final.bam
6. F38ConR2_TTAGGC_merged_aligned_pe.final.bam
7. F38SelR1_ACTTGA_merged_aligned_pe.final.bam
8. F38SelR2_GATCAG_merged_aligned_pe.final.bam
9 .F77ConR1_ATGTCA_merged_aligned_pe.final.bam
10. F77ConR2_ATTCCT_merged_aligned_pe.final.bam
11. F77SelR1_TTAGGC_merged_aligned_pe.final.bam
12. F77SelR2_GATCAG_merged_aligned_pe.final.bam
13. MGD2_SO_CAGATC_merged_aligned_pe.final.bam
14. MGD_SO_CAGATC_merged_aligned_pe.final.bam
```

###CMH Test


?? check values and comparisons?
comparisons are not between same replicates (i.e R1 and R2) but rather based on generations.
Comparisons:
1-5,1-9,1-13,1-14
2-6,2-10,2-13,2-14
3-7,3-11,3-13,3-14
4-8,4-12,4-13,4-14*

* Can only have one value a time: need to add (copy) base populations
* Just copy MGD2 (13.) 3 times at end: should now be 13,15,16,17 that are compared to per generation (4 samples at each time point)
Code is this: added below to script to copy and paste
```
cat ${mpileup_dir}/${project_name}.sync | awk 'BEGIN{OFS="\t"}{print $0,$13,$13,$13}' > ${mpileup_dir}/${project_name}_MGD2.sync
```


From Script Directory
```
nano sync_add_base
```
copy into nano
```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

cat ${mpileup_dir}/${project_name}.sync | awk 'BEGIN{OFS="\t"}{print $0,$13,$13,$13}' > ${mpileup_dir}/${project_name}_MGD2.sync
```
change permissions
```
chmod +x sync_add_base
```
Run on screen with a log
```
screen -r
script catMGD2Screen.log
sync_add_base
```
Stop the screen log
```
exit
```

Error = 
```
Script started, file is catMGD2Screen.log
sync_add_base[paul@info114 scripts]$ sync_add_base
awk: ’BEGIN{OFS=\t}{print
awk: ^ invalid char '�' in expression
[paul@info114 scripts]$ exi
bash: exi: command not found
[paul@info114 scripts]$ exit
exit
Script done, file is catMGD2Screen.log
```
cat p1_p2.sync|awk 'BEGIN{OFS="\t"}{print $0,$4,$5}' > cmh/p1_p2_p1_p2.sync (difference is the quotes: changed above....)
- Rerun all above
(FYI this is old version -- cat ${mpileup_dir}/${project_name}.sync | awk ’BEGIN{OFS="\t"}{print $0,$13,$13,$13}’ > ${mpileup_dir}/${project_name}_MGD2.sync)


Now should be able to run the cmh test with comparisions per generation (B = base, B-115 (x4), B-77 (x4), B-38 (X40)

Comparisons Proper (with Sarah discussion)
```
MGD2 X F115ConR1 & MGD2 X F115ConR2
MGD2 X F115SelR1 & MGD2 X F115SelR2
MGD2 X F38ConR1 & MGD2 X F38ConR2
MGD2 X F38SelR1 & MGD2 X F38SelR2
MGD2 X F77ConR1 & MGD2 X F77ConR2
MGD2 X F77SelR1 & MGD2 X F77SelR2
```
Equvilant to:
```
13 X 1 & 15 X 2
13 X 3 & 15 X 4
13 X 5 & 15 X 6
13 X 7 & 15 X 8
13 X 9 & 15 X 10
13 X 11 & 15 X 12
```
Meaning too many MGD2 copied with cat / awk above..
- just need the two per run:

The script:
-- Check flags
```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

# Can change here to other comparisons

population=11-13,12-15
cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl
mkdir ${mpileup_dir}/${population}
pop_dir=${mpileup_dir}/${population}

perl ${cmh_test} --min-count 3 --min-coverage 10 --max-coverage 250 --population ${population} --input ${mpileup_dir}/${project_name}_MGD2.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-20
```
Flags:
  --min-count 3 
  --min-coverage 10 
  --max-coverage 250
: numbers decided with Ian last term!

What is run in terminal (general):

```
nano cmh_test_${population}
```
Copy script from above and change population to aproptiote groupings (and change for the nano above also)
```
chmod +x cmh_test_${population}
```
screen in scripts is best
```
screen -r
script cmh_${population}_screen.log
cmh_test_${population}
```

```
exit
```

What is run in terminal (example):

```
nano cmh_test_1-13,2-15
```
Copy script from above and change population to aproptiote groupings (and change for the nano above also)
```
chmod +x cmh_test_1-13,2-15
```
screen in scripts is best
```
screen -r
script cmh_1-13,2-15_screen.log
cmh_test_1-13,2-15
```

```
exit
```

```
script cmh_3-13,4-15_screen.log
cmh_test_3-13,4-15

script cmh_5-13,6-15_screen.log
cmh_test_5-13,6-15

script cmh_7-13,8-15_screen.log
cmh_test_7-13,8-15

script cmh_9-13,10-15_screen.log
cmh_test_9-13,10-15

script cmh_11-13,12-15_screen.log
cmh_test_11-13,12-15
```


Try to make script to move all to /Users/paulknoops/Sequence_analysis_2016
1st need to bring all GWAS together into one directory (mkdir all_cmh) and cp all .gwas into it
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/all_cmh/* /Users/paulknoops/Sequence_analysis_2016
```
Load IGV (java -Xmx750m -jar igv.jar) from igv direcotry on local machine
Load all files to IGV
```
1. F115ConR1_TAGCTT_merged_aligned_pe.final.bam
2. F115ConR2_GGCTAC_merged_aligned_pe.final.bam
3. F115SelR1_GTTTCG_merged_aligned_pe.final.bam
4. F115SelR2_GTGGCC_merged_aligned_pe.final.bam
5. F38ConR1_ATCACG_merged_aligned_pe.final.bam
6. F38ConR2_TTAGGC_merged_aligned_pe.final.bam
7. F38SelR1_ACTTGA_merged_aligned_pe.final.bam
8. F38SelR2_GATCAG_merged_aligned_pe.final.bam
9 .F77ConR1_ATGTCA_merged_aligned_pe.final.bam
10. F77ConR2_ATTCCT_merged_aligned_pe.final.bam
11. F77SelR1_TTAGGC_merged_aligned_pe.final.bam
12. F77SelR2_GATCAG_merged_aligned_pe.final.bam
13. MGD2_SO_CAGATC_merged_aligned_pe.final.bam
14. MGD_SO_CAGATC_merged_aligned_pe.final.bam
15. MGD2_SO_CAGATC_merged_aligned_pe.final.bam
16. MGD2_SO_CAGATC_merged_aligned_pe.final.bam
17. MGD2_SO_CAGATC_merged_aligned_pe.final.bam
```



###FST values?

Check pool size?
etc.

```
nano fst_test
```
```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
mpileup_dir=${project_dir}/mpileup_dir
fst_test=/usr/local/popoolation/fst-sliding.pl

perl ${fst_test} --window-size 500 --step-size 500 --suppress-noninformative --input ${mpileup_dir}/${project_name}.sync --min-covered-fraction 1.0 --min-coverage 10 --max-coverage 250 --min-count 3 --output ${mpileup_dir}/${project_name}.fst.txt --pool-size 60
```
```
chmod +x fst_test
```
```
screen -r
script fstScreen.log
fst_test
```
```
exit
```

To view in IGV

```
nano fst_igv
```
```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
mpileup_dir=${project_dir}/mpileup_dir
fst_igv=/usr/local/popoolation/export/pwc2igv.pl

perl ${fst_igv} --input ${mpileup_dir}/${project_name}.fst.txt --output ${mpileup_dir}/${project_name}.fst.igv 
```

```
chmod +x fst_igv
```
```
screen -r
script fst_igv_screen.log
fst_igv
```

```
exit
```

_________________________________________________________________________________________________________
### Re-run with different mapping software
```
cd /home/paul/episodicData
mkdir bowtie
```
remake necessary files in bowtie location (new project dir)
``` 
#! /bin/bash

project_dir=/home/paul/episodicData/bowtie
mkdir ${project_dir}/sam_dir
mkdir ${project_dir}/bam_dir
mkdir ${project_dir}/merged
mkdir ${project_dir}/sort_dir
mkdir ${project_dir}/tmp
mkdir ${project_dir}/rmd_dir
mkdir ${project_dir}/final_bam
mkdir ${project_dir}/mpileup_dir
```

All necessary directories from here
```
project_name=episodic_data_bowtie

#project dir 1 for trim outputs and referenece genome
project_dir1=/home/paul/episodicData

project_dir=/home/paul/episodicData/bowtie

pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir1}/trim_dir

bowtie2_dir=/usr/local/bowtie2/2.2.2

sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir
```


Run bowtie 2

need to run twice; once for paired and once for single (than merge)....https://www.biostars.org/p/165333/
-- Skipped SE on last run through (skip here too? keep consistant.....)
 
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

#Check if correct bowtie2 name?
# -x indicates "base name" of index..... might need to move index?
  # try with just ${ref_genome} -1 ${trim_dir}/${base}_ -2 ${trim_dir}/${base}_ -S ${sam_dir}/${base}_bowtie.sam

files=(${trim_dir}/*_R1_PE.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
${bowtie2_dir}/bowtie2 -x ${ref_genome} -1 ${trim_dir}/${base}_R1_PE.fastq.gz -2 ${trim_dir}/${base}_R2_PE.fastq.gz -S ${sam_dir}/${base}_bowtie_pe.sam
done
```

Check everything makes sense:

###Trial with Bowtie2
Trial directory make all dir's
```
#! /bin/bash

project_dir=/home/paul/episodicData/1misc/Trial
cd ${project_dir}
mkdir ${project_dir}/trim_dir

mkdir ${project_dir}/bwa_dir
mkdir ${project_dir}/sam_dir
mkdir ${project_dir}/bam_dir
mkdir ${project_dir}/merged
mkdir ${project_dir}/sort_dir
mkdir ${project_dir}/tmp
mkdir ${project_dir}/rmd_dir
mkdir ${project_dir}/final_bam
mkdir ${project_dir}/mpileup_dir
```

```
cp /home/paul/episodicData/trim_dir/F115ConR1_TAGCTT_L001_* /home/paul/episodicData/1misc/Trial/trim_dir
```

Run Bowtie
```
nano bowtie_trial
```
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/1misc/Trial
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${trim_dir}/*_R1_PE.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
${bowtie2_dir}/bowtie2 -x ${ref_genome} -1 ${trim_dir}/${base}_R1_PE.fastq.gz -2 ${trim_dir}/${base}_R2_PE.fastq.gz -S ${sam_dir}/${base}_bowtie_pe.sam
done
```

```
chmod +x bowtie_trial
```

```
screen -r
script bowtie_trial.log
bowtie_trial
```

```
exit
```
 Did not work: 
 ```
Could not locate a Bowtie index corresponding to basename "/home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57.fasta.gz"
Error: Encountered internal Bowtie 2 exception (#1)
Command: /usr/local/bowtie2/2.2.2/bowtie2-align-s --wrapper basic-0 -x /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57.fasta.gz -S /home/paul/episodicData/1misc/Trial/sam_dir/F115ConR1_TAGCTT_L001_bowtie_pe.sam -1 /tmp/20100.inpipe1 -2 /tmp/20100.inpipe2
(ERR): bowtie2-align exited with value 1
 ```
 remove the -x (associated with basename of index seqeunce) and test again
 
 ...  
 
 Same error
 Put the -x back in...
 Need to use bowtie2-build for the reference?
 ...
 bowtie2-build hg19.fa hg19
...

make a script for it
- might need it for real trial; so script will be in general episodic direcotry

```
nano bowtie_build_trial
```
```
#! /bin/bash

project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/1misc/Trial

index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta.gz
trim_dir=${project_dir}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 

#using hg19 (just common throughout examples online)

${bowtie2_dir}/bowtie2-build ${ref_genome} hg19
```
Error
```
Settings:
  Output files: "hg19.*.bt2"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 4 (one in 16)
  FTable chars: 10
  Strings: unpacked
  Max bucket size: default
  Max bucket size, sqrt multiplier: default
  Max bucket size, len divisor: 4
  Difference-cover sample period: 1024
  Endianness: little
  Actual local endianness: little
  Sanity checking: disabled
  Assertions: disabled
  Random seed: 0
  Sizeofs: void*:8, int:4, long:8, size_t:8
Input files DNA, FASTA:
  /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57.fasta.gz
Building a SMALL index
Reading reference sizes
Warning: Encountered reference sequence with only gaps
Warning: Encountered reference sequence with only gaps
Warning: Encountered reference sequence with only gaps
Warning: Encountered empty reference sequence
Warning: Encountered reference sequence with only gaps
.... repeated....
Time reading reference sizes: 00:00:01
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
Reference file does not seem to be a FASTA file
  Time to join reference sequences: 00:00:00
Total time for call to driver() for forward index: 00:00:01
Error: Encountered internal Bowtie 2 exception (#1)
Command: bowtie2-build --wrapper basic-0 /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57.fasta.gz hg19 
Deleting "hg19.3.bt2" file written during aborted indexing attempt.
Deleting "hg19.4.bt2" file written during aborted indexing attempt.
Deleting "hg19.1.bt2" file written during aborted indexing attempt.
Deleting "hg19.2.bt2" file written during aborted indexing attempt.
```
b/c the Fasta files are zipped????

```
# copy to retain origional
cp /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57.fasta.gz /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta.gz
gunzip -c /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta.gz
```
no change in files: still .gz
try without -c
```
gunzip /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta.gz
```
worked: new index = no .gz at end

nano to bowtie_build_trial and edit to new index seqeunece to:
```
#! /bin/bash

project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/1misc/Trial

index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta
trim_dir=${project_dir}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 

#Don't use hg19 (examples were just for Humans, no real context

${bowtie2_dir}/bowtie2-build ${ref_genome} dmel-all-chromosome-r5.57_2
```

Try it out - looks to be working!!!

How does the script for bowtie work now with basename index scripts
---> nano into bowtie_trial
ref_genome_base needs to be defined now!
  = /home/paul/episodicData/bowtie/bowtie_indexes
  
  * may need to move bowtie_indexes to index_dir or other places.... change script as needed....
This is for the trial directory...
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/1misc/Trial
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=/home/paul/episodicData/bowtie/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir


files=(${trim_dir}/*_R1_PE.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
${bowtie2_dir}/bowtie2 -x ${ref_genome_base} -1 ${trim_dir}/${base}_R1_PE.fastq.gz -2 ${trim_dir}/${base}_R2_PE.fastq.gz -S ${sam_dir}/${base}_bowtie_pe.sam
done
```
```
screen -r
script bowtie_trial#2.log
bowtie_trial
```

```
exit
```

Worked: 
Run through for full script:
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir


files=(${trim_dir}/*_R1_PE.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
${bowtie2_dir}/bowtie2 -x ${ref_genome_base} -1 ${trim_dir}/${base}_R1_PE.fastq.gz -2 ${trim_dir}/${base}_R2_PE.fastq.gz -S ${sam_dir}/${base}_bowtie_pe.sam
done
```
```
nano bowtie2_script
```

```
chmod +x bowtie2_script
```
Make sure screen -r is in correct script dir
```
screen -r
script bowtie2.log
bowtie2_script
```
```
exit
```


### Runthrough
Run through rest of scripts the same
Change the project_dir and make sure names will be different in the end...
Try with generic script???

 Sam to Bam
 - need to change SAM to sam
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${sam_dir}/*.sam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sam`
samtools view -b -S -q 20 ${sam_dir}/${base}.sam | samtools sort - ${bam_dir}/${base}
done
```

```
nano bowtie_samToBam
```

```
chmod +x bowtie_samToBam
```

```
screen -r
script bowtie2_SamBam.log
bowtie_samToBam
```
```
exit
```

Merge:

```
#!/bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${bam_dir}/*_L001_bowtie_pe.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_bowtie_pe.bam`
samtools merge ${merged}/${base}_merged_bowtie_pe.bam ${bam_dir}/${base}_L001_bowtie_pe.bam ${bam_dir}/${base}_L002_bowtie_pe.bam
done
```

```
nano bowtie_merge
```
```
chmod +x bowtie_merge
```

```
screen -r
script bowtie_merge.log
bowtie_merge
```
```
exit
```
 Picard Sort
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${merged}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
java -Xmx2g -Djava.io.tmpdir=${tmp} -jar ${pic} SortSam I= ${merged}/${base}.bam O= ${sort_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${tmp}
done
```

```
nano bowtie_picsort
```
```
chmod +x bowtie_picsort
screen -r
script bowtie_picsort.log
bowtie_picsort
```

```
exit
```

Remove Duplicates

```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${sort_dir}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${sort_dir}/${base}.sort.bam O= ${rmd_dir}/${base}.rmd.sort.bam M= ${rmd_dir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```
```
nano bowtie_rmd
```

```
chmod +x bowtie_rmd
screen -r
script bowtie_rmd.log
bowtie_rmd
```

```
exit
```

Remove low quality reads
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

files=(${rmd_dir}/*.rmd.sort.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .rmd.sort.bam`
samtools view -q 20 -F 0x0004 -b ${rmd_dir}/${base}.rmd.sort.bam > ${final_bam}/${base}.final.bam
done
```

```
nano bowtie_quality
```

```
chmod +x bowtie_quality
screen -r
script bowtie_quality.log
bowtie_quality
```

```
exit
```
Layout of final bam dir
```
F115ConR1_TAGCTT_merged_bowtie_pe.final.bam
F115ConR2_GGCTAC_merged_bowtie_pe.final.bam
F115SelR1_GTTTCG_merged_bowtie_pe.final.bam
F115SelR2_GTGGCC_merged_bowtie_pe.final.bam
F38ConR1_ATCACG_merged_bowtie_pe.final.bam
F38ConR2_TTAGGC_merged_bowtie_pe.final.bam
F38SelR1_ACTTGA_merged_bowtie_pe.final.bam
F38SelR2_GATCAG_merged_bowtie_pe.final.bam
F77ConR1_ATGTCA_merged_bowtie_pe.final.bam
F77ConR2_ATTCCT_merged_bowtie_pe.final.bam
F77SelR1_TTAGGC_merged_bowtie_pe.final.bam
F77SelR2_GATCAG_merged_bowtie_pe.final.bam
MGD2_SO_CAGATC_merged_bowtie_pe.final.bam
MGD_SO_CAGATC_merged_bowtie_pe.final.bam
```


mpileup
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

samtools mpileup -B -Q 0 -f ${ref_genome} ${final_bam}/*.bam > ${mpileup_dir}/${project_name}.mpileup
```

```
nano bowtie_mpileup
```

```
chmod +x bowtie_mpileup
screen -r
script bowtie_mpileup.log
bowtie_mpileup
```

```
exit
```

sync
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

java -ea -Xmx7g -jar ${sync} --input ${mpileup_dir}/${project_name}.mpileup --output ${mpileup_dir}/${project_name}.sync --fastq-type sanger --min-qual 20 --threads 2
```

```
nano bowtie_sync
```

```
chmod +x bowtie_sync
screen -r
script bowtie_sync.log
bowtie_sync
```

```
exit
```

Add MGD2's (cp from 1st run through and edit the project_dir
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

cat ${mpileup_dir}/${project_name}.sync | awk 'BEGIN{OFS="\t"}{print $0,$13,$13,$13}' > ${mpileup_dir}/${project_name}_MGD2.sync

```

```
screen -r
script bowtie_addMGD.log
bowtie_sync_add_base
```
```
exit
```

CMH TEST:

```
MGD2 X F115ConR1 & MGD2 X F115ConR2
MGD2 X F115SelR1 & MGD2 X F115SelR2
MGD2 X F38ConR1 & MGD2 X F38ConR2
MGD2 X F38SelR1 & MGD2 X F38SelR2
MGD2 X F77ConR1 & MGD2 X F77ConR2
MGD2 X F77SelR1 & MGD2 X F77SelR2
```
Equvilant to:
```
13 X 1 & 15 X 2
13 X 3 & 15 X 4
13 X 5 & 15 X 6
13 X 7 & 15 X 8
13 X 9 & 15 X 10
13 X 11 & 15 X 12
```


The script:

```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

# Can change here to other comparisons (populations)

population=1-13,2-15
cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl
mkdir ${mpileup_dir}/${population}
pop_dir=${mpileup_dir}/${population}

perl ${cmh_test} --min-count 3 --min-coverage 10 --max-coverage 250 --population ${population} --input ${mpileup_dir}/${project_name}_MGD2.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-20
```

What is run in terminal (general):

```
nano cmh_test_${population}
```
Copy script from above and change population to aproptiote groupings (and change for the nano above also)
```
chmod +x cmh_test_${population}
screen -r
script cmh_${population}_screen.log
cmh_test_${population}
```
```
exit
```
Running all
```
script cmh_bowtie_1-13,2-15_screen.log
cmh_test_bowtie_1-13,2-15


script cmh_bowtie_3-13,4-15_screen.log
cmh_test_bowtie_3-13,4-15

script cmh_bowtie_5-13,6-15_screen.log
cmh_test_bowtie_5-13,6-15


script cmh_bowtie_7-13,8-15_screen.log
cmh_test_bowtie_7-13,8-15

script cmh_bowtie_9-13,10-15_screen.log
cmh_test_bowtie_9-13,10-15

script cmh_bowtie_11-13,12-15_screen.log
cmh_test_bowtie_11-13,12-15

```



Move to local machine: all need to be put into relevant dir (mkdir all_bowtie_cmh) in mpileup dir (cp all .gwas from the dirctories to all_bowtie_cmh

```
scp paul@info.mcmaster.ca:/home/paul/episodicData/bowtie/mpileup_dir/all_bowtie_cmh/* /Users/paulknoops/Sequence_analysis_2016/bowtie2
```

To open IGV...

```
java -Xmx2g -jar /Users/paulknoops/igv/IGV_2.3.67/igv.jar
```


### Check the F77 zip files
* The files from Ian had those that failed zipped, and the ones I zipped worked:
* will try to unzip and rezip all together and test the md5 again


```
mkdir Gen77
mv F77* Gen77
mv MDG_* Gen77
mv md5_F77* Gen77
mv mdg_MGD* Gen77
```

```
screen -r
gunzip *.gz
```

```
gzip *.fastq
```

```
screen -r
script md5sum_F77_MGD.log
md5sum -c md5_F77.txt
```
```
md5sum -c md5_MGD_SO.txt
```
```
exit
```
All failed after re-zipping... 


MGD from Synthetic outbred (most used file for ancestor) all worked!

New md5sum file from Ian;

All working from the head!
```
[paul@infoserv ~]$ cd /scratch/dworkinlab_backups/drosophilaPredation/gen0_gen77_sequence/
[paul@infoserv gen0_gen77_sequence]$ ls
checklist.chk				 Con_R2_F77_ATTCCT_L004_R2_001.fastq.gz*  Sel_R1_F77_TTAGGC_L004_R1_001.fastq*
Con_R1_F77_ATGTCA_L003_R1_001.fastq.gz*  MGD_SO_CAGATC_L005_R1_001.fastq.gz*	  Sel_R1_F77_TTAGGC_L004_R2_001.fastq*
Con_R1_F77_ATGTCA_L003_R2_001.fastq.gz*  MGD_SO_CAGATC_L005_R2_001.fastq*	  Sel_R2_F77_GATCAG_L003_R1_001.fastq.gz*
Con_R1_F77_ATGTCA_L004_R1_001.fastq.gz*  MGD_SO_CAGATC_L006_R1_001.fastq*	  Sel_R2_F77_GATCAG_L003_R2_001.fastq.gz*
Con_R1_F77_ATGTCA_L004_R2_001.fastq.gz*  MGD_SO_CAGATC_L006_R2_001.fastq*	  Sel_R2_F77_GATCAG_L004_R1_001.fastq.gz*
Con_R2_F77_ATTCCT_L003_R1_001.fastq.gz*  readme.txt				  Sel_R2_F77_GATCAG_L004_R2_001.fastq.gz*
Con_R2_F77_ATTCCT_L003_R2_001.fastq.gz*  Sel_R1_F77_TTAGGC_L003_R1_001.fastq*
Con_R2_F77_ATTCCT_L004_R1_001.fastq.gz*  Sel_R1_F77_TTAGGC_L003_R2_001.fastq*
[paul@infoserv gen0_gen77_sequence]$ md5sum -c checklist.chk 
Con_R1_F77_ATGTCA_L003_R1_001.fastq.gz: OK
Con_R1_F77_ATGTCA_L003_R2_001.fastq.gz: OK
Con_R1_F77_ATGTCA_L004_R1_001.fastq.gz: OK
Con_R1_F77_ATGTCA_L004_R2_001.fastq.gz: OK
Con_R2_F77_ATTCCT_L003_R1_001.fastq.gz: OK
Con_R2_F77_ATTCCT_L003_R2_001.fastq.gz: OK
Con_R2_F77_ATTCCT_L004_R1_001.fastq.gz: OK
Con_R2_F77_ATTCCT_L004_R2_001.fastq.gz: OK
MGD_SO_CAGATC_L005_R1_001.fastq.gz: OK
MGD_SO_CAGATC_L005_R2_001.fastq: OK
MGD_SO_CAGATC_L006_R1_001.fastq: OK
MGD_SO_CAGATC_L006_R2_001.fastq: OK
Sel_R1_F77_TTAGGC_L003_R1_001.fastq: OK
Sel_R1_F77_TTAGGC_L003_R2_001.fastq: OK
Sel_R1_F77_TTAGGC_L004_R1_001.fastq: OK
Sel_R1_F77_TTAGGC_L004_R2_001.fastq: OK
Sel_R2_F77_GATCAG_L003_R1_001.fastq.gz: OK
Sel_R2_F77_GATCAG_L003_R2_001.fastq.gz: OK
Sel_R2_F77_GATCAG_L004_R1_001.fastq.gz: OK
Sel_R2_F77_GATCAG_L004_R2_001.fastq.gz: OK
```


So now have a set working files: no need to worry about them now (still odd with the asterix with them & the files sizes)

###Now:
Additional comparisons???

CMH TEST:

```
F115SelR1 X F115ConR1 & F115SelR2 X F115ConR2
F38SelR1 X F38ConR1 & F38SelR2 X F38ConR2
F77SelR1 X F77ConR1 & F77SelR2 X F77ConR2
```
Equvilant to:
```
1 x 3, 2 x 4
5 x 7, 6 x 8
9 x 11, 10 x 12
```

The script from before for Bowtie

- change population (and change start to match bowtie vs. BWA):

### Possibly need to change coverages etc.
For now match with before than worry about changes

```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

# Can change here to other comparisons (populations)

population=9-11,10-12

cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl
mkdir ${mpileup_dir}/${population}
pop_dir=${mpileup_dir}/${population}

perl ${cmh_test} --min-count 3 --min-coverage 10 --max-coverage 250 --population ${population} --input ${mpileup_dir}/${project_name}_MGD2.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-20
```

General:

```
nano cmh_test_${population}
```
Copy script from above and change population to aproptiote groupings (and change for the nano above also)
```
chmod +x cmh_test_${population}
screen -r
script cmh_${population}_screen.log
cmh_test_${population}
```
```
exit
```
Running all
```
script cmh_bowtie_1-3,2-4_screen.log
cmh_test_bowtie_1-3,2-4

script cmh_bowtie_5-7,6-8_screen.log
cmh_test_bowtie_5-7,6-8

script cmh_bowtie_9-11,10-12_screen.log
cmh_test_bowtie_9-11,10-12


```

###On info113

Repeat for non-bowtie

```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

# Can change here to other comparisons

population=1-3,2-4
cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl
mkdir ${mpileup_dir}/${population}
pop_dir=${mpileup_dir}/${population}

perl ${cmh_test} --min-count 3 --min-coverage 10 --max-coverage 250 --population ${population} --input ${mpileup_dir}/${project_name}_MGD2.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-20

```

```
script cmh_1-3,2-4_screen.log
cmh_test_1-3,2-4

script cmh_5-7,6-8_screen.log
cmh_test_5-7,6-8

script cmh_9-11,10-12_screen.log
cmh_test_9-11,10-12


```
Check logs for: bowtie_1-3,2-4 -- Oddity
- once opened = not complete file!!!!!!

From Screen log
```
Script started on Wed 03 Aug 2016 12:27:26 PM EDT
cmh_test_bowtie_1-3,2-4ESC]0;paul@info113:~/episodicData/bowtie/scriptsESC\[paul@info113 scripts]$ cmh_test_bowtie_1-3,2-4
Reading sync file and writing temporary R output file
Calling R, to calculate the Cochran-Mantel-Haenszel test statistic
Error in mantelhaen.test(array(c(0, 0, 1, 0, 3, 18, 8, 0), dim = c(2,  : 
  sample size in each stratum must be > 1
Execution halted
Parsing R-output and writing output file
Use of uninitialized value $pvalue in scalar chomp at /usr/local/popoolation/cmh-test.pl line 132, <$ifh> line 3777105.
Use of uninitialized value $pvalue in substitution (s///) at /usr/local/popoolation/cmh-test.pl line 137, <$ifh> line 3777105.
Use of uninitialized value $pvalue in numeric gt (>) at /usr/local/popoolation/cmh-test.pl line 138, <$ifh> line 3777105.
Use of uninitialized value $pvalue in concatenation (.) or string at /usr/local/popoolation/cmh-test.pl line 139, <$ifh> line 3777105.
Done
Argument "0:0:73:59:0:0" isn't numeric in numeric lt (<) at /usr/local/popoolation/export/cmh2gwas.pl line 39, <$ifh> line 1888553.

```
is it the p-value/???? -- error was before that....

Re run-  (log = with retry) = another fail..... same error
=- But ran with not these two and ran with other (non-bowtie) comparison
--> ??

pwd each new directory with the file: copy to start of scp (after the "scp paul@info.mcmaster.ca:" and before "/*.gwas /Users/paulknoops/Sequence_analysis_2016"

```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/9-11,10-12/*.gwas /Users/paulknoops/Sequence_analysis_2016
```


To open IGV...

```
java -Xmx2g -jar /Users/paulknoops/igv/IGV_2.3.67/igv.jar
```


Changing the parameters:
-- may want to change details of cmh test (namely min/max coverage and min p-value)

- Move all scripts with 1st parameters to new dir (likely cmh_10_250_20 ; For each level of coverage and pvalue used (cmh_minCov_maxCov_minP-value)

-Make new with changed parameters (cmh_bowtie_5_100_40)?

The script with new parameters: just change population and add something different to names for script (i.e #2 likely for second run)
-- need to change the new mkdir with populations
Bowtie:
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
pic=/usr/local/picard-tools-1.131/picard.jar
sync=/usr/local/popoolation/mpileup2sync.jar
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

# Can change here to other comparisons (populations)

population=5-13,6-15

cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl
mkdir ${mpileup_dir}/${population}_2
pop_dir=${mpileup_dir}/${population}_2

perl ${cmh_test} --min-count 3 --min-coverage 5 --max-coverage 100 --population ${population} --input ${mpileup_dir}/${project_name}_MGD2.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-40
```

```
nano cmh_test2_bowtie_11-13,12-15

nano cmh_test2_bowtie_5-7,6-8

nano cmh_test2_bowtie_1-13,2-15

nano cmh_test2_bowtie_7-13,8-15

nano cmh_test2_bowtie_1-3,2-4

nano cmh_test2_bowtie_9-11,10-12

nano cmh_test2_bowtie_3-13,4-15

nano cmh_test2_bowtie_9-13,10-15

nano cmh_test2_bowtie_5-13,6-15

```
In the dir
```
chmod +x cmh*
```
for each seperatly (on a different screen)
```
screen
script cmh_2_bowtie_${population}_screen.log
cmh_test2_bowtie_${population}
```
```
exit
```

```
screen
script cmh_2_bowtie_11-13,12-15_screen.log
cmh_test2_bowtie_11-13,12-15

screen
script cmh_2_bowtie_1-13,2-15_screen.log
cmh_test2_bowtie_1-13,2-15

screen
script cmh_2_bowtie_1-3,2-4_screen.log
cmh_test2_bowtie_1-3,2-4

screen
script cmh_2_bowtie_3-13,4-15_screen.log
cmh_test2_bowtie_3-13,4-15

screen
script cmh_2_bowtie_5-13,6-15_screen.log
cmh_test2_bowtie_5-13,6-15

screen
script cmh_2_bowtie_5-7,6-8_screen.log
cmh_test2_bowtie_5-7,6-8

screen
script cmh_2_bowtie_7-13,8-15_screen.log
cmh_test2_bowtie_7-13,8-15

screen
script cmh_2_bowtie_9-11,10-12_screen.log
cmh_test2_bowtie_9-11,10-12

screen
script cmh_2_bowtie_9-13,10-15_screen.log
cmh_test2_bowtie_9-13,10-15
```



Non-Bowtie

- could try a loop; runs much longer becuase all in sequence, but easier to make one script?
```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

# Can change here to other comparisons

pop[0]=11-13,12-15
pop[1]=1-13,2-15
pop[2]=1-3,2-4
pop[3]=3-13,4-15
pop[4]=5-13,6-15
pop[5]=5-7,6-8
pop[6]=7-13,8-15
pop[7]=9-11,10-12
pop[8]=9-13,10-15

cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl

for population in ${pop[@]}
do
mkdir ${mpileup_dir}/${population}_2
pop_dir=${mpileup_dir}/${population}_2

perl ${cmh_test} --min-count 3 --min-coverage 5 --max-coverage 100 --population ${population} --input ${mpileup_dir}/${project_name}_MGD2.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-40

done
```

Forget these individuals (for non-for loop if needed when/If this fails:
```
nano cmh_test2_11-13,12-15

nano cmh_test2_1-13,2-15

nano cmh_test2_1-3,2-4

nano cmh_test2_3-13,4-15

nano cmh_test2_5-13,6-15

nano cmh_test2_5-7,6-8

nano cmh_test2_7-13,8-15

nano cmh_test2_9-11,10-12

nano cmh_test2_9-13,10-15
```


Script for the for loop:
in the appropriote dir

```
nano cmh_test2_FORLOOPS
```
```
chmod +x cmh_test2_FORLOOPS
```
```
screen
script cmh_test2_FORLOOPS.log
cmh_test2_FORLOOPS
```
dettach and make a note on which screen this is* 

```
59579.pts-5.info113
```
-- seems to be working

```
exit
```

Noticed: in script: make sure output is different///// will need to mv to make sure on that (change names)
Changing names once all in same dir
for Bowtie:
```
mv episodic_data_bowtie_11-13,12-15.cmh.gwas episodic_data_bowtie_2_11-13,12-15.cmh.gwas
mv episodic_data_bowtie_5-7,6-8.cmh.gwas episodic_data_bowtie_2_5-7,6-8.cmh.gwas
mv episodic_data_bowtie_1-13,2-15.cmh.gwas episodic_data_bowtie_2_1-13,2-15.cmh.gwas
mv episodic_data_bowtie_7-13,8-15.cmh.gwas episodic_data_bowtie_2_7-13,8-15.cmh.gwas
mv episodic_data_bowtie_1-3,2-4.cmh.gwas episodic_data_bowtie_2_1-3,2-4.cmh.gwas
mv episodic_data_bowtie_9-11,10-12.cmh.gwas episodic_data_bowtie_2_9-11,10-12.cmh.gwas
mv episodic_data_bowtie_3-13,4-15.cmh.gwas episodic_data_bowtie_2_3-13,4-15.cmh.gwas
mv episodic_data_bowtie_9-13,10-15.cmh.gwas episodic_data_bowtie_2_9-13,10-15.cmh.gwas
mv episodic_data_bowtie_5-13,6-15.cmh.gwas episodic_data_bowtie_2_5-13,6-15.cmh.gwas
```
For BWA:
```
mv episodic_data_11-13,12-15.cmh.gwas episodic_data_2_11-13,12-15.cmh.gwas
mv episodic_data_5-7,6-8.cmh.gwas episodic_data_2_5-7,6-8.cmh.gwas
mv episodic_data_1-13,2-15.cmh.gwas episodic_data_2_1-13,2-15.cmh.gwas
mv episodic_data_7-13,8-15.cmh.gwas episodic_data_2_7-13,8-15.cmh.gwas
mv episodic_data_1-3,2-4.cmh.gwas episodic_data_2_1-3,2-4.cmh.gwas
mv episodic_data_9-11,10-12.cmh.gwas episodic_data_2_9-11,10-12.cmh.gwas
mv episodic_data_3-13,4-15.cmh.gwas episodic_data_2_3-13,4-15.cmh.gwas
mv episodic_data_9-13,10-15.cmh.gwas episodic_data_2_9-13,10-15.cmh.gwas
mv episodic_data_5-13,6-15.cmh.gwas episodic_data_2_5-13,6-15.cmh.gwas
```

Redefine the gwas_dir
For Bowtie:
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/bowtie/mpileup_dir/all_bowtie_cmh_2 /Users/paulknoops/Sequence_analysis_2016/bowtie_igv_2nd
```
For BWA:
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/all_cmh_2/*.gwas /Users/paulknoops/Sequence_analysis_2016/bwa_igv_2nd
```


###Getting smaller Sync file to practice:

```
grep '3R' episodic_data_MGD2.sync > episodic_data_3R.sync
```
Has 3RHet:

-v = opposite?
renames files!
```
grep -v 'Het' episodic_data_3R_Het.sync > episodic_data_3R.sync
```
IDEA::
Can remove all with Het and Extra here --> something on the lines of:
```
grep -v 'Het' episodic_data.sync > episodic_data_lessHet.sync
```
::

```
wc -l episodic_data_3R.sync
```
```
27293799 episodic_data_3R.sync
```

Example:
```
sed -n '16224,16482 p' orig-data-file > new-file
```
For me: random selection --
```
sed -n '11111111,11112111 p' episodic_data_3R.sync > episodic_data_Ian_subset.sync
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/episodic_data_Ian_subset.sync /Users/paulknoops/Sequence_analysis_2016
```

#DUMB DUMB Realization: 
when doing Cat and Awk adding in MGD2:
Copied $13 BUT that was actually row 10 (one from F77.....)
  - do the merge and re-run a bunch of stuff
  - For Ian, give from non-Cat/Awk stuff (redo some of above....)

### Practice with .txt file (CMH done)

- to remove unneeded areas:

```
grep -v 'Het' practice_episodic_data_1-3,2-4.cmh.txt > practice_episodic_data_1-3,2-4_2.cmh.txt
grep -v 'U' practice_episodic_data_1-3,2-4_2.cmh.txt > practice_episodic_data_1-3,2-4_3.cmh.txt
```


### Merge MGD and MGD2
- need to do with final bam files

1st = bwa -mem
```
MGD2_SO_CAGATC_merged_aligned_pe.final.bam
MGD_SO_CAGATC_merged_aligned_pe.final.bam
```
from final bam directory:
```
samtools merge MGD3_SO_CAGATC_merged_aligned_pe.final.bam MGD2_SO_CAGATC_merged_aligned_pe.final.bam MGD_SO_CAGATC_merged_aligned_pe.final.bam
```

2nd = bowtie2
```
MGD2_SO_CAGATC_merged_bowtie_pe.final.bam
MGD_SO_CAGATC_merged_bowtie_pe.final.bam
```
from bowtie final bam 
```
samtools merge MGD3_SO_CAGATC_merged_bowtie_pe.final.bam MGD2_SO_CAGATC_merged_bowtie_pe.final.bam MGD_SO_CAGATC_merged_bowtie_pe.final.bam
```

For both (before rerunning the mpileup and such scripts) = remove old non merged bases
  - will change numbers for mpileup......

As long as old is removed, delete old mpileup and sync and rerun the scripts.

Layout for BWA
```
F115ConR1_TAGCTT_merged_aligned_pe.final.bam
F115ConR2_GGCTAC_merged_aligned_pe.final.bam
F115SelR1_GTTTCG_merged_aligned_pe.final.bam
F115SelR2_GTGGCC_merged_aligned_pe.final.bam
F38ConR1_ATCACG_merged_aligned_pe.final.bam
F38ConR2_TTAGGC_merged_aligned_pe.final.bam
F38SelR1_ACTTGA_merged_aligned_pe.final.bam
F38SelR2_GATCAG_merged_aligned_pe.final.bam
F77ConR1_ATGTCA_merged_aligned_pe.final.bam
F77ConR2_ATTCCT_merged_aligned_pe.final.bam
F77SelR1_TTAGGC_merged_aligned_pe.final.bam
F77SelR2_GATCAG_merged_aligned_pe.final.bam
MGD3_SO_CAGATC_merged_aligned_pe.final.bam
```
Run script in screen

```
mpileup_episodic
```

Layout for Bowtie
```
F115ConR1_TAGCTT_merged_bowtie_pe.final.bam
F115ConR2_GGCTAC_merged_bowtie_pe.final.bam
F115SelR1_GTTTCG_merged_bowtie_pe.final.bam
F115SelR2_GTGGCC_merged_bowtie_pe.final.bam
F38ConR1_ATCACG_merged_bowtie_pe.final.bam
F38ConR2_TTAGGC_merged_bowtie_pe.final.bam
F38SelR1_ACTTGA_merged_bowtie_pe.final.bam
F38SelR2_GATCAG_merged_bowtie_pe.final.bam
F77ConR1_ATGTCA_merged_bowtie_pe.final.bam
F77ConR2_ATTCCT_merged_bowtie_pe.final.bam
F77SelR1_TTAGGC_merged_bowtie_pe.final.bam
F77SelR2_GATCAG_merged_bowtie_pe.final.bam
MGD3_SO_CAGATC_merged_bowtie_pe.final.bam
```

Run script 

```
bowtie_mpileup
```

Run sync scripts:
```
sync_episodic
```

```
sync_bowtie
```

Possibly do not need to cat/awk Base generation again, can just call same population in Popoolation CMH test

Remove unneeded areas (with 'Het' and 'U')
BWA:
```
grep -v 'Het' episodic_data.sync > episodic_data_less_het.sync
grep -v 'U' episodic_data_less_het.sync > episodic_data_removed_U_Het.sync
grep -v 'dmel_mitochondrion_genome' episodic_data_removed_U_Het.sync > episodic_data_main.sync
```

For Bowtie
```
grep -v 'Het' episodic_data_bowtie.sync > episodic_data_bowtie_less_het.sync
grep -v 'U' episodic_data_bowtie_less_het.sync > episodic_data_bowtie_removed_U_Het.sync
grep -v 'dmel_mitochondrion_genome' episodic_data_bowtie_removed_U_Het.sync > episodic_data_bowtie_main.sync
```

Remove all old CMH tests that compared with the base generation: copied the incorrect one remember with cat/awk!!!


Using the subsample synchronized script help::::

```
perl subsample_synchronized.pl --help
```
```
SUBSAMPLE-SYNCHRONIZEDUser Contributed Perl DocumentaSUBSAMPLE-SYNCHRONIZED(1)



NAME
       perl subsample-synchronized.pl - Reduce the coverage of a synchronized file, by random subsampling, to the given target coverage

SYNOPSIS
       perl subsample-synchronized.pl --input input.sync --output output.sync --target-coverage 50 --max-coverage 2%  --method
       withoutreplacement

OPTIONS
       --input
           The input file in the synchronized format; Mandatory.

       --output
           The output file, will be a synchronized file  Mandatory.

       --target-coverage
           Reduce the coverage of the pileup-file to the here provided value; The target coverage also acts as minimum coverage, i.e.: if
           the coverage in any population is smaller than the targetcoverage the whole pileup entry is discarded. Mandatory

       --max-coverage
           The maximum coverage; All populations are required to have coverages lower or equal than the maximum coverage; Mandatory The
           maximum coverage may be provided as one of the following:
              '500' a maximum coverage of 500 will be used for all populations
              '300,400,500' a maximum coverage of 300 will be used for the first population, a maximum coverage of 400 for the second population and so on
            '2%' the 2% highest coverages will be ignored, this value is independently estimated for every population

       --method
           Specify the method for subsampling of the synchronized file. Either: withreplace, withoutreplace, fraction; Mandatory

            withreplace: subsample with replacement
            withoutreplace: subsample without replacement
            fraction: calculate the exact fraction of the allele frequencies and linearly scale them to the C<--target-coverage> and rounding to an integer;

       --help
           Display help for this script

Details

 Input synchronized
       A synchronized file

   Output
       The output will be a reduced coverage synchronized file
```
subsample_synchronized.pl

Targer coverage? 50?????
Ian with a shrug " Do 40".


Max-coverage? = 200 ??

Methods?
```
perl /usr/local/popoolation/subsample-synchronized.pl --input episodic_data_main.sync --output episodic_data_subsample.sync --target-coverage 40 --max-coverage 200 --method withreplace
```

Test this with screen
From mpileup dir for BWA files!
```
screen -r
script subsample_sync_log.log
perl /usr/local/popoolation/subsample-synchronized.pl --input episodic_data_main.sync --output episodic_data_subsample.sync --target-coverage 40 --max-coverage 200 --method withreplacement
```
error
unknown method for subsampling withreplacement at /usr/local/popoolation/Modules/Subsample.pm line 39.
withoutrepacement is wrong: should be withreplace (mistake on help page)

```
screen -r
script subsample_sync_log.log
perl /usr/local/popoolation/subsample-synchronized.pl --input episodic_data_main.sync --output episodic_data_subsample.sync --target-coverage 40 --max-coverage 200 --method withreplace
```

Possibly to high with 40? looks like everything is constant???
Try a CMH Test: 
Located in thing_test



```
#! /bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
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
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

# Can change here to other comparisons

pop[0]=1-3,2-4

cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl

for population in ${pop[@]}
do
mkdir ${mpileup_dir}/${population}_test
pop_dir=${mpileup_dir}/${population}_test

perl ${cmh_test} --min-count 3 --min-coverage 5 --max-coverage 100 --population ${population} --input ${mpileup_dir}/episodic_data_subsample_target40.sync --output ${pop_dir}/episodic_data_target40_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/episodic_data_target40_${population}.cmh.txt --output ${pop_dir}/episodic_data_target40_${population}.cmh.gwas --min-pvalue 1.0e-40

done
```

```
screen -r
script CMH_target40_13_24.log
CMH_13_24_target40_script
```

```
exit
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/1-3,2-4_test/episodic_data_target40_1-3,2-4.cmh.gwas /Users/paulknoops/Sequence_analysis_2016

```
```
java -Xmx2g -jar /Users/paulknoops/igv/IGV_2.3.67/igv.jar

```

Idea: make many sync files for comparisons of interst: likely need to do after subsample.syn

$0 prints all again.. just want to prink populations of interst
- will need to change print...., the output... and possibly the OFS=.....

```
cat ${mpileup_dir}/${project_name}.sync | awk 'BEGIN{OFS="\t"}{print $0,$13,$13,$13}' > ${mpileup_dir}/${project_name}_MGD2.sync
```

Not sure how to work with it yet.

Just run CMH tests (or Fst) with what I have


# It has been a few months: back at it:
### Recap:
1) have many sync files with merged base generation, removed unneeded areas and can can call the base generation multiple times

BWA
```
episodic_data_main.sync
```
Bowtie
```
episodic_data_bowtie_main.sync
```
Layout for BWA
```
F115ConR1_TAGCTT_merged_aligned_pe.final.bam
F115ConR2_GGCTAC_merged_aligned_pe.final.bam
F115SelR1_GTTTCG_merged_aligned_pe.final.bam
F115SelR2_GTGGCC_merged_aligned_pe.final.bam
F38ConR1_ATCACG_merged_aligned_pe.final.bam
F38ConR2_TTAGGC_merged_aligned_pe.final.bam
F38SelR1_ACTTGA_merged_aligned_pe.final.bam
F38SelR2_GATCAG_merged_aligned_pe.final.bam
F77ConR1_ATGTCA_merged_aligned_pe.final.bam
F77ConR2_ATTCCT_merged_aligned_pe.final.bam
F77SelR1_TTAGGC_merged_aligned_pe.final.bam
F77SelR2_GATCAG_merged_aligned_pe.final.bam
MGD3_SO_CAGATC_merged_aligned_pe.final.bam
```
Layout for Bowtie
```
F115ConR1_TAGCTT_merged_bowtie_pe.final.bam
F115ConR2_GGCTAC_merged_bowtie_pe.final.bam
F115SelR1_GTTTCG_merged_bowtie_pe.final.bam
F115SelR2_GTGGCC_merged_bowtie_pe.final.bam
F38ConR1_ATCACG_merged_bowtie_pe.final.bam
F38ConR2_TTAGGC_merged_bowtie_pe.final.bam
F38SelR1_ACTTGA_merged_bowtie_pe.final.bam
F38SelR2_GATCAG_merged_bowtie_pe.final.bam
F77ConR1_ATGTCA_merged_bowtie_pe.final.bam
F77ConR2_ATTCCT_merged_bowtie_pe.final.bam
F77SelR1_TTAGGC_merged_bowtie_pe.final.bam
F77SelR2_GATCAG_merged_bowtie_pe.final.bam
MGD3_SO_CAGATC_merged_bowtie_pe.final.bam
```

2) need to rerun all CMH tests and possibly Fst values:

- Should removed / move all old CMH outputs and start fresh to be sure everything is working properly.
- to clean up - moving all those with "target_40" from BWA diretory (keep script in scripts, CMH outputs in "old CMH" dir, )
- old scripts work??? possibly, but edit existing anyway to be sure on directories etc.

____________________
First set up smaller sync file for R practice analysis:
- only 3R (personal favourite)

- Also do 2R -- looks smaller region for BWA

- Need to do with both eventually, starting with just BWA -mem (better according to Scholotterer)

- Doing with more edited and proper .sync files

- Sync files have all het, U and mitochondria removed

From mpileup dir
```
grep '3R' episodic_data_main.sync > episodic_data_3R.sync
```

```
grep '2R' episodic_data_main.sync > episodic_data_2R.sync
```

Move to local (then copy/move to git-linked dir)
  - Check the size to make sure not to large for analysis!
  - Doing 2R (smaller by ~ a Gig == 3.6 G)
  
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/episodic_data_2R.sync /Users/paulknoops/Bioinformatics/Sequence_analysis_2016
```
Need to make smaller!

```
wc -l episodic_data_2R.sync
```
```
20547525 episodic_data_2R.sync
```
For me: random selection: 
Try different sizes

10000 bp
```
sed -n ' 10268762, 10278762 p' episodic_data_2R.sync > episodic_data_2R_subset.sync
```
File size = 1.8MB

1000bp
```
sed -n ' 10273262, 10274262 p' episodic_data_2R.sync > episodic_data_2R_1000_subset.sync
```
file size = 184 K
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/episodic_data_2R_subset.sync /Users/paulknoops/Bioinformatics/Sequence_analysis_2016 
```

```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/episodic_data_2R_1000_subset.sync /Users/paulknoops/Bioinformatics/Sequence_analysis_2016
```
Playing on R: see R script "sync_file.R"



Break the 2R for R
```
wc -l episodic_data_2R.sync
```
```
20547525 episodic_data_2R.sync
```

10000 bp == File size = 1.8MB
example:
```
sed -n ' 10268762, 10278762 p' episodic_data_2R.sync > episodic_data_2R_subset.sync
```
Try 100000
20547525/20 == 1027376

0,1027376,2054752,3082128,4109504,5136880,6164256,7191632,8219008,9246384,10273760,11301136,12328512,13355888,14383264,15410640,16438016,17465392,18492768,19520144,20547525

input = ../episodic_data_2R.sync
Run from outputs: /home/paul/episodicData/mpileup_dir/episodic_data2R_subsetting/
```
sed -n ' 1, 1027376 p' ../episodic_data_2R.sync > episodic_data_2R_1.sync
```
== 179 MB

Try One up ==360MB (remove the one above)
0,2054752,4109504,6164256,8219008,10273760,12328512,14383264,16438016,18492768,20547525
```
sed -n ' 1, 2054752 p' ../episodic_data_2R.sync > episodic_data_2R_1.sync
```

```
#scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/episodic_data2R_subsetting/episodic_data_2R_1.sync ~/MyDocuments/episodic_local_work
```
After filtering: still 22 Mb

#Next section:
```
sed -n ' 2054753, 4109504 p' ../episodic_data_2R.sync > episodic_data_2R_2.sync
```

The rest:
```
sed -n ' 4109505, 6164256 p' ../episodic_data_2R.sync > episodic_data_2R_3.sync
sed -n ' 6164257, 8219008 p' ../episodic_data_2R.sync > episodic_data_2R_4.sync
sed -n ' 8219009, 10273760 p' ../episodic_data_2R.sync > episodic_data_2R_5.sync
sed -n ' 10273761, 12328512 p' ../episodic_data_2R.sync > episodic_data_2R_6.sync
sed -n ' 12328513, 14383264 p' ../episodic_data_2R.sync > episodic_data_2R_7.sync
sed -n ' 14383265, 16438016 p' ../episodic_data_2R.sync > episodic_data_2R_8.sync
sed -n ' 16438017, 18492768 p' ../episodic_data_2R.sync > episodic_data_2R_9.sync
sed -n ' 18492769, 20547525 p' ../episodic_data_2R.sync > episodic_data_2R_10.sync
```
From local
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/mpileup_dir/episodic_data2R_subsetting/*.sync ~/MyDocuments/episodic_local_work
```


Want to look at uniform coverage! says slide 16: http://drrobertkofler.wikispaces.com/file/view/pooledAnalysis_part2.pdf/489488298/pooledAnalysis_part2.pdf



Viewing in IGV on brians machine:
Load IGV (java -Xmx750m -jar igv.jar) from igv direcotry on machine
"current" in brians machine??/usr/local/igv/current

```
java -Xmx2g -jar /usr/local/igv/current/igv.jar
```
Need to index bam files:
```
samtools index *.bam
#Do individually
```
Odd error with F115ConR2.....? something with SAM file....

Same with F77: fixed with 115, had another copy that could be opened
samtools index: failed to open "F77ConR2_ATTCCT_merged_aligned_pe.final.bam": Exec format error
Error looked to be small file size, when edditing and moving, must have removed the real file: run through again ?
Had an extra copy of F115, but not F77; need to run F77 again (from Trim_dir)
Run Regular scripts; just make a new trim dir with a copu of F77ConR2.



### Test for running in paralle: have two data sets for F115ConR1 and F115ConR2 bam files to be merged
Run each in parallel:
Make each into seperate scripts (F115ConR1_merge.sh and F115ConR2_merge.sh)
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_merge=${project_dir}/novo_merge

samtools merge ${novo_merge}/F115ConR1_TAGCTT_novo_merge.bam ${novo_bam}/F115ConR1_TAGCTT_L001_novo.bam ${novo_bam}/F115ConR1_TAGCTT_L002_novo.bam
```

```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_merge=${project_dir}/novo_merge

samtools merge ${novo_merge}/F115ConR2_GGCTAC_novo_merge.bam ${novo_bam}/F115ConR2_GGCTAC_L001_novo.bam ${novo_bam}/F115ConR2_GGCTAC_L002_novo.bam
```

Make new script == Parallel_merge.sh
```
F115ConR1_merge.sh &
F115ConR2_merge.sh &

wait
echo "both merged"
```


NOTE: editing names (remove A_ and add it to the back): for file in *; do mv "${file}" "${file/A_/}_A"; done

### Running GATK on final.bam files for BWA-mem?

For BWA-mem; can add a readgroup at that time; (look into adding for new run through)

1) Need an unzipped version of reference genome (make sure unzipped -- gunzip)

2) make a gatk directory (mkdir gatk_dir)

3) need to make sure the index directory has a .dict (done already for novoalign practice == location of script below) (novo_dict_index.sh*)

	- Note: copied reference to _2 at end to have one zipped copy and one unzipped
```
#! /bin/bash

pic=/usr/local/picard-tools-1.131/picard.jar
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

java -jar ${pic} CreateSequenceDictionary R=${ref_genome} O=${index_dir}/dmel-all-chromosome-r5.57_2.dict

```

4) Need Read Groups for GATK to run: So far as I can tell, the read groups can be anything, they just need to be there; Can edit them after the fact
  - RGID --Read Group Identifier; for Illumina, are composed using the flowcell + lane name and number [using Lanes L001_L002 for now]
  - RGLB -- DNA Preperation Library Identifier [library1 as place holder]
  - RGPL - platform/technology used to produce the read [Illumina]
  - RGPU -- Platform Unit; details on the sequencing unit (i.e run barcode) [None, used for practice]
  - RGSM -- Sample [Using the basename which is each unique sequence]
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to .bam files
novo_final=${project_dir}/final_bam

files=(${novo_final}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${pic} AddOrReplaceReadGroups I=${novo_final}/${base}.bam O=${novo_final}/${base}_RG.bam RGID=L001_L002 RGLB=library1 RGPL=illumina RGPU=None RGSM=${base}

done
```

5) Run GATK indelrealigner:

```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData

#Path to input directory
final_bam=${project_dir}/final_bam

#Path to output directory
gatk_dir=${project_dir}/gatk_dir

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

java -Xmx8g -jar ${gatk} -I ${final_bam}/${base}_RG.bam -R ${ref_genome} -T RealignerTargetCreator -o ${gatk_dir}/${base}.intervals

java -Xmx8g -jar ${gatk} -I ${final_bam}/${base}_RG.bam -R ${ref_genome} -T IndelRealigner -targetIntervals ${gatk_dir}/${base}.intervals -o ${gatk_dir}/${base}_realigned.bam

done
```
### Delay running GATK for now; may either change parameters, may not be working, or may ned to use Unified Genotyper

## Testing some variant callers:
__________________
### CRISP
Getting CRISP onto Brians Machine:
1) download crisp onto local machine (.tar.gz): download from https://bansal-lab.github.io/software/crisp.html
2) SCP to remote location (``` scp CRISP-122713.tar.gz paul@info.mcmaster.ca:/home/paul ```)
3) Unpack file (```  tar xvzf CRISP-122713.tar.gz ```) 

Running CRISP: 

For Novoalign files test:

https://bansal-lab.github.io/software/crisp.html

https://github.com/vibansal/crisp

Example from Bergland:
 
 ```
  	baseDir=/mnt/Alan_Backup/bam
 	for i in 2L 2R 3L 3R X; do
 		for j in $(seq 1 1500000 30000000); do
 			CRISP --bam $baseDir/FL_rep1.sort.rmdup.realign.bam \
 				--bam $baseDir/FL_rep2.sort.rmdup.realign.bam \
 				--bam $baseDir/GA.sort.rmdup.realign.bam \
 				--bam $baseDir/SC.sort.rmdup.realign.bam \
 				--bam $baseDir/NC.sort.rmdup.realign.bam \
 				--bam $baseDir/ME_rep1.sort.rmdup.realign.bam \
 				--bam $baseDir/ME_rep2.sort.rmdup.realign.bam \
 				--bam $baseDir/PA_7_2009_run2.sort.rmdup.realign.bam \
 				--bam $baseDir/PA_11_2009.sort.rmdup.realign.bam \
 				--bam $baseDir/PA_7_2010.sort.rmdup.realign.bam \
 				--bam $baseDir/PA_11_2010.sort.rmdup.realign.bam \
 				--bam $baseDir/PA_7_2011.merged_run1_run2.rmdup.realign.bam \
 				--bam $baseDir/PA_10_2011.sort.rmdup.realign.bam \
 				--bam $baseDir/PA_11_2011.merged_run1_run2.rmdup.realign.bam \
 				--ref /mnt/Alan/dmel_reference/all_dmel.fasta \
 				--poolsize 100 \
 				--perms 1000 \
 				--filterreads 0 \
 				--regions $i:$j-$((j+1500000-1)) \
 				--qvoffset 33 \
 				--mbq 10 \
 				--mmq 10 \
 				--minc 4 \
 				--VCF /mnt/Alan/new_alignments/vcf2/6d_v7.2.crisp.vcf_$i\_$j > 6d_v7.2.log_$i\_$j & 
 		done
 		wait
 	done
 ```
```
./CRISP [options] --bams file_bam_paths --ref reference.fasta --VCF variantcalls.VCF --poolsize poolsize --bed targets.bed > variantcalls.log
```
Need a direcory for outputs:

mkdir novo_crisp

Flags when running CRISP:
```
--bams/bam: bams == textfile with list of bam file paths (one for each pool), bam == bam file for one pool, specify filename for each pool using --bam pool1.bam --bam pool2.bam .... --bam pooln.bam

--ref: Indexed Reference Sequence file (fasta)

--poolsize: poolsize (number of haploid genomes in each pool), for diploid genomes: 2 x # individuals [120]

--perms: maximum number of permutations for calculating contingency table p-value, default 20,000

--filterreads: filter reads with excessive number of mismatches/gaps compared to the reference sequence, Default==1, 0 to disable

--regions: region(s) in which variants will be called, e.g chr1:654432-763332. BAM files should be indexed for this option with the index file pooln.bam.bai 

--qvoffset: quality value offset, 33 for Sanger format, 64 for Illumina 1.3+ format

--mbq: minimum base quality to consider a base for variant calling, default 10

--mmq: minimum read mapping quality to consider a read for variant calling, default 20

--minc: minimum number of reads with alternate allele required for calling a variant, default 4

--ctpval: threshold on the contingency table p-value for calling position as variant (specified as log10), default is -3.5, increase this threshold as number of pools increases

--qvpval: threshold on the quality values based p-value for calling position as variant (specified as log10), default is -5

--bed: bedfile with intervals in which variants will be called

--VCF: VCF file to which the variant calls will be output 
```

Need a list of all bam files with path names: script to call files from .bam file dir (realigned with GATK better)

```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .bam files from GATK
novo_gatk=${project_dir}/novo_GATK


files=(${novo_gatk}/*_realigned.bam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _realigned.bam`
echo "${novo_gatk}/${base}_realigned.bam" >> ${novo_gatk}/novo_list.bam

done
```

Running CRISP:

```
#! /bin/bash

#Variable for project name (file name)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to CRISP
crisp=/home/paul/CRISP-122713/CRISP

#Variable for reference genome (non-zipped) ## .faidx indexed with samtools? .fai present
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

#Path to .bam files from GATK
novo_gatk=${project_dir}/novo_GATK

#Output
novo_crisp=${project_dir}/novo_crisp

#files=(${novo_gatk/*_novo_merge_novo_final_realigned.bam)
#for file in ${files[@]}
#do
#name=${file}
#base=`basename ${name} _novo_merge_novo_final_realigned.bam`

${crisp} --bams ${novo_gatk}/novo_list.bam \
			--ref ${ref_genome} \
 				--poolsize 120 \
 				--perms 1000 \
 				--filterreads 0 \
 				--qvoffset 33 \
 				--mbq 10 \
 				--mmq 10 \
 				--minc 4 \
 				--VCF ${novo_crisp}/${project_name}.vcf > ${project_name}_variantcalls.log
```

Example output in VCF file

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  F115ConR1_TAGCTT_novo_merge_novo_final_realigned        F115ConR2_GGCTAC_novo_merge_novo_final_realigned

2L      891237  .       C       A       1844    LowDepth        NP=2;DP=54,61,2;VT=SNV;CT=-inf;VP=2;VF=EMpass;AC=121;AF=0.50515;EMstats=184.47:-89.69;HWEstats=-0.0;MQS=0,0,0,117;FLANKSEQ=tactaatctc:C:atatcaacat      MLAC:GQ:DP:ADf:ADr:ADb  .:0:72:12,23:17,19:1,0  .:0:45:11,8:16,9:0,1   
```



### Varscan> 
-- 
1) download latest version of Varscan (VarScan.v2.3.9.jar --> https://sourceforge.net/projects/varscan/files/)
2) scp onto Brians machine: (``` scp VarScan.v2.3.9.jar paul@info.mcmaster.ca:/home/paul   ```)
3) scp the txt file onto machine: ( ``` scp VarScan.v2.3.9.description.txt paul@info.mcmaster.ca:/home/paul   ```)
4) Can now run in command line with java -jar VarScan.jar

-- Manual http://dkoboldt.github.io/varscan/doc/index.html

mpileup2snp (following Huang et al. Evaluation of variant detection software for pooled next-generation sequence data)


### SnpEff>
--



### LoFreq>

https://github.com/CSB5/lofreq

http://csb5.github.io/lofreq/commands/

Downloading LoFreq:
1) download most recent loFreq onto local machine (lofreq_star-2.1.2_linux-x86-64.tgz): download from https://sourceforge.net/projects/lofreq/files/
2) SCP to remote location (``` scp lofreq_star-2.1.2_linux-x86-64.tgz paul@info.mcmaster.ca:/home/paul ```)
3) Unpack file (```  tar xvzf lofreq_star-2.1.2_linux-x86-64.tgz	 ```) 

Running LoFreq:


> Huang et al. “Evaluation of variant detection software for pooled next-generation sequence data”:

*Lofreq version 2.0.0-rc-1 (commit 2.0.0-rc-1-3-g63449f7) was installed and run in “call” mode.*

"One program, LoFreq, only permitted the analysis of one pool of samples at a time, and another, CRISP, would only run on groups of pools."

> From LoFreq Paper (Wilm et al. 2012)

LoFreq

LoFreq takes a samtools pileup as input (samtools mpileup; Version 0.1.18). By default samtools applies a coverage cap and we set this to be sufficiently high to avoid filtering reads in a sample (-d 100000). Whenever indels were not allowed for read mapping, we switched off samtools BAQ computation (-B). SNVs were called with a Bonferroni-corrected P-value threshold of 0.05 and the same threshold was applied for calling somatic variants with the binomial test. Unless stated otherwise, we removed variant positions with a significant strand bias (Holm–Bonferroni-corrected P-value < 0.05) from LoFreq predictions.

> Confusing: most documentation uses input as a .bam file (sems to be a single .bam), but the paper (and only from the paper) is mpileup mentioned....https://github.com/CSB5/lofreq/issues/32

Need to run one per .bam file?
Create script to create individual scripts per Bam file (similar to Novoalign scripts to run novoalign seperate)

### Snape>
--


### Dropping the packages Test:
Just using the GATK indel Realigner and CRISP for now: get them going on bowtie and BWA-mem packages:

## BWA GatK:
Copied from Above:

Running GATK on final.bam files for BWA-mem?

For BWA-mem; can add a readgroup at that time; (look into adding for new run through)

1) Need an unzipped version of reference genome (make sure unzipped -- gunzip)

2) make a gatk directory (mkdir gatk_dir)

3) need to make sure the index directory has a .dict (done already for novoalign practice == location of script below)
```
#! /bin/bash

pic=/usr/local/picard-tools-1.131/picard.jar
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

java -jar ${pic} CreateSequenceDictionary R=${ref_genome} O=${index_dir}/dmel-all-chromosome-r5.57_2.dict

```

4) Need Read Groups for GATK to run: So far as I can tell, the read groups can be anything, they just need to be there; Can edit them after the fact
  - RGID --Read Group Identifier; for Illumina, are composed using the flowcell + lane name and number [using Lanes L001_L002 for now]
  - RGLB -- DNA Preperation Library Identifier [library1 as place holder]
  - RGPL - platform/technology used to produce the read [Illumina]
  - RGPU -- Platform Unit; details on the sequencing unit (i.e run barcode) [None, used for practice]
  - RGSM -- Sample [Using the basename which is each unique sequence]
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to .bam files
input=${project_dir}/final_bam

files=(${input}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${pic} AddOrReplaceReadGroups I=${input}/${base}.bam O=${input}/${base}_RG.bam RGID=L001_L002 RGLB=library1 RGPL=illumina RGPU=None RGSM=${base}

done
```
5) 
Index RG.bams
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData

#Path to input directory
input=${project_dir}/final_bam

files=(${input}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`
samtools index ${input}/${base}_RG.bam
done
```

6) Run GATK indelrealigner:

```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData

#Path to input directory
input=${project_dir}/final_bam

#Path to output directory
output=${project_dir}/gatk_dir

#Variable for reference genome (non-zipped)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar


files=(${input}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

java -Xmx8g -jar ${gatk} -I ${input}/${base}_RG.bam -R ${ref_genome} -T RealignerTargetCreator -o ${output}/${base}.intervals

java -Xmx8g -jar ${gatk} -I ${input}/${base}_RG.bam -R ${ref_genome} -T IndelRealigner -targetIntervals ${output}/${base}.intervals -o ${output}/${base}_realigned.bam

done
```

## Bowtie2 GatK

1) Create GATK output dir2
2) Read Groups
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/bowtie

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to .bam files
input=${project_dir}/final_bam

files=(${input}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${pic} AddOrReplaceReadGroups I=${input}/${base}.bam O=${input}/${base}_RG.bam RGID=L001_L002 RGLB=library1 RGPL=illumina RGPU=None RGSM=${base}

done
```
3) Index RG.bams
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/bowtie

#Path to input directory
input=${project_dir}/final_bam

files=(${input}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`
samtools index ${input}/${base}_RG.bam
done
```

4) Run Gatk
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/bowtie

#Path to input directory
input=${project_dir}/final_bam

#Path to output directory
output=${project_dir}/gatk_bowtie

#Variable for reference genome (non-zipped)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar


files=(${input}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

java -Xmx8g -jar ${gatk} -I ${input}/${base}_RG.bam -R ${ref_genome} -T RealignerTargetCreator -o ${output}/${base}.intervals

java -Xmx8g -jar ${gatk} -I ${input}/${base}_RG.bam -R ${ref_genome} -T IndelRealigner -targetIntervals ${output}/${base}.intervals -o ${output}/${base}_realigned.bam

done
```

### CRISP
Change the hash between the two;

```
#! /bin/bash

#Variable for project name (file name)

#project_name=episodic_data
project_name=episodic_data_bowtie



#Variable for project:

#project_dir=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie


#Path to CRISP
crisp=/home/paul/CRISP-122713/CRISP

#Variable for reference genome
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

#Path to .bam files from GATK

#input=${project_dir}/gatk_dir
input=${project_dir}/gatk_bowtie


#Output

#output=${project_dir}/CRISP
output=${project_dir}/CRISP_Bowtie


files=(${input}/*.bam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
echo "${input}/${base}.bam" >> ${input}/${project_name}_BAMlist.txt

done

${crisp} --bams ${input}/${project_name}_BAMlist.txt \
			--ref ${ref_genome} \
 			--poolsize 120 \
 			--perms 1000 \
 			--filterreads 0 \
 			--qvoffset 33 \
 			--minc 4 \
 			--VCF ${output}/${project_name}.vcf > ${output}/${project_name}_variantcalls.log
```
Flags:
```
--bams/bam: bams == textfile with list of bam file paths (one for each pool)
--ref: Indexed Reference Sequence file (fasta)
--poolsize: poolsize (number of haploid genomes in each pool), for diploid genomes: 2 x # individuals [pool = 60, == 120]
--perms: maximum number of permutations for calculating contingency table p-value, default 20,000 [Bergland example used 1000, keeps run times lower for this many samples]
--filterreads: filter reads with excessive number of mismatches/gaps compared to the reference sequence, Default==1, 0 to disable [0, filtered enough before, keep what is left]
--qvoffset: quality value offset, 33 for Sanger format, 64 for Illumina 1.3+ format [33, Illumina (6) far above failed?]
--mbq: minimum base quality to consider a base for variant calling, default 10 [default < previous filtering]
--mmq: minimum read mapping quality to consider a read for variant calling, default 20 [default == previous filtering]
--minc: minimum number of reads with alternate allele required for calling a variant, default 4 [stick with default]
--ctpval: threshold on the contingency table p-value for calling position as variant (specified as log10), default is -3.5, increase this threshold as number of pools increases [not specified in Bergland, leave at deafault? == Bergland Pool == 100 (close enough)]
--qvpval: threshold on the quality values based p-value for calling position as variant (specified as log10), default is -5 [stick with default]
--VCF: VCF file to which the variant calls will be output 
```

In appropriate scripts directory
```
nano bwa_Crisp.sh
nano bowtie_Crisp.sh
```
Run on own screen with logs. (Running on Info115)
________________________
Using VCFtools???
```
/usr/local/vcftools --vcf episodic_data.vcf
/usr/local/vcftools --vcf episodic_data_bowtie.vcf
#Top two do not work: make sure path has the vcftools within (this is src, also a perl option??)

/usr/local/vcftools/src/cpp/vcftools --vcf episodic_data_bowtie.vcf
/usr/local/vcftools/src/cpp/vcftools --vcf episodic_data.vcf
```

Bergland Script for After .VCF made with CRISP

```
## merge subvcf files

	# extract header
	head -n200 $(ls 6d_v7.2.crisp* | head -n1) | grep "#" > vcf.header
	
	# strip headers
	nHeaderLines=$(wc -l vcf.header | grep -Eo "[0-9]+")
	sed -i "1,$nHeaderLines d" 6d_v7.2.crisp*
	
	# merge
	cat vcf.header 6d_v7.2.crisp* > 6d_v7.2.crisp.vcf
	
	# move sub bits for backup
	mkdir vcfBits
	mv 6d_v7.2.crisp.vcf_* vcfBits/
	mv 6d_v7.2.log* vcfBits/


## remove SNPs in repetative regions
	/home/alan/vcftools_0.1.7/bin/vcftools --vcf 6d_v7.2.crisp.vcf --out 6d_v7.2.crisp.norepeats --exclude-bed /mnt/Alan/new_alignments/vcf2/repetative/rm.edit.bed --recode --recode-INFO-all
	mv 6d_v7.2.crisp.norepeats.recode.vcf 6d_v7.2.crisp.norepeats.vcf
	
# ## run through awk converstion script; this strips out triallelic sites, calculated average allele freq and converts population columns to AF:DP format
	### here is contents of crisp_vcf_conv.awk:
	
	# 		{
# 			# Copy file header unchanged
# 				if(substr($1, 1, 1)=="#") {
# 					print $0
# 				} else {
# 					if(match($5, ",")==0) {
# 						# Maintains first 7 columns (pos, chr, etc.)			
# 						for(i=1; i<=7; i++) {
# 							printf $i"\t"
# 						}
# 						
# 						# get total allele frequency
# 						np = 0
# 						af = 0
# 						for(i=10; i<=NF; i++) {
# 							split($i, alleleCounts, ":")
# 							split(alleleCounts[3], forwardCounts, ",")
# 							split(alleleCounts[4], reverseCounts, ",")
# 							
# 							altCounts[i-9] = forwardCounts[2] + reverseCounts[2]
# 							rd[i-9] = forwardCounts[1] + reverseCounts[1] + forwardCounts[2] + reverseCounts[2]
# 							if(rd[i-9]>0) {
# 								af += altCounts[i-9]/rd[i-9]
# 								np++
# 							}
# 						}
# 						
# 						# print info column
# 						printf "AF="af/np";"$8"\t"
# 						
# 						# print format column
# 						printf "AF:DP\t"
# 						
# 						#print population columns
# 						for(i=10; i<=NF; i++) {
# 							if(rd[i-9] > 0) printf altCounts[i-9]/rd[i-9]":"rd[i-9]
# 							if(rd[i-9] == 0) printf "0:0"
# 							if(i<NF) printf "\t"
# 							if(i==NF) printf "\n"
# 						}
# 					}
# 				}
# 			}
				
	awk -f /mnt/Alan/new_alignments/vcf2/helperScripts/crisp_vcf_conv.awk < 6d_v7.2.crisp.norepeats.vcf > 6d_v7.2.norepeats.vcf
	wait
	
## split SNPs and INDELs, and flag SNPs wihtin 5bp of indel
	awk -F ';VT=SNV;' '{if(NF==2) print > "6d_v7.2.norepeats.snps.vcf"; if(NF==1) print > "6d_v7.2.norepeats.indels.vcf" }' < 6d_v7.2.norepeats.vcf
	sort -k1,1 -k2,2n -T /mnt/Alan_Backup --parallel 10 6d_v7.2.norepeats.snps.vcf > 6d_v7.2.norepeats.snps.vcf.sort
	rm 6d_v7.2.norepeats.snps.vcf 
	mv 6d_v7.2.norepeats.snps.vcf.sort 6d_v7.2.norepeats.snps.vcf

	sort -k1,1 -k2,2n -T /mnt/Alan_Backup --parallel 10 6d_v7.2.norepeats.indels.vcf > 6d_v7.2.norepeats.indels.vcf.sort
	rm 6d_v7.2.norepeats.indels.vcf 
	mv 6d_v7.2.norepeats.indels.vcf.sort 6d_v7.2.norepeats.indels.vcf

	
	/usr/bin/closestBed -D ref -t first -a 6d_v7.2.norepeats.snps.vcf -b 6d_v7.2.norepeats.indels.vcf | awk '{for(i=1; i<=7; i++) printf $i"\t"; printf "D2I="$NF";"$8"\t"$9"\t"; for(i=10; i<=23; i++) {printf $i; if(i<23) printf "\t"; if(i==23) printf "\n" }}' > 6d_v7.2.d2i.norepeats.vcf
	wait
```


My Modifactions: 

1) awk converstion script; this strips out triallelic sites, calculated average allele freq and converts population columns to AF:DP format (Completely from Bergland)
```	
{
 	# Copy file header unchanged
			if(substr($1, 1, 1)=="#") {
 				print $0
 			} else {
 				if(match($5, ",")==0) {
 					# Maintains first 7 columns (pos, chr, etc.)			
 					for(i=1; i<=7; i++) {
						printf $i"\t"
 					}
 						
 			# get total allele frequency
 				np = 0
				af = 0
 				for(i=10; i<=NF; i++) {
 					split($i, alleleCounts, ":")
 					split(alleleCounts[3], forwardCounts, ",")
 					split(alleleCounts[4], reverseCounts, ",")
 						
 						altCounts[i-9] = forwardCounts[2] + reverseCounts[2]
 						rd[i-9] = forwardCounts[1] + reverseCounts[1] + forwardCounts[2] + reverseCounts[2]
 						if(rd[i-9]>0) {
 						af += altCounts[i-9]/rd[i-9]
						np++
							}
						}
				
						# print info column
 						printf "AF="af/np";"$8"\t"
 						
 						# print format column
 						printf "AF:DP\t"
 						
 						#print population columns
 						for(i=10; i<=NF; i++) {
 							if(rd[i-9] > 0) printf altCounts[i-9]/rd[i-9]":"rd[i-9]
 							if(rd[i-9] == 0) printf "0:0"
 							if(i<NF) printf "\t"
 							if(i==NF) printf "\n"
 						}
 					}
 				}
 			}

```

2) Run Awk script:
```

# awk -f /mnt/Alan/new_alignments/vcf2/helperScripts/crisp_vcf_conv.awk < 6d_v7.2.crisp.norepeats.vcf > 6d_v7.2.norepeats.vcf

#BWA
awk -f /home/paul/episodicData/CRISP/crisp_vcf_conv.awk < episodic_data.vcf  > episodic_data.norepeats.vcf

#Bowtie
awk -f /home/paul/episodicData/CRISP/crisp_vcf_conv.awk < episodic_data_bowtie.vcf  > episodic_data_bowtie.norepeats.vcf

```

5) split SNPs and INDELs, and flag SNPs wihtin 5bp of indel

```
awk -F ';VT=SNV;' '{if(NF==2) print > "episodic_data.norepeats.snps.vcf"; if(NF==1) print > "episodic_data.norepeats.indels.vcf" }' < episodic_data.norepeats.vcf
	
	#SNPs
	#-T == use temporary directory???
	
	sort -k1,1 -k2,2n -T /home/paul/episodicData/CRISP/Paul_Backup episodic_data.norepeats.snps.vcf > episodic_data.norepeats.snps.vcf.sort
	
	#rm episodic_data.norepeats.snps.vcf
	
	#mv episodic_data.norepeats.snps.vcf.sort episodic_data.norepeats.snps.vcf

	#Indels
	
	sort -k1,1 -k2,2n -T /home/paul/episodicData/CRISP/Paul_Backup episodic_data.norepeats.indels.vcf > episodic_data.norepeats.indels.vcf.sort
	
	#rm episodic_data.norepeats.indels.vcf
	
	#mv episodic_data.norepeats.indels.vcf.sort episodic_data.norepeats.indels.vcf

	#Does not work?
	
	/usr/bin/closestBed -D ref -t first -a episodic_data.norepeats.snps.vcf -b episodic_data.norepeats.indels.vcf | awk '{for(i=1; i<=7; i++) printf $i"\t"; printf "D2I="$NF";"$8"\t"$9"\t"; for(i=10; i<=23; i++) {printf $i; if(i<23) printf "\t"; if(i==23) printf "\n" }}' > episodic_data.norepeats.vcf

```

Have 2 .vcf files: one for SNPS and one for Indels:

Open Java (NEW VERSION)
```
java -Xmx2g -jar /Users/paulknoops/Bioinformatics/IGV_2.3.94.app/Contents/Java/igv.jar
```

Can make a mpileup and sync with the .bam realigned with GATK!

BWA/Bowtie -- Output to gatk_dir -- 

Change Hash between two
```
#! /bin/bash

## Variable for project name (file name)

#project_name=episodic_data
project_name=episodic_data_bowtie



## Variable for project:

#project_dir=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie

## Variable for reference genome

index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

## Path to .bam files from GATK

#gatk=${project_dir}/gatk_dir
gatk=${project_dir}/gatk_bowtie


sync=/usr/local/popoolation/mpileup2sync.jar

samtools mpileup -B -Q 0 -f ${ref_genome} ${gatk}/*.bam > ${gatk}/${project_name}.gatk.mpileup

```

Sync: Hash between!
```
#! /bin/bash

## Variable for project name (file name)

project_name=episodic_data
#project_name=episodic_data_bowtie

## Variable for project:

project_dir=/home/paul/episodicData
#project_dir=/home/paul/episodicData/bowtie

## Variable for reference genome

index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

## Path to .bam files from GATK

gatk=${project_dir}/gatk_dir
#gatk=${project_dir}/gatk_bowtie

sync=/usr/local/popoolation/mpileup2sync.jar

java -ea -Xmx7g -jar ${sync} --input ${gatk}/${project_name}.gatk.mpileup --output ${gatk}/${project_name}.gatk.sync --fastq-type sanger --min-qual 20 --threads 2
```

Remove unneeded areas (with 'Het' and 'U')

-- One script: in BWA Script dir

```
#! /bin/bash

## Variable for project name (file name)
project_name=episodic_data
project_name_2=episodic_data_bowtie

## Variable for project:
project_dir=/home/paul/episodicData
project_dir_2=/home/paul/episodicData/bowtie

## Path to .bam files from GATK
gatk=${project_dir}/gatk_dir
gatk_2=${project_dir_2}/gatk_bowtie

grep -v 'Het' ${gatk}/${project_name}.gatk.sync > ${gatk}/${project_name}_less_het.sync

wait

grep -v 'U' ${gatk}/${project_name}_less_het.sync > ${gatk}/${project_name}_removed_U_Het.sync

wait

grep -v 'dmel_mitochondrion_genome' ${gatk}/${project_name}_removed_U_Het.sync > ${gatk}/${project_name}_main.gatk.sync

wait

rm -f ${gatk}/${project_name}_less_het.sync

rm -f ${gatk}/${project_name}_removed_U_Het.sync

grep -v 'Het' ${gatk_2}/${project_name_2}.gatk.sync > ${gatk_2}/${project_name_2}_less_het.sync

wait

grep -v 'U' ${gatk_2}/${project_name_2}_less_het.sync > ${gatk_2}/${project_name_2}_removed_U_Het.sync

wait

grep -v 'dmel_mitochondrion_genome' ${gatk_2}/${project_name_2}_removed_U_Het.sync > ${gatk_2}/${project_name_2}_main.gatk.sync

wait


rm -f ${gatk_2}/${project_name_2}_less_het.sync

rm -f ${gatk_2}/${project_name_2}_removed_U_Het.sync

```

Random subset for work:
```
grep '2R' episodic_data_main.gatk.sync > episodic_data_2R.gatk.sync
```

```
grep '2R' episodic_data_bowtie_main.gatk.sync > episodic_data_bowtie_2R.gatk.sync
```

Move to test directory:
```
#mkdir /home/paul/episodicData/subsetting

mv /home/paul/episodicData/gatk_dir/episodic_data_2R.gatk.sync /home/paul/episodicData/subsetting
mv /home/paul/episodicData/bowtie/gatk_bowtie/episodic_data_bowtie_2R.gatk.sync /home/paul/episodicData/subsetting
```

Need to make smaller!
```
wc -l episodic_data_2R.gatk.sync
```
== 20547525
```
wc -l episodic_data_bowtie_2R.gatk.sync
```
== 20420127


For me: random selection: 10000 bp (different b/w bowtie and bwa???)
```
sed -n ' 10268762, 10278762 p' episodic_data_2R.gatk.sync > episodic_data_2R_subset.gatk.sync
sed -n ' 10268762, 10278762 p' episodic_data_bowtie_2R.gatk.sync > episodic_data_bowtie_2R_subset.gatk.sync
```
Test everything from here (both full and subset??)

Test CMH:
```
! /bin/bash

## Variable for project name (file name)

#project_name=episodic_data_2R_subset
project_name=episodic_data_bowtie_2R_subset

## Variable for project:

project_dir=/home/paul/episodicData/subsetting

# Can change here to other comparisons

pop[0]=11-13,12-13
pop[1]=1-13,2-13
pop[2]=1-3,2-4
pop[3]=3-13,4-13
pop[4]=5-13,6-13
pop[5]=5-7,6-8
pop[6]=7-13,8-13
pop[7]=9-11,10-12
pop[8]=9-13,10-13

cmh_test=/usr/local/popoolation/cmh-test.pl
cmh_gwas=/usr/local/popoolation/export/cmh2gwas.pl

for population in ${pop[@]}
do
mkdir ${project_dir}/${project_name}_${population}
pop_dir=${project_dir}/${project_name}_${population}

perl ${cmh_test} --min-count 3 --min-coverage 5 --max-coverage 100 --population ${population} --input ${project_dir}/${project_name}.gatk.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-40

mv ${pop_dir}/${project_name}_${population}.cmh.gwas ${project_dir}/gwas_dir

done
```
(in a pinch with two many directories (which I will have --> rm -rf ${dir} to force remove full directories)

```
scp paul@info.mcmaster.ca:/home/paul/episodicData/subsetting/gwas_dir/*.gwas /Users/paulknoops/Bioinformatics/2Rsubset_cmhGWAS
```
```
java -Xmx2g -jar /Users/paulknoops/Bioinformatics/IGV_2.3.94.app/Contents/Java/igv.jar
```

Subset files to SCP:
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/subsetting/episodic_data_2R_subset.gatk.sync /Users/paulknoops/Bioinformatics/episodic_practice
scp paul@info.mcmaster.ca:/home/paul/episodicData/subsetting/episodic_data_bowtie_2R_subset.gatk.sync /Users/paulknoops/Bioinformatics/episodic_practice

```

Getting Haploreconstruct onto brians machine:
```
#install.packages(path_to_file, repos = NULL, type="source")
install.packages("/home/paul/R-packages/haploReconstruct_0.1.2.tar.gz", repos=NULL,type="source") 
#Does not work

#try from command line
R CMD INSTALL haploReconstruct_0.1.2.tar.gz

#e-mailed brian
```

###FST values? old:
```
#! /bin/bash

## Variable for project name (file name)

#project_name=episodic_data_2R_subset
project_name=episodic_data_bowtie_2R_subset

## Variable for project:

project_dir=/home/paul/episodicData/subsetting



fst_test=/usr/local/popoolation/fst-sliding.pl
fst_igv=/usr/local/popoolation/export/pwc2igv.pl

perl ${fst_test} --window-size 500 --step-size 500 --suppress-noninformative --input ${project_dir}/${project_name}.gatk.sync --min-covered-fraction 1.0 --min-coverage 10 --max-coverage 250 --min-count 3 --output ${project_dir}/${project_name}_gatk.fst.txt --pool-size 60

#To view in IGV
perl ${fst_igv} --input ${project_dir}/${project_name}_gatk.fst.txt --output ${project_dir}/${project_name}_gatk.fst.igv
```

```
scp paul@info.mcmaster.ca:/home/paul/episodicData/subsetting/episodic_data_2R_subset_gatk.fst.igv /Users/paulknoops/Bioinformatics/episodic_practice
scp paul@info.mcmaster.ca:/home/paul/episodicData/subsetting/episodic_data_bowtie_2R_subset_gatk.fst.igv /Users/paulknoops/Bioinformatics/episodic_practice

java -Xmx2g -jar /Users/paulknoops/Bioinformatics/IGV_2.3.94.app/Contents/Java/igv.jar
```

#Fishers; old (failed) from Bio 720 (2 years ago!)

```
#! /bin/bash

## Variable for project name (file name)

project_name=episodic_data_2R_subset
#project_name=episodic_data_bowtie_2R_subset

## Variable for project:

project_dir=/home/paul/episodicData/subsetting

fisher=/usr/local/popoolation/fisher-test.pl
fisher_igv=/usr/local/popoolation/export/pwc2igv.pl

perl ${fisher} --input ${project_dir}/${project_name}.gatk.sync --output ${project_dir}/${project_name}.gatk.fet --min-count 3 --min-coverage 10 --max-coverage 250 --suppress-noninformative

perl ${fisher_igv} --input ${project_dir}/${project_name}.gatk.fet --output ${project_dir}/${project_name}.gatk.fet.igv

# Load to IGV java -Xmx2g -jar /usr/local/igv/IGV_2.1.21/igv.jar
```

Same subset??
```
[paul@info115 subsetting]$ sed -n ' 10268762 p' episodic_data_2R.gatk.sync 
2R      10635452        C       0:0:55:0:0:0    0:0:26:0:0:0    0:0:37:0:0:0    0:0:45:0:0:0    0:0:48:0:0:0    0:0:43:0:0:0    0:0:69:0:0:0  0:0:53:0:0:0    0:0:107:0:0:0   0:0:130:0:0:0   0:0:152:0:0:0   0:0:111:0:0:0   0:0:41:0:0:0
[paul@info115 subsetting]$ sed -n ' 10268762 p' episodic_data_bowtie_2R.gatk.sync 
2R      10755361        C       0:0:15:0:0:0    0:0:21:0:0:0    0:0:15:0:0:0    0:0:18:0:0:0    0:0:20:0:0:0    0:0:23:0:0:0    0:0:27:0:0:0  0:0:26:0:0:0    0:0:158:0:0:0   0:0:143:0:0:0   0:0:177:0:0:0   0:0:152:1:0:0   0:0:49:0:0:0
```
119909 diff == 10148853
```
[paul@info115 subsetting]$ sed -n ' 10148853 p' episodic_data_bowtie_2R.gatk.sync
2R      10635445        C       0:0:9:0:0:0     0:0:5:0:0:0     0:0:5:0:0:0     0:0:3:0:0:0     0:0:11:0:0:0    0:0:6:0:0:0     0:0:15:0:0:0  0:0:7:0:0:0     0:0:19:0:0:0    0:0:39:0:0:0    0:0:36:0:0:0    0:0:22:0:0:0    0:0:11:0:0:0
```
A little less but okay none the less: so re subset bowtie (delete old)
```
#Old:
#sed -n ' 10268762, 10278762 p' episodic_data_bowtie_2R.gatk.sync > episodic_data_bowtie_2R_subset.gatk.sync

#New:
sed -n ' 10148853, 10158853 p' episodic_data_bowtie_2R.gatk.sync > episodic_data_bowtie_2R_subset.gatk.sync

scp paul@info.mcmaster.ca:/home/paul/episodicData/subsetting/episodic_data_bowtie_2R_subset.gatk.sync /Users/paulknoops/Bioinformatics/episodic_practice
```
All the fishers, Fst and CMH are NOT the same place as BWA


### R stuff

Moved the episodic_main.gatk.sync to R_dir for creating a .csv file for R use

2 scripts: sync_to_counts_bowtie.R and sync_to_counts_bwa.R

Open screen and start a script,

Start R
```
R
```
and source the file (should require("tidyr") first to make sure it is in there).

Did not work: 
Do by chrsomosome

```
grep '3R' episodic_data_main.gatk.sync > episodic_data_3R.gatk.sync &
grep '2R' episodic_data_main.gatk.sync > episodic_data_2R.gatk.sync &
grep '3L' episodic_data_main.gatk.sync > episodic_data_3L.gatk.sync &
grep '2L' episodic_data_main.gatk.sync > episodic_data_2L.gatk.sync &
grep '^4' episodic_data_main.gatk.sync > episodic_data_4.gatk.sync &
grep 'X' episodic_data_main.gatk.sync > episodic_data_X.gatk.sync
```

```
grep '3R' episodic_data_bowtie_main.gatk.sync > episodic_data_bowtie_3R.gatk.sync &
grep '2R' episodic_data_bowtie_main.gatk.sync > episodic_data_bowtie_2R.gatk.sync &
grep '3L' episodic_data_bowtie_main.gatk.sync > episodic_data_bowtie_3L.gatk.sync &
grep '2L' episodic_data_bowtie_main.gatk.sync > episodic_data_bowtie_2L.gatk.sync &
grep '^4' episodic_data_bowtie_main.gatk.sync > episodic_data_bowtie_4.gatk.sync &
grep 'X' episodic_data_bowtie_main.gatk.sync > episodic_data_bowtie_X.gatk.sync
```
Make a script for each arm of the chromosome to write a .csv file:



Open screen and source package: do each in unison

Start R
```
R
```
and source the file (should require("tidyr") first to make sure it is in there).

```
require("tidyr")
source("sync_to_counts_bowtie_3R.R")
#Possible fail; not outright kill?

require("tidyr")
source("sync_to_counts_bowtie_2R.R")
#Killed?

require("tidyr")
source("sync_to_counts_bowtie_3L.R")
#Killed?


require("tidyr")
source("sync_to_counts_bowtie_2L.R")
#Killed

require("tidyr")
source("sync_to_counts_bowtie_4.R")
#Worked: Weird warning - could not allocate memory -- but have a .csv

require("tidyr")
source("sync_to_counts_bowtie_X.R")
# Long run time: Warning at end:
#Error in unlist(x, recursive = FALSE) :
#long vectors not supported yet: memory.c:1648
#No output
```

Rerunning with removing intermediate steps (Change all scripts to remove intermediate data.frames etc. and only have the current form worked on)

```
require("tidyr")
source("sync_to_counts_bowtie_3R.R")

#This far without crashing (using ~36% data)
#[1] "Read in Data"
#[1] "Make long"
#[1] "Remove episodic_counts"




require("tidyr")
source("sync_to_counts_bowtie_2R.R")


require("tidyr")
source("sync_to_counts_bowtie_3L.R")

require("tidyr")
source("sync_to_counts_bowtie_2L.R")

require("tidyr")
source("sync_to_counts_bowtie_X.R")
```



### BWA
Screen 
Start R
```
R
```
and source the file (should require("tidyr") first to make sure it is in there).

source package: do each in unison
```
source("sync_to_counts_3R.R")
source("sync_to_counts_2R.R")
source("sync_to_counts_3L.R")
source("sync_to_counts_2L.R")
source("sync_to_counts_4.R")
source("sync_to_counts_X.R")
```

How to break up the chromosomes into smaller pieces to work with>
Break the 2R for R
```
wc -l episodic_data_2R.sync
```

wc -l episodic_data_2R.sync/ X == Managable sizes

Ex for 2R split every ~20,000:
0,2054752,4109504,6164256,8219008,10273760,12328512,14383264,16438016,18492768,20547525
```
sed -n ' 1, 2054752 p' ../episodic_data_2R.sync > episodic_data_2R_1.sync
sed -n ' 2054753, 4109504 p' ../episodic_data_2R.sync > episodic_data_2R_2.sync
sed -n ' 4109505, 6164256 p' ../episodic_data_2R.sync > episodic_data_2R_3.sync
sed -n ' 6164257, 8219008 p' ../episodic_data_2R.sync > episodic_data_2R_4.sync
sed -n ' 8219009, 10273760 p' ../episodic_data_2R.sync > episodic_data_2R_5.sync
sed -n ' 10273761, 12328512 p' ../episodic_data_2R.sync > episodic_data_2R_6.sync
sed -n ' 12328513, 14383264 p' ../episodic_data_2R.sync > episodic_data_2R_7.sync
sed -n ' 14383265, 16438016 p' ../episodic_data_2R.sync > episodic_data_2R_8.sync
sed -n ' 16438017, 18492768 p' ../episodic_data_2R.sync > episodic_data_2R_9.sync
sed -n ' 18492769, 20547525 p' ../episodic_data_2R.sync > episodic_data_2R_10.sync
```

R script to run model test: episodic_data_bowtie_4.csv
- doing on local to have a look at it;
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/bowtie/R_bowtie/R_scripts/episodic_data_bowtie_4.csv /Users/paulknoops/Bioinformatics/episodic_practice
```
Looking at the data; looks weird and repeated similarites across positions ...

Script to break apart each large file into smaller files:
```
#! /bin/bash

#Variable location
SyncFiles=/home/paul/episodicData/bowtie/R_bowtie

#.sync file to be broken apart:

sync[0]=${SyncFiles}/episodic_data_bowtie_3R.gatk.sync
sync[1]=${SyncFiles}/episodic_data_bowtie_2R.gatk.sync
sync[2]=${SyncFiles}/episodic_data_bowtie_3L.gatk.sync
sync[3]=${SyncFiles}/episodic_data_bowtie_2L.gatk.sync
sync[4]=${SyncFiles}/episodic_data_bowtie_4.gatk.sync
sync[5]=${SyncFiles}/episodic_data_bowtie_X.gatk.sync

for syncro in ${sync[@]}
do
#mkdir ${SyncFiles}/Dir_${syncro}
#syncrodir=${SyncFiles}/Dir_${syncro}

length=($(wc -l ${syncro}))
echo "${length}" >> ${SyncFiles}/wordcounts.txt

done


```
```
#! /bin/bash

#Variable location
SyncFiles=/home/paul/episodicData/bowtie/R_bowtie


sync[0]=${SyncFiles}/episodic_data_bowtie_3R.gatk.sync
#length[0]=27224824

sync[1]=${SyncFiles}/episodic_data_bowtie_2R.gatk.sync
#length[1]=20420127

sync[2]=${SyncFiles}/episodic_data_bowtie_3L.gatk.sync
#length[2]=23647807

sync[3]=${SyncFiles}/episodic_data_bowtie_2L.gatk.sync
#length[3]=22085857

sync[4]=${SyncFiles}/episodic_data_bowtie_X.gatk.sync
#length[4]=21432322

#Removing b/c complete and short
#sync[5]=${SyncFiles}/episodic_data_bowtie_4.gatk.sync
#length[5]=1235697


for file in ${sync[@]}
do
name=${file}
base=`basename ${name} .gatk.sync`

mkdir ${SyncFiles}/${base}_dir
basedir=${SyncFiles}/${base}_dir

length=($(wc -l ${SyncFiles}/${base}.gatk.sync))

sed -n ' 1, 2054752 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_1.sync
sed -n ' 2054753, 4109504 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_2.sync
sed -n ' 4109505, 6164256 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_3.sync
sed -n ' 6164257, 8219008 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_4.sync
sed -n ' 8219009, 10273760 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_5.sync
sed -n ' 10273761, 12328512 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_6.sync
sed -n ' 12328513, 14383264 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_7.sync
sed -n ' 14383265, 16438016 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_8.sync
sed -n ' 16438017, 18492768 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_9.sync
sed -n " 18492769, ${length} p" ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_10.sync

done

```
would need many R scripts to run: possibly change the constant change of variables
probably made to many: each b/w 300-600

Using a variation of sync to counts; creating 10 .csv files for analysis with 2 for loops in R
first
```
mkdir subsettingDirectories
```
Then from file with all the directories made above 
```
all the ${base}_dir in ${SyncFiles}/
```
```
mv *_dir subsettingDirectories/
```
Then run the new R script
```
source("R_loop_sync_to_counts.R")
```

which can be seen with BWA below (which is edited for BWA

For this run: most of 3R not done (the _2 was killed).
	- means the X failed as well
	- move to own dir in that one and rerun script with changed dir call.

Will need to rename .csv files; will have 50 .csv files to combine later as well....

should have removed all intermediate things at the end again.. but oh well


BWA:
Check the lenths
```
#! /bin/bash

#Variable location
SyncFiles=/home/paul/episodicData/R_dir

#.sync file to be broken apart:

sync[0]=${SyncFiles}/episodic_data_3R.gatk.sync
sync[1]=${SyncFiles}/episodic_data_2R.gatk.sync
sync[2]=${SyncFiles}/episodic_data_3L.gatk.sync
sync[3]=${SyncFiles}/episodic_data_2L.gatk.sync
sync[4]=${SyncFiles}/episodic_data_X.gatk.sync
sync[5]=${SyncFiles}/episodic_data_4.gatk.sync

for syncro in ${sync[@]}
do
length=($(wc -l ${syncro}))
echo "${length}" >> ${SyncFiles}/wordcounts.txt
done
```
27293799
20547525
23769939
22169930
21629162
1277434
Will be one very large one and a few small ones (2R == so very small..
```
#! /bin/bash

#Variable location
SyncFiles=/home/paul/episodicData/R_dir

#output dir:
subsets=/home/paul/episodicData/R_dir/bwa_subsetDirectories

sync[0]=${SyncFiles}/episodic_data_3R.gatk.sync
sync[1]=${SyncFiles}/episodic_data_2R.gatk.sync
sync[2]=${SyncFiles}/episodic_data_3L.gatk.sync
sync[3]=${SyncFiles}/episodic_data_2L.gatk.sync
sync[4]=${SyncFiles}/episodic_data_X.gatk.sync

#Removing b/c complete and short
#sync[5]=${SyncFiles}/episodic_data_4.gatk.sync
# Run on its own

for file in ${sync[@]}
do
name=${file}
base=`basename ${name} .gatk.sync`

mkdir ${subsets}/${base}_dir
basedir=${subsets}/${base}_dir

length=($(wc -l ${SyncFiles}/${base}.gatk.sync))

sed -n ' 1, 2054752 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_1.sync
sed -n ' 2054753, 4109504 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_2.sync
sed -n ' 4109505, 6164256 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_3.sync
sed -n ' 6164257, 8219008 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_4.sync
sed -n ' 8219009, 10273760 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_5.sync
sed -n ' 10273761, 12328512 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_6.sync
sed -n ' 12328513, 14383264 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_7.sync
sed -n ' 14383265, 16438016 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_8.sync
sed -n ' 16438017, 18492768 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_9.sync
sed -n ' 18492769, 20547500 p' ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_10.sync
sed -n " 20547501, ${length} p" ${SyncFiles}/${base}.gatk.sync > ${basedir}/${base}_11.sync
done
```

copy the 4th chromo into subset directory

```
mkdir episodic_data_4_dir
cp episodic_data_4.gatk.sync episodic_data_4_dir
```
```
#For loop sync to counts

# To run: open R (> R) and source this script (be sure to edit based on file used). 

## Convert a .sync file into long format, filter somewhat, and have only position, treatment, Cage, Generation and Maj/Min counts

#Source a packages script with necessary packages needed
#source("packages.R")

### Packages source code: only need these two for this script
require('dplyr')

#tidyr may not work: require on own before running
require('tidyr')


#1) Need to change details as needed above and below string of #####

#2) Needs a .sync file made by popoolation2

#3) Need to change most importantly for analysis the read in and read out names 

# Read in Data: Big Data Sets
#Pwd a direcotry containing only the directories of interest (made with other sed -n script)

mydirs <- list.dirs(path = "/home/paul/episodicData/R_dir/bwa_subsetDirectories", recursive = FALSE)

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

Running for BOWTIE_4.csv a model test with script and .csv of coefficients

This version will include the just run 4th chromosome, maybe change that up?
```
#Chromosome 4 practice:
#Episodic data analysis:


episodic_long <- read.csv("/home/paul/episodicData/bowtie/R_bowtie/R_scripts/episodic_data_bowtie_4.csv", h=TRUE)

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
  #Don't need to print this out necessarily, just nice to track
  print(paste("Running entity:", i, "which is", which(position==i), "out of", no.pos))
  
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
  
  #Add to data frame (total for the whole data set == coeffs_df + the newly made X)
  coeffs_df <- rbind(coeffs_df, x)
  
  #Remove i for safety and it starts over
  rm(i)
}

#Change column names to workable
colnames(coeffs_df) <- c("Estimate", "Standard_error", "z-value", "p-value", "position")

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

#coeffs_df <- coeffs_df[grep("Treatment", coeffs_df$Effects), ]
#coeffs_df <- coeffs_df[which(coeffs_df$Effects=="TreatmentSel"),]
#Save data frame???
write.csv(coeffs_df, file="episodic_bowtie_4_coeffs.csv", row.names = FALSE)

```





Worked; need to make a loop for each directory: change working directory?

```
#Episodic data analysis: loop .csv files to run model:

#change to directory holding all directories:

#mydirs <- list.dirs(path = "/home/paul/episodicData/R_dir/bwa_subsetDirectories", recursive = FALSE)


#for (dir in mydirs){
  
#setwd('/home/paul/episodicData/bowtie/R_bowtie/bowtie_subset_Ranalysis/')

# make script into '/home/paul/episodicData/bowtie/R_bowtie/bowtie_subset_Ranalysis/' and open from same direcotry screen/R
#setwd("episodic_data_bowtie_2R_dir")
#setwd("episodic_data_bowtie_3L_dir")
setwd("episodic_data_bowtie_3R_dir")
#setwd("episodic_data_bowtie_4_dir")
#setwd("episodic_data_bowtie_X_dir")
  
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
#}

```
Running Bowtie is a large loop (should stop): moved the other directories (Not 2L which is first) to own directory and making own loop (above)
	-- "Running entity: 20598073 which is 59410 out of 343990 file= episodic_data_bowtie_3R_10.sync.csv"
	-- Stuck here; .... used ctrl X and it restarted....?

running BWA per chromosome in a loop
	-- FWI: 2R for bwa number 11 was so small and nothing was left after filtering (so removed)

name screens:
```
screen -S NAME
```

the 4th chromsome finished first: bring to local to practice:
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/bowtie/R_bowtie/bowtie_subset_Ranalysis/episodic_data_bowtie_4_dir/episodic_data_bowtie_4.csv.coeffs.csv /Users/paulknoops/Bioinformatics/episodic_practice

scp paul@info.mcmaster.ca:/home/paul/episodicData/R_dir/bwa_subsetDirectories/episodic_data_4_dir/episodic_data_4.gatk.sync.csv.coeffs.csv /Users/paulknoops/Bioinformatics/episodic_practice

```


Ancestral Pi (run from new Ancestral directory)
```
#Make pileup first for only MGD3
samtools mpileup MGD3_SO_CAGATC_merged_aligned_pe.final_realigned.bam  > MGD3.pileup
```
run variance sliding
```
perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input MGD3.pileup --output MGD3.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120
```

Visualize:
```
perl /home/paul/popoolation_1.2.2/Visualise-output.pl --input MGD3.pi --output MGD3.pi.pdf --ylab pi --chromosomes "X 2L 2R 3L 3R 4"
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/Ancestor/MGD3.pi.pdf /Users/paulknoops/Bioinformatics/episodic_practice
scp paul@info.mcmaster.ca:/home/paul/episodicData/Ancestor/MGD3.pi /Users/paulknoops/Bioinformatics/episodic_practice
```

Bowtie: full script:
```
#! /bin/bash

input=/home/paul/episodicData/bowtie/MDG3_bow
samtools mpileup ${input}/MGD3_SO_CAGATC_merged_bowtie_pe.final_realigned.bam > ${input}/MGD3_bowtie.pileup

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${input}/MGD3_bowtie.pileup --output ${input}/MGD3_bowtie.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

perl /home/paul/popoolation_1.2.2/Visualise-output.pl --input ${input}/MGD3_bowtie.pi --output ${input}/MGD3_bowtie.pi.pdf --ylab pi --chromosomes "X 2L 2R 3L 3R 4"
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/bowtie/MDG3_bow/MGD3_bowtie.pi.pdf /Users/paulknoops/Bioinformatics/episodic_practice
scp paul@info.mcmaster.ca:/home/paul/episodicData/bowtie/MDG3_bow/MGD3_bowtie.pi /Users/paulknoops/Bioinformatics/episodic_practice
```

```

scp paul@info.mcmaster.ca:/home/paul/Chromosomes/*.csv /Users/paulknoops/Bioinformatics/episodic_practice/DATA

scp paul@info.mcmaster.ca:/home/paul/Chromosomes/*plots.R /Users/paulknoops/Bioinformatics/episodic_practice/DATA
```


```
#! /bin/bash

input=/home/paul/episodicData/gatk_dir
output=/home/paul/episodicData/115
#samtools mpileup ${input}/F115ConR1_TAGCTT_merged_aligned_pe.final_realigned.bam > ${output}/F115ConR1bwa.pileup
#samtools mpileup ${input}/F115ConR2_GGCTAC_merged_aligned_pe.final_realigned.bam > ${output}/F115ConR2bwa.pileup
#samtools mpileup ${input}/F115SelR1_GTTTCG_merged_aligned_pe.final_realigned.bam > ${output}/F115SelR1bwa.pileup
samtools mpileup ${input}/F115SelR2_GTGGCC_merged_aligned_pe.final_realigned.bam > ${output}/F115SelR2bwa.pileup

#perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115ConR1bwa.pileup --output ${output}/F115ConR1bwa.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

#perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115ConR2bwa.pileup --output ${output}/F115ConR2bwa.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

#perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115SelR1bwa.pileup --output ${output}/F115SelR1bwa.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115SelR2bwa.pileup --output ${output}/F115SelR2bwa.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120
```


```
scp paul@info.mcmaster.ca:/home/paul/episodicData/115/*.pi /Users/paulknoops/Bioinformatics/episodic_practice/F115
```

ALL Failed: Try Again (change coverage/quality?)
Try Bowtie:

```
#! /bin/bash

input=/home/paul/episodicData/bowtie/gatk_bowtie
output=/home/paul/episodicData/bowtie/115_bow

#samtools mpileup ${input}/F115ConR1_TAGCTT_merged_bowtie_pe.final_realigned.bam > ${output}/F115ConR1bowtie.pileup
#samtools mpileup ${input}/F115ConR2_GGCTAC_merged_bowtie_pe.final_realigned.bam > ${output}/F115ConR2bowtie.pileup
#samtools mpileup ${input}/F115SelR1_GTTTCG_merged_bowtie_pe.final_realigned.bam > ${output}/F115SelR1bowtie.pileup
samtools mpileup ${input}/F115SelR2_GTGGCC_merged_bowtie_pe.final_realigned.bam > ${output}/F115SelR2bowtie.pileup

#perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115ConR1bowtie.pileup --output ${output}/F115ConR1bowtie.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

#perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115ConR2bowtie.pileup --output ${output}/F115ConR2bowtie.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

#perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115SelR1bowtie.pileup --output ${output}/F115SelR1bowtie.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${output}/F115SelR2bowtie.pileup --output ${output}/F115SelR2bowtie.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120
```
Check Encoding: could be sanger again issue????

Lots of NA's ?? mpileup issue, or coverage issue? or what....

```

#! /bin/bash

input=/home/paul/episodicData/bowtie/gatk_bowtie
output=/home/paul/episodicData/bowtie/bowtie_pileups
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

files=(${input}/*_pe.final_realigned.bam)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _pe.final_realigned.bam`

samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_pe.final_realigned.bam > ${output}/${base}.pileup

done
```

```

#! /bin/bash

input=/home/paul/episodicData/bowtie/bowtie_pileups
output=/home/paul/episodicData/bowtie/bowtie_pi

files=(${input}/*_merged_bowtie.pileup)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _merged_bowtie.pileup`

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input ${input}/${base}_merged_bowtie.pileup --output ${output}/${base}.bowtie.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120

done
```
Still na's throughout:

Issue advice from https://sourceforge.net/p/popoolation/wiki/Manual/ comments:

	"Values of 0 and na in my case came from using a bam file generated from reads using phred33 scores when the default of the Variance-sliding.pl script is to use phred64 encoding. Specify "--fastq-type sanger" to use the phred33 scoring scheme."
	
Sanger Test:
--fastq-type sanger
```
perl /home/paul/popoolation_1.2.2/Variance-sliding.pl \
	--input /home/paul/episodicData/bowtie/bowtie_pileups/F115ConR1_TAGCTT_merged_bowtie.pileup \
	--output /home/paul/episodicData/bowtie/bowtie_pi/Test_F115ConR1_TAGCTT_merged_bowtie.pi \
	--measure pi \
	--window-size 10000 \
	--step-size 10000 \
	--min-count 2 \
	--min-coverage 4 \
	--max-coverage 400 \
	--min-qual 20 \
	--pool-size 120 \
	--fastq-type sanger
```
WORKED!

```
#! /bin/bash

input=/home/paul/episodicData/bowtie/bowtie_pileups
output=/home/paul/episodicData/bowtie/bowtie_pi

files=(${input}/*_merged_bowtie.pileup)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _merged_bowtie.pileup`

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl \
	--input ${input}/${base}_merged_bowtie.pileup \
	--output ${output}/${base}.bowtie.pi \
	--measure pi \
	--window-size 10000 \
	--step-size 10000 \
	--min-count 2 \
	--min-coverage 4 \
	--max-coverage 400 \
	--min-qual 20 \
	--pool-size 120 \
	--fastq-type sanger

done

```

Move it to personal:
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/bowtie/bowtie_pi/*.pi /Users/paulknoops/Bioinformatics/episodic_practice/Bowtie_Pi
```

Read each seperate into function below to create plot ( R function )
```
Pi_PlotFunction <- function(x) {
  require(ggplot2)
  x2 <- gsub("\\_.*","",x)
#Read in the data:
  Datt <- read.table(x)
  colnames(Datt) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  #Remove unnecessary regions: Not necessary based on later steps
  Datt$chr <- as.character(Datt$chr)
  Datt2 <- Datt
  
  #Datt2 <- Datt[-which(Datt$chr=="YHet"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="2RHet"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="2LHet"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="3LHet"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="3RHet"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="U"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="XHet"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="dmel_mitochondrion_genome"),]
  #Datt2 <- Datt2[-which(Datt2$chr=="Uextra"),]
  
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
  
  Pi_plot <- ggplot(DattFull, aes(x = number, y= Pi, colour = chr)) 
  
  Pi_plot_2 <- Pi_plot + 
    geom_point(size=0.3, show.legend = F) +
    scale_y_continuous(limits=c(0, 0.02), breaks=seq(0, 0.02, 0.005)) + 
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    theme(text = element_text(size=20), 
          axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) +
    ggtitle(x2)
  
  return(Pi_plot_2)
}

```


______________________________________________________________________________________________________________________________________

NOTE: Running in parallel (in background) and screen is closed: list the jobs running still with pgrep PROGRAM (ex. pgrep novoalign)
https://www.digitalocean.com/community/tutorials/how-to-use-bash-s-job-control-to-manage-foreground-and-background-processes


Run Trimmomatic with different MAXINFO 0.3 -- 0.8

For Variant callers (or MPILEUP), don't use the reference sequence, but use the ancestor (turn into a .fasta.... file?????)

Present on the analysis in two weeks!
	-- Present Tested Trim Files
	-- Present CRISP or other variant callers?
	-- Thomas Taus Scripts Testing????
	
	
______________

Running Novoalign full data set: error

Remidied here (Merged to early)

Script: novo_MGD_merge.sh
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_merge=${project_dir}/novo_merge

samtools merge ${novo_merge}/MGD3_SO_CAGATC_novo_merge.bam \
${novo_merge}/MGD2_SO_CAGATC_novo_merge.bam \
${novo_merge}/MGD_SO_CAGATC_novo_merge.bam

```

Need to move the ancestor unmerged away (so not read for later steps)
```
mkdir ancestorUnmerged
mv MGD2_SO_CAGATC_merged_aligned_pe.final.bam ancestorUnmerged
mv MGD_SO_CAGATC_merged_aligned_pe.final.bam ancestorUnmerged
```

Note: the process fails for the ancestor b/c the merging early disrupted the ID's for reading with picard: rerun here (in one script

need the ancsetor Unmerged dir and a new picard directory for ancestor to be put into (mkdir pic_ancestor)

```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_merge=${project_dir}/novo_merge/ancestorUnmerged

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to output directory
novo_pic=${project_dir}/novo_pic/pic_ancestor

#Path to tmp
novo_tmp=${project_dir}/novo_tmp

#Path to output directory
novo_rmd=${project_dir}/novo_rmd


files=(${novo_merge}/*.bam)
for file in ${files[@]}
do
name=${file}

base=`basename ${name} .bam`

java -Xmx2g -Djava.io.tmpdir=${novo_tmp} -jar ${pic} SortSam \
I= ${novo_merge}/${base}.bam \
O= ${novo_pic}/${base}_novo_sort.bam \
VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${novo_tmp}

java -Xmx2g -jar ${pic} MarkDuplicates \
I= ${novo_pic}/${base}_novo_sort.bam \
O= ${novo_rmd}/${base}_novo_rmd.bam \
M= ${novo_rmd}/dupstat_anc.txt \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES= true

done
```






__Because I merged the ancestor early (before Picard), the markduplicates fails, maybe try picard merge__

__PICARD MERGING:__
```
java -jar picard.jar MergeSamFiles \
      I=input_1.bam \
      I=input_2.bam \
      O=merged_files.bam
```
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_pic_merge=${project_dir}/novo_pic_merge

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

files=(${novo_bam}/*_L001_novo.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_novo.bam`
java -jar ${pic} MergeSamFiles I=${novo_bam}/${base}_L001_novo.bam I=${novo_bam}/${base}_L002_novo.bam O=${novo_pic_merge}/${base}_novo_pic_merge.bam
done
```

```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_pic_merge=${project_dir}/novo_pic_merge

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

java -jar ${pic} MergeSamFiles \
	I=${novo_pic_merge}/MGD2_SO_CAGATC_novo_pic_merge.bam\
	I=${novo_pic_merge}/MGD_SO_CAGATC_novo_pic_merge.bam \
	O=${novo_pic_merge}/MGD3_SO_CAGATC_novo_merge.bam
	
mkdir ${novo_pic_merge}/Anc_unmerged
mv ${novo_pic_merge}/MGD2_SO_CAGATC_novo_pic_merge.bam ${novo_pic_merge}/Anc_unmerged
mv ${novo_pic_merge}/MGD_SO_CAGATC_novo_pic_merge.bam ${novo_pic_merge}/Anc_unmerged
```

Comparing Samtools merge and Picard Merge: 

	-- Picard Merge has slightly larger file sizes (better at merging???)
	-- Potentially more b/c Read Group (RG) information is preserved??
	-- Attempt the rest of the pipeline in duplicate:
	

```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_pic_merge=${project_dir}/novo_pic_merge

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to output directory
novo_pic=${project_dir}/novo_picMerge_sort

#Path to tmp
novo_tmp=${project_dir}/novo_tmp

#Path to rmd directory
novo_rmd=${project_dir}/novo_pic_rmd

files=(${novo_pic_merge}/*.bam)
for file in ${files[@]}
do
name=${file}

base=`basename ${name} .bam`

java -Xmx2g -Djava.io.tmpdir=${novo_tmp} -jar ${pic} SortSam \
I= ${novo_pic_merge}/${base}.bam \
O= ${novo_pic}/${base}_novo_sort.bam \
VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${novo_tmp}

java -Xmx2g -jar ${pic} MarkDuplicates \
I= ${novo_pic}/${base}_novo_sort.bam \
O= ${novo_rmd}/${base}_novo_rmd.bam \
M= ${novo_rmd}/dupstat_anc.txt \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES= true

done
```
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_rmd=${project_dir}/novo_pic_rmd

#Path to output directory
novo_final=${project_dir}/novo_pic_final


files=(${novo_rmd}/*_novo_rmd.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _novo_rmd.bam`
samtools view -q 20 -F 0x0004 -b ${novo_rmd}/${base}_novo_rmd.bam > ${novo_final}/${base}_novo_final.bam
done
```

```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to final directory
novo_final=${project_dir}/novo_pic_final

samtools merge ${novo_final}/MGD3_SO_CAGATC_novo_pic_merge_novo_final.bam \
${novo_final}/MGD2_SO_CAGATC_novo_pic_merge_novo_final.bam \
${novo_final}/MGD_SO_CAGATC_novo_pic_merge_novo_final.bam

mkdir ${novo_final}/Anc_unmerged
mv ${novo_final}/MGD2_SO_CAGATC_novo_pic_merge_novo_final.bam ${novo_final}/Anc_unmerged
mv ${novo_final}/MGD_SO_CAGATC_novo_pic_merge_novo_final.bam ${novo_final}/Anc_unmerged
```
All the final files look the same: this is a possible alternative method, but does not really matter it looks like....

### Subsample.sync script attempts

__The first attempt done far above__

Can also be done on individual pileup files

Using the subsample synchronized script help::::
```
perl subsample_synchronized.pl --help
```
```
SUBSAMPLE-SYNCHRONIZEDUser Contributed Perl DocumentaSUBSAMPLE-SYNCHRONIZED(1)

NAME
       perl subsample-synchronized.pl - Reduce the coverage of a synchronized file, by random subsampling, to the given target coverage

SYNOPSIS
       perl subsample-synchronized.pl --input input.sync --output output.sync --target-coverage 50 --max-coverage 2%  --method
       withoutreplacement

OPTIONS
       --input
           The input file in the synchronized format; Mandatory.

       --output
           The output file, will be a synchronized file  Mandatory.

       --target-coverage
           Reduce the coverage of the pileup-file to the here provided value; The target coverage also acts as minimum coverage, i.e.: if
           the coverage in any population is smaller than the targetcoverage the whole pileup entry is discarded. Mandatory

       --max-coverage
           The maximum coverage; All populations are required to have coverages lower or equal than the maximum coverage; Mandatory The
           maximum coverage may be provided as one of the following:
              '500' a maximum coverage of 500 will be used for all populations
              '300,400,500' a maximum coverage of 300 will be used for the first population, a maximum coverage of 400 for the second population and so on
            '2%' the 2% highest coverages will be ignored, this value is independently estimated for every population

       --method
           Specify the method for subsampling of the synchronized file. Either: withreplace, withoutreplace, fraction; Mandatory

            withreplace: subsample with replacement
            withoutreplace: subsample without replacement
            fraction: calculate the exact fraction of the allele frequencies and linearly scale them to the C<--target-coverage> and rounding to an integer;

       --help
           Display help for this script

Details

 Input synchronized
       A synchronized file

   Output
       The output will be a reduced coverage synchronized file
```

Info on fastqc
```
https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/
```

### Steps after pileups made (testing on Novoalign files)
 First run of Tajima Pi for Novoalign: some files are empty? Not sure why?
```
#!/bin/bash

samtools mpileup -B -Q 0 -f /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta /home/paul/episodicData/novoalign/novo_GATK/F115SelR2_GTGGCC_novo_merge_novo_final_realigned.bam > /home/paul/episodicData/novoalign/novo_pileup/F115SelR2_GTGGCC_novo.pileup

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input /home/paul/episodicData/novoalign/novo_pileup/F115SelR2_GTGGCC_novo.pileup --output /home/paul/episodicData/novoalign/novo_pileup/F115SelR2_GTGGCC_novo.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output /home/paul/episodicData/novoalign/novo_pileup/F115SelR2_GTGGCC_novo.snps --min-covered-fraction 0.5
```
```
#!/bin/bash

samtools mpileup -B -Q 0 -f /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta /home/paul/episodicData/novoalign/novo_GATK/F115SelR1_GTTTCG_novo_merge_novo_final_realigned.bam > /home/paul/episodicData/novoalign/novo_pileup/F115SelR1_GTTTCG_novo.pileup

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input /home/paul/episodicData/novoalign/novo_pileup/F115SelR1_GTTTCG_novo.pileup --output /home/paul/episodicData/novoalign/novo_pileup/F115SelR1_GTTTCG_novo.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output /home/paul/episodicData/novoalign/novo_pileup/F115SelR1_GTTTCG_novo.snps --min-covered-fraction 0.5
```

```
#!/bin/bash

samtools mpileup -B -Q 0 -f /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta /home/paul/episodicData/novoalign/novo_GATK/F38SelR1_ACTTGA_novo_merge_novo_final_realigned.bam > /home/paul/episodicData/novoalign/novo_pileup/F38SelR1_ACTTGA_novo.pileup

perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input /home/paul/episodicData/novoalign/novo_pileup/F38SelR1_ACTTGA_novo.pileup --output /home/paul/episodicData/novoalign/novo_pileup/F38SelR1_ACTTGA_novo.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output /home/paul/episodicData/novoalign/novo_pileup/F38SelR1_ACTTGA_novo.snps --min-covered-fraction 0.5

```
F38SelR1_ACTTGA_novo.pi

```
perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input /home/paul/episodicData/novoalign/novo_pileup/F38ConR1_ATCACG_novo.pileup --output /home/paul/episodicData/novoalign/novo_pileup/F38ConR1_ATCACG_novo.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output /home/paul/episodicData/novoalign/novo_pileup/F38ConR1_ATCACG_novo.snp --min-covered-fraction 0.5
```

```
perl /home/paul/popoolation_1.2.2/Variance-sliding.pl --input /home/paul/episodicData/novoalign/novo_pileup/F38ConR2_TTAGGC_novo.pileup --output /home/paul/episodicData/novoalign/novo_pileup/F38ConR2_TTAGGC_novo.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output /home/paul/episodicData/novoalign/novo_pileup/F38ConR2_TTAGGC_novo.snp --min-covered-fraction 0.5
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_pileup/*.pi /Users/paulknoops/Bioinformatics/R-projects_git/episodicSequenceData/R_scripts/Pi_Analysis_Novo
```
Then move them back over to novo_pi

1) Could subsample to uniform coverage (but why for Popoolation1 looking within one population)
	--> from slides _"Several population genetic estimators are sensitive to sequencing errors. For example a very low Tajima’s D, usually indicative of a selective sweep, may be, as an artifact, frequently be found in highly covered regions because these regions have just more sequencing errors. To avoid these kinds of biases we recommend to subsample to an uniform coverage."_

2) Filter Indels -- 
	--> From Slide: _"FILTERING INDELS
	perl ˜/programs/popoolation/basic-pipeline/identifygenomic-indel-regions.pl --indel-window 5 --mincount --input pe.mpileup --output indels.gtf
	I –indel-window how many bases surrounding indels should be ignored
	I –min-count minimum count for calling an indel. Note that indels may be sequencing errors as well
	
	_perl ˜/programs/popoolation/basic-pipeline/filterpileup-by-gtf.pl --input pe.mpileup --gtf indels.gtf --output pe.idf.mpileup

	_Note: the filter-pileup script could also be used to remove entries overlapping with transposable elements (RepeatMasker produces a gtf as well)."_

	-- Necessary if IndelRealignment Done???

3) Now Run analysis of Tajima's Pi


4) Gene by Gene Tajima's Pi

- Download most recent .gff file --> Go to the FlyBase homepage (http://flybase.org/) and get the annotation: dmel-all-r6.18.gff.gz
	 (Should get the matching to current index (that was mapped with) - r5.57

- Move to Brians Machine
```
scp /Users/paulknoops/Downloads/dmel-all-r5.57.gff.gz paul@info.mcmaster.ca:/home/paul/episodicData/index_dir
```
- Unzip the file.
```
gunzip dmel-all-r5.57.gff.gz
```
.... No space left on device....

- Filter for exons and convert it into a gtf file:
```
cat /home/paul/episodicData/index_dir/dmel-all-r5.57.gff| awk '$2=="FlyBase" && $3=="exon"'| perl -pe 's/ID=([^:;]+)([^;]+)?;.*/gene_id "$1"; transcript_id "$1:1";/'> /home/paul/episodicData/novoalign/novo_exons/exons.gtf
```
From Popoolation:

cat dmel-all-r5.22.gff| awk '$2=="FlyBase" && $3=="exon"'| perl -pe 's/ID=([^:;]+)([^;]+)?;.*/gene_id "$1"; transcript_id "$1:1";/'> exons.gtf

Does not work properly...not a properly formated gtf...

TO check whats going on ....
```
cat /home/paul/episodicData/index_dir/dmel-all-r5.57.gff| awk '$2=="FlyBase" && $3=="exon"'> exons.gff
```
Looks the same as exons.gtf....


Trying cufflinks (spoiler.. does not work for this..)
```
cat dmel-all-r5.57.gff| awk '$2=="FlyBase" && $3=="exon"'| /usr/local/cufflinks1.1.0/gffread -T -o my.gtf
```
/usr/local/cufflinks1.1.0/gffread dmel-all-r5.57.gff -T -o my.gtf

gffread my.gff3 -T -o my.gtf

Issue seems to be with the /ID from popoolation script:

```
cat /home/paul/episodicData/index_dir/dmel-all-r5.57.gff| awk '$2=="FlyBase" && $3=="exon"'| perl -pe 's/Name=([^:;]+)([^;]+)?;.*/gene_id "$1"; transcript_id "$1:1";/'> /home/paul/episodicData/novoalign/novo_exons/exons.gtf
```
THIS WORKS!



- run Variance at position:

```
perl /home/paul/popoolation_1.2.2/Variance-at-position.pl --pool-size 120 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 400 --pileup /home/paul/episodicData/novoalign/novo_pileup/MGD3_SO_CAGATC_novo.pileup --gtf /home/paul/episodicData/novoalign/novo_exons/exons.gtf --output /home/paul/episodicData/novoalign/novo_exons/MGD3_SO_CAGATC_novo.genes.pi --measure pi --fastq-type sanger
```
Stopped early?
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_exons/MGD3_SO_CAGATC_novo.genes.pi /Users/paulknoops/Bioinformatics/episodic_practice/MGD3
```


```
perl /home/paul/popoolation_1.2.2/Variance-at-position.pl --pool-size 120 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 400 --pileup /home/paul/episodicData/novoalign/novo_pileup/MGD3_SO_CAGATC_novo.pileup --gtf /home/paul/episodicData/novoalign/novo_exons/2R-exons.gtf --output /home/paul/episodicData/novoalign/novo_exons/MGD3_SO_CAGATC_novo.2Rgenes.pi --measure pi --fastq-type sanger
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_exons/MGD3_SO_CAGATC_novo.2Rgenes.pi /Users/paulknoops/Bioinformatics/episodic_practice/MGD3
```
Compare to a F115 Selection: NEED TO RUN STILL!

```
perl /home/paul/popoolation_1.2.2/Variance-at-position.pl --pool-size 120 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 400 --pileup /home/paul/episodicData/novoalign/novo_pileup/F115SelR1_GTTTCG_novo.pileup --gtf /home/paul/episodicData/novoalign/novo_exons/2R-exons.gtf --output /home/paul/episodicData/novoalign/novo_exons/F115SelR1_GTTTCG_novo.2Rgenes.pi --measure pi --fastq-type sanger
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_exons/F115SelR1_GTTTCG_novo.2Rgenes.pi /Users/paulknoops/Bioinformatics/episodic_practice/MGD3
```

### 2R Tajima's Pi script:

```
#! /bin/bash

# Path to PoPoolation1 (Currently in Paul's Home directory)
popoolation=/home/paul/popoolation_1.2.2

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_pileup

# Path to output Tajima Pi files
output=${project_dir}/novo_exons/2R_dir

#gtf file
gtf=${project_dir}/novo_exons


files=(${input}/*.pileup)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} .pileup`
perl ${popoolation}/Variance-at-position.pl \
	--pool-size 120 \
	--min-qual 20 \
	--min-coverage 4 \
	--min-count 2 \
	--max-coverage 400 \
	--pileup ${input}/${base}.pileup \
	--gtf ${gtf}/2R-exons.gtf \
	--output ${output}/${base}.2Rgenes.pi \
	--measure pi \
	--fastq-type sanger

done
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_exons/2R_dir/*.pi /Users/paulknoops/Bioinformatics/episodic_practice/2R_GenePi
```

```
middle=$((`wc -l < file` / 2))
```

#Done but Dumb: Forget it \/ \/  \/
```
perl /home/paul/popoolation_1.2.2/Visualise-output.pl --input /home/paul/episodicData/novoalign/novo_pi/MGD3_SO_CAGATC_novo.pi --output /home/paul/episodicData/novoalign/novo_pi/MGD3_SO_CAGATC_novo.pi.pdf --ylab pi --chromosomes "X 2L 2R 3L 3R 4"
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_pi/MGD3_SO_CAGATC_novo.pi.pdf /Users/paulknoops/Bioinformatics/episodic_practice
```
#Done but Dumb: Forget it /\  /\  /\

Rscript /home/paul/episodicData/novoalign/novo_mpileup/RSCriptingTest.R '/home/paul/episodicData/novoalign/novo_mpileup'
```
args <- commandArgs(trailingOnly = TRUE)
# rnorm(n=as.numeric(args[1]), mean=as.numeric(args[2]))
# Rscript myScript.R 5 100

require('tidyr')
require('dplyr')

mydirs <- list.dirs(path = args[1], recursive = FALSE)

for (dir in mydirs){

    setwd(dir)
  
  mysyncs <- list.files(pattern=".sync")
  
  for (sync in mysyncs){
    
    print(sync)
x3 <- gsub("\\..*","",sync)
J3 <- gsub('(.*)_\\w+', '\\1', x3)

X <- args[1]

file=paste(X,"/",J3,"_chromo.csv", sep="")  
print(file)  
  }
}

messy <- data.frame(
  name = c("Wilbur", "Petunia", "Gregory"),
  a = c(67, 80, 64),
  b = c(56, 90, 50)
)
print(messy)

messy2 <- messy %>%
  gather(drug, heartrate, a:b)
  
  print(messy2)
 
```
WORKS!
Issue with tidyr: 
-works with interactive (Running R) -- 3.2.2
-Does not work with Rscript:
	-- different order

```
#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup

Rscript RSCriptingTest.R ${SyncFiles}

```
```
#! /bin/bash
length=20547500

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
echo '$((${cut_6}+1)))'
echo '${cut_10}
```

```
#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup

#Output dir:
subsets=${project_dir}/ChromoSubsets

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts

sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 
sync[5]=${SyncFiles}/novo_episodic_4.sync 

for file in ${sync[@]}
	do
	(name=${file}
	base=`basename ${name} .sync`
	basedir=${subsets}/${base}_dir
	Rscript ${Rscripts}/Counts_to_model_2.R ${basedir}) &	
done
wait

```

Copy of model script that goes through directories (Chromos) in series:
```
### Script: Counts_to_model.R
#	- Script for running each chromosome in unison
#To run on own:

# Rscript Counts_to_model.R 'DIRECTORY'

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

#mkdir ${project_dir}/novo_coeffs
coeff_dir=${project_dir}/novo_coeffs


Rscript ${Rscripts}/Combine_chromo.R ${subsets} ${coeff_dir}
```

Redone for BWA and Bowtie
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

subsets=/home/paul/episodicData/R_dir/bwa_subsetDirectories

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts

mkdir /home/paul/episodicData/R_dir/bwa_coeffs
coeff_dir=/home/paul/episodicData/R_dir/bwa_coeffs


Rscript ${Rscripts}/Combine_chromo.R ${subsets} ${coeff_dir}
```
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

subsets=/home/paul/episodicData/bowtie/R_bowtie/bowtie_subset_Ranalysis

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts

mkdir /home/paul/episodicData/bowtie/R_bowtie/bowtie_coeffs
coeff_dir=/home/paul/episodicData/bowtie/R_bowtie/bowtie_coeffs

Rscript ${Rscripts}/Combine_chromo.R ${subsets} ${coeff_dir}
```


### Subsampling:
withreplacement
```
perl /usr/local/popoolation/subsample-synchronized.pl --input novo_episodic_main.sync --output novo_episodic_main_subsample.sync --target-coverage 40 --max-coverage 200 --method withreplace

```
```
perl /usr/local/popoolation/subsample-synchronized.pl --input novo_episodic_main.sync --output novo_episodic_sub_w.o_replace.sync --target-coverage 40 --max-coverage 200 --method withoutreplace

```

Changing .sync file for selection vs. controls and duplicate ancestor! (ALL FOR TAUS ET AL. 2017== PoolSeq package)
```
F115ConR1_TAGCTT_novo_merge_novo_final_realigned.bam
F115ConR2_GGCTAC_novo_merge_novo_final_realigned.bam
F115SelR1_GTTTCG_novo_merge_novo_final_realigned.bam
F115SelR2_GTGGCC_novo_merge_novo_final_realigned.bam
F38ConR1_ATCACG_novo_merge_novo_final_realigned.bam
F38ConR2_TTAGGC_novo_merge_novo_final_realigned.bam
F38SelR1_ACTTGA_novo_merge_novo_final_realigned.bam
F38SelR2_GATCAG_novo_merge_novo_final_realigned.bam
F77ConR1_ATGTCA_novo_merge_novo_final_realigned.bam
F77ConR2_ATTCCT_novo_merge_novo_final_realigned.bam
F77SelR1_TTAGGC_novo_merge_novo_final_realigned.bam
F77SelR2_GATCAG_novo_merge_novo_final_realigned.bam
MGD3_SO_CAGATC_novo_merge_novo_final_realigned.bam
```
FOR SELECTION LINEAGES (WITH *2* replicates of ancestor.) == still >2G
```
cat novo_episodic_2R.sync | awk '{print $1,$2,$3,$6,$7,$10, $11, $14, $15, $16, $16}' > novo_episodic_2R_Sel.sync
```
```
cat novo_episodic_2R.sync | awk '{print $1,$2,$3,$4,$5,$8, $9, $12, $13, $16, $16}' > novo_episodic_2R_Con.sync
```

Not needed:
```
cat novo_episodic_2R.sync | awk '{print $1,$2,$3,$6,$10,$14,$16}' > novo_episodic_2R_SelR1.sync
```
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_2R_SelR1.sync /Users/paulknoops/Bioinformatics/episodic_practice/
```


```
#! /bin/bash

length=($(wc -l /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync))
echo ${length}

#Split length into 8 segements (8th == length) (can extend this if to large)
cut=$((${length}/8))
cut_2=$((${cut}*2))
cut_3=$((${cut}*3))
cut_4=$((${cut}*4))
cut_5=$((${cut}*5))
cut_6=$((${cut}*6))
cut_7=$((${cut}*7))

sed -n " 1, ${cut} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_1.sync

sed -n " $((${cut} + 1)), ${cut_2} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_2.sync

sed -n " $((${cut_2} + 1)), ${cut_3} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_3.sync

sed -n " $((${cut_3} + 1)), ${cut_4} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_4.sync

sed -n " $((${cut_4} + 1)), ${cut_5} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_5.sync

sed -n " $((${cut_5} + 1)), ${cut_6} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_6.sync

sed -n " $((${cut_6} + 1)), ${cut_7} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_7.sync

sed -n " $((${cut_7} + 1)), ${length} p" /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1.sync > /Users/paulknoops/Bioinformatics/episodic_practice/novo_episodic_2R_SelR1_8.sync

```


In case I need to regrep:
```
grep '3R' /home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_main.sync > /home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_3R.sync
```

## PoolSeq

need high version of R


### Loop for running Poolseq (Taus) -- breaking apart files: 

Wrong but keeping as reference: do in two steps:
```
#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup

# The seperated .sync files
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
	
	mkdir ${SyncFiles}/splitsync_dir
	splitSync=${SyncFiles}/splitsync_dir
	
	mkdir ${splitSync}/${base}_SelSync
	SelSync=${splitSync}/${base}_SelSync
	
	mkdir ${splitSync}/${base}_ConSync
	ConSync=${splitSync}/${base}_ConSync
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$6,$7,$10, $11, $14, $15, $16, $16}' > ${splitSync}/${base}_Sel.sync
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$4,$5,$8, $9, $12, $13, $16, $16}' > ${splitSync}/${base}_Con.sync
	
	SplitSync[1]=${splitSync}/${base}_Sel.sync
	SplitSync[2]=${splitSync}/${base}_Con.sync
	
	for file2 in ${SplitSync[@]}
		do
		name2=${file2}
		base2=`basename ${name2} .sync`
		
		length=($(wc -l ${base2}.sync))
		#echo ${length}

		#Split length into 8 segements (8th == length) (can extend this if to large)
		cut=$((${length}/8))
		cut_2=$((${cut}*2))
		cut_3=$((${cut}*3))
		cut_4=$((${cut}*4))
		cut_5=$((${cut}*5))
		cut_6=$((${cut}*6))
		cut_7=$((${cut}*7))

		sed -n " 1, ${cut} p" ${base2}.sync > ${base2}_1.sync

		sed -n " $((${cut} + 1)), ${cut_2} p" ${base2}.sync > ${base2}_2.sync

		sed -n " $((${cut_2} + 1)), ${cut_3} p" ${base2}.sync > ${base2}_3.sync

		sed -n " $((${cut_3} + 1)), ${cut_4} p" ${base2}.sync > ${base2}_4.sync

		sed -n " $((${cut_4} + 1)), ${cut_5} p" ${base2}.sync > ${base2}_5.sync

		sed -n " $((${cut_5} + 1)), ${cut_6} p" ${base2}.sync > ${base2}_6.sync

		sed -n " $((${cut_6} + 1)), ${cut_7} p" ${base2}.sync > ${base2}_7.sync

		sed -n " $((${cut_7} + 1)), ${length} p" ${base2}.sync > ${base2}_8.sync
		
		rm -f ${base2}.sync
	
	done

done
```

Split Treatments 1st:
```
#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup
	
mkdir ${SyncFiles}/splitsync_dir
splitSync=${SyncFiles}/splitsync_dir
	
# The seperated .sync files
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
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$6,$7,$10, $11, $14, $15, $16, $16}' > ${splitSync}/${base}_Sel.sync
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$4,$5,$8, $9, $12, $13, $16, $16}' > ${splitSync}/${base}_Con.sync

done
```


```
#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup

#Path to treatment split sync files
splitSync=${SyncFiles}/splitsync_dir

#All files in splitSync need to be split further:

files=(${splitSync}/*.sync)

for file in ${files[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	
	mkdir ${splitSync}/${base}_Split
	split_sync=${splitSync}/${base}_Split
	
	length=($(wc -l ${splitSync}${base}.sync))
	#echo ${length}
		
	#Split length into 8 segements (8th == length) (can extend this if to large)
	cut=$((${length}/8))
	cut_2=$((${cut}*2))
	cut_3=$((${cut}*3))
	cut_4=$((${cut}*4))
	cut_5=$((${cut}*5))
	cut_6=$((${cut}*6))
	cut_7=$((${cut}*7))
	
	sed -n " 1, ${cut} p"  ${splitSync}/${base}.sync > ${split_sync}${base}_1.sync

	sed -n " $((${cut} + 1)), ${cut_2} p"  ${splitSync}/${base}.sync >  ${split_sync}/${base}_2.sync

	sed -n " $((${cut_2} + 1)), ${cut_3} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_3.sync
	
	sed -n " $((${cut_3} + 1)), ${cut_4} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_4.sync

	sed -n " $((${cut_4} + 1)), ${cut_5} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_5.sync

	sed -n " $((${cut_5} + 1)), ${cut_6} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_6.sync

	sed -n " $((${cut_6} + 1)), ${cut_7} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_7.sync

	sed -n " $((${cut_7} + 1)), ${length} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_8.sync
	
done

```

```
#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup

#Path to treatment split sync files
splitSync=${SyncFiles}/splitsync_dir

#Output dir:
mkdir ${project_dir}/novo_PoolSeq
poolSeq=${project_dir}/novo_PoolSeq

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts


sync=(${splitSync}/*.sync)

for file in ${sync[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	basedir=${subsets}/${base}_Split
	Chromo=$(cat ${file} | awk '{print $1; exit}')
	/usr/bin/Rscript ${Rscripts}/PoolSeq_SelCoeff.R ${basedir} ${Chromo} ${poolSeq}	
done

```

 
for one chromosome:
novo_episodic_2L_Sel_Split
```
Rscript /home/paul/episodicData/novoalign/novo_Rscripts/PoolSeq_SelCoeff_one.R novo_episodic_2L_Sel_Split 2L /home/paul/episodicData/novoalign/novo_PoolSeq
```
/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus

R SCRIPT TO BE USED:

3 Arguments:
	1 == input (location of files)
	2 == Chromosome (identifies above)
	3 == output


### Script == PoolSeq_SelCoeff.R
```
### Running the PoolSeq Package to run on .sync files for estimates of selection coefficeints per position

### Requires: R (>= 3.3.1), data.table (>= 1.9.4), foreach (>= 1.4.2), stringi (>= 0.4-1), matrixStats (>= 0.14.2)

### 

### Required Packages:

  #install.packages("/home/paul/poolSeq_0.3.2.tar.gz", repos=NULL, type="source"

  require(poolSeq)

### Possible helpful help files:
  
  #?read.sync
  #?estimateSH
  #?af
  #?af.traj
  #??`poolSeq-package`
  
### These are part of the dependencies for poolSeq

  #require(data.table)
  #require(foreach)
  #require(stringi)
  #require(matrixStats)

### Not actually required: But used for plotting (ggplot) 

  #require(tidyverse)


### Possibly need custom function to read in manipulated .sync files:
  
  source("/home/paul/episodicData/novoalign/novo_Rscripts/Taus_ReadSync.R")
  
  args <- commandArgs(trailingOnly = TRUE)

### Reading in multiple sync file:
  
### Location of split .sync files:
  
  setwd(args[1])
  
  print(paste("Running directory:", dir))

### List the .sync files:
  
  SyncList <- list.files(pattern=".sync")
  
### Create empty data frame to hold all chromosomal information:
  
  DFULL <- data.frame(NULL)
  
### Loop through list of files:
  
  for (SyncFile in SyncList){

    mySync <- read.sync_Personal(file=SyncFile, gen=c(115, 115, 38, 38, 77, 77, 0, 0), repl=c(1,2,1,2,1,2,1,2), polarization = "rising")

    
### Make data.frame of just alleles information to sort out relevent positions:
    
  ff <- as.data.frame(mySync@alleles)
  pst <- as.numeric(ff$pos)
  pst2 <- sort(pst)
  rm(pst)
  rm(ff)

### Create empty data frame to read into for estiamting S:
  
  DF <- data.frame(NULL)
  ccc <- c(0,38,77,115)
  
### For Test:
  #len <- length(pst2)
  #pst2 <- sample(pst2[1]:pst2[len], 50)
  
for (i in pst2) {
  Traj115 <- af.traj(mySync, args[2], repl=c(1,2), pos=i)
  Bfsf <- estimateSH(Traj115, Ne=150, t=ccc, h=0.5, haploid = FALSE, simulate.p.value=TRUE)
  Fd <- data.frame(Bfsf$s, Bfsf$p0, Bfsf$p.value)
  Fd$pos <- i
  DF <- rbind(DF, Fd)
  DF <- na.omit(DF)
  #print(paste("Running entity:", i, "of", END))
  rm(i)
  
}
  
  DFULL <- rbind(DFULL, DF)
  
  rm(DF)
  rm(mySync)
  rm(ccc)
  rm(pst2)
  }
  
  
  X <- args[1]
  X3 <- gsub('.{7}$', '', X)
  
  write.csv(DFULL, file=paste(args[3], "/", X3, "SelCoeff.csv", sep=""), row.names=FALSE)

  
```

for one chromosome:
novo_episodic_2L_Sel_Split
```
Rscript /home/paul/episodicData/novoalign/novo_Rscripts/PoolSeq_SelCoeff_one.R /home/paul/episodicData/novoalign/novo_mpileup/splitsync_dir/novo_episodic_2L_Sel_Split 2L /home/paul/episodicData/novoalign/novo_PoolSeq
```
/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus


### Script == PoolSeq_SelCoeff_one.R
```
### Running the PoolSeq Package to run on .sync files for estimates of selection coefficeints per position

### Requires: R (>= 3.3.1), data.table (>= 1.9.4), foreach (>= 1.4.2), stringi (>= 0.4-1), matrixStats (>= 0.14.2)

  args <- commandArgs(trailingOnly = TRUE)

### Required Packages:

  #install.packages("/home/paul/poolSeq_0.3.2.tar.gz", repos=NULL, type="source")
  #install.packages("/home/paul/matrixStats_0.53.0.tar.gz", repos=NULL, type="source")
  
  #Not available: so source seperate:
  #require(poolSeq)
  
  ### These are part of the dependencies for poolSeq
  
  require(methods)
  require(data.table)
  require(foreach)
  require(stringi)
  require(matrixStats)
  
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/loadaf.R')  
  #estne.R Fails.
  #source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estne.R')
  
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estsh.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/idsel.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/simaf.R')


### Possible helpful help files:
  
  #?read.sync
  #?estimateSH
  #?af
  #?af.traj
  #??`poolSeq-package`
  

### Possibly need custom function to read in manipulated .sync files:
  
  source("/home/paul/episodicData/novoalign/novo_Rscripts/Taus_ReadSync.R")
  
### Reading in multiple sync file:
  
### Location of split .sync files:
  

### List the .sync files:
  
  SyncList <- list.files(path = args[1], pattern=".sync")
### Create empty data frame to hold all chromosomal information:
  
  DFULL <- data.frame(NULL)
  
### Loop through list of files:
  
  for (SyncFile in SyncList){

    mySync <- read.sync_Personal(file=Pasttg/SyncFile, gen=c(115, 115, 38, 38, 77, 77, 0, 0), repl=c(1,2,1,2,1,2,1,2), polarization = "rising")

    
### Make data.frame of just alleles information to sort out relevent positions:
    
  ff <- as.data.frame(mySync@alleles)
  pst <- as.numeric(ff$pos)
  pst2 <- sort(pst)
  rm(pst)
  rm(ff)

### Create empty data frame to read into for estiamting S:
  
  DF <- data.frame(NULL)
  ccc <- c(0,38,77,115)
  
### For Test:
  #len <- length(pst2)
  #pst2 <- sample(pst2[1]:pst2[len], 50)
  
for (i in pst2) {
  Traj115 <- af.traj(mySync, args[2], repl=c(1,2), pos=i)
  Bfsf <- estimateSH(Traj115, Ne=150, t=ccc, h=0.5, haploid = FALSE, simulate.p.value=TRUE)
  Fd <- data.frame(Bfsf$s, Bfsf$p0, Bfsf$p.value)
  Fd$pos <- i
  DF <- rbind(DF, Fd)
  DF <- na.omit(DF)
  #print(paste("Running entity:", i, "of", END))
  rm(i)
  
}
  
  DFULL <- rbind(DFULL, DF)
  
  rm(DF)
  rm(mySync)
  rm(ccc)
  rm(pst2)
  }
  
  
  X <- args[1]
  X3 <- gsub('.{7}$', '', X)
  
  write.csv(DFULL, file=paste(args[3], "/", X3, "SelCoeff.csv", sep=""), row.names=FALSE)


```


Run Fst and create plots

FST Output trial One: Min coverage == 50, why there is all na..
```
2L      3810250 3       0.798   58.2    1:2=na  1:3=na  1:4=na  1:5=na  1:6=na  1:7=na  1:8=na  1:9=na  1:10=na 1:11=na 1:12=na 1:13=na 2:3=na  2:4=na  2:5=na  2:6=na  2:7=na  2:8=na  2:9=na  2:10=na 2:11=na 2:12=na 2:13=na 3:4=na  3:5=na  3:6=na  3:7=na  3:8=na  3:9=na  3:10=na 3:11=na 3:12=na 3:13=na 4:5=na  4:6=na  4:7=na  4:8=na  4:9=na  4:10=na 4:11=na 4:12=na 4:13=na 5:6=na  5:7=na  5:8=na  5:9=na  5:10=na 5:11=na 5:12=na 5:13=na 6:7=na  6:8=na  6:9=na  6:10=na 6:11=na 6:12=na 6:13=na 7:8=na  7:9=na  7:10=na 7:11=na 7:12=na 7:13=na 8:9=na  8:10=na 8:11=na 8:12=na 8:13=na 9:10=na 9:11=na 9:12=na 9:13=na 10:11=na 10:12=na 10:13=na 11:12=na 11:13=na 12:13=na

2L      3810750 4       0.602   54.2    1:2=na  1:3=na  1:4=na  1:5=na  1:6=na  1:7=na  1:8=na  1:9=na  1:10=na 1:11=na 1:12=na 1:13=na 2:3=na  2:4=na  2:5=na  2:6=na  2:7=na  2:8=na  2:9=na  2:10=na 2:11=na 2:12=na 2:13=na 3:4=na  3:5=na  3:6=na  3:7=na  3:8=na  3:9=na  3:10=na 3:11=na 3:12=na 3:13=na 4:5=na  4:6=na  4:7=na  4:8=na  4:9=na  4:10=na 4:11=na 4:12=na 4:13=na 5:6=na  5:7=na  5:8=na  5:9=na  5:10=na 5:11=na 5:12=na 5:13=na 6:7=na  6:8=na  6:9=na  6:10=na 6:11=na 6:12=na 6:13=na 7:8=na  7:9=na  7:10=na 7:11=na 7:12=na 7:13=na 8:9=na  8:10=na 8:11=na 8:12=na 8:13=na 9:10=na 9:11=na 9:12=na 9:13=na 10:11=na        10:12=na        10:13=na        11:12=na        11:13=na        12:13=na
```
