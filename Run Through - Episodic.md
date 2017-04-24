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



Test for running in paralle: have two data sets for F115ConR1 and F115ConR2 bam files to be merged
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


Getting CRISP onto Brians Machine:
1) download crisp onto local machine (.tar.gz): download from https://bansal-lab.github.io/software/crisp.html
2) SCP to remote location (``` scp CRISP-122713.tar.gz paul@info.mcmaster.ca:/home/paul ```)
3) Unpack file (```  tar xvzf CRISP-122713.tar.gz ```) 

Running CRISP: For Novoalign files test:
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

mkdir novo_crisp

FLAGS:
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
Need a list of all bam files with path names:
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
