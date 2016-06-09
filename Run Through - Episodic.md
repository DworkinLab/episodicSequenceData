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
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/1misc/Trial
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=/home/paul/episodicData/bowtie/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
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




### Runthrough
Run through rest of scripts the same
Change the project_dir and make sure names will be different in the end...










