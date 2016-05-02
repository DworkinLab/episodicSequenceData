# Full run through of all sequence data - Trimming to mpileup file format
## 
## Change parameters here at top (before scripts), rest should fall in place, unless changes need to be made for quality purposes
____________________________________________________________________________________________________
  
###Step 1: make project_dir and put all files are in raw_dir=${project_dir}/raw_dir

###Step 2: md5sum all raw files: changes depending on the file name (example below)
```
md5sum - c md5.txt
```

###Step 3: Fastqc; run as a quality control and view

###Step 4: Create all working Directories and bringing in reference sequence and indexing


```
#! /bin/bash

#Change name of project directory (made above with raw_dir)

project_dir=/home/paul/episodicData

# Can change the index sequence here

mkdir ${project_dir}/index_dir
index_dir=${project_dir}/index_dir
cd ${index_dir}
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

bwa index dmel-all-chromosome-r5.57.fasta.gz

ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz



cd ${project_dir}

raw_dir=${project_dir}/raw_dir
mkdir ${project_dir}/trim_dir
mkdir ${project_dir}/index_dir
mkdir ${project_dir}/sam_dir
mkdir ${project_dir}/bam_dir
mkdir ${project_dir}/merged
mkdir ${project_dir}/sort_dir
mkdir ${project_dir}/tmp
mkdir ${project_dir}/rmd_dir
mkdir ${project_dir}/final_bam
mkdir ${project_dir}/mpileup_dir

```


###Def_Dir
Defining all directories (copy to start of all scripts)

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
sam_dir=${project_dir}/sam_dir
bam_dir=${project_dir}/bam_dir 
merged=${project_dir}/merged
sort_dir=${project_dir}/sort_dir
tmp=${project_dir}/tmp
rmd_dir=${project_dir}/rmd_dir
final_bam=${project_dir}/final_bam
mpileup_dir=${project_dir}/mpileup_dir

```

In order for easier run through, change names to common names (ending in _L001_RX_001.fastq.gz or _L002_RX_001.fastq.gz)
  - Change using mv function: easiest method would be to copy and paste with changes already made for ones you want
  - Change to match F38SelR2_GATCAG_L001_R2_001.fastq.gz style
  - should be same start
  - all L003 = L001
  - all L004 = L002 (etc.)

##Scripts: copy the file/directory paths above into each script
###Nano in scripts directory for each step below
__________________________________________________________________________________

### Trimmomatic
Flags:
  -phred33 = may not need to be specified
  -trimlog = log of trim outputs
  -IlluminaClip = adapter removal (${adapter})
  -LEADING & TRAILING = 3; removal at start end end if below quality
  -MINLEN = minimum length of 36
  -MAXINFO = adaptive quality (balance b/w length and quality) = 0.5


```
#! /bin/bash

files=(${raw_dir}/*_R1_001.fastq.gz)
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`
java -jar ${trim} PE -phred33 -trimlog ${trim_dir}/trimlog.txt ${raw_dir}/${base}_R1_001.fastq.gz ${raw_dir}/${base}_R2_001.fastq.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R1_SE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz ${trim_dir}/${base}_R2_SE.fastq.gz ILLUMINACLIP:${adapter} LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done
```

### BWA mapping
Flags:
  -t 8 = number of processors
  -M = Mark shorter split hits as secondary (for Picard compatibility)

```
#!/bin/bash

cd ${bwa_path}
files=(${trim_dir}/*_R1_PE.fastq.gz)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
bwa mem -t 8 -M ${ref_genome} ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${sam_dir}/${base}_aligned_pe.SAM
done
```

### convert SAM to BAM
Flags:
  -b = Bam files
  -S = Sam files
  -q 20 = minimum quality score 

```
#! /bin/bash


files=(${sam_dir}/*.SAM)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .SAM`
samtools view -b -S -q 20 ${sam_dir}/${base}.SAM | samtools sort - ${bam_dir}/${base}
done
```
### Merge files
Only works if all lanes are L001/L002

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
Flags: 
  -Xmx2g = allocated Java 2 Gb of memory
  -SO = sort order
  -VALIDATION_STRINGENCY = surpress validation messages completely
  -I = input file
  -O = output file
  -Djava.io.tmpdir=${tmp} = to be sure for java programs in general (may not be needed and TMP_DIR could be sufficient)
  -TMP_DIR = temporary scratch directory for large files
  
```
#! /bin/bash

files=(${merged}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
java -Xmx2g -Djava.io.tmpdir=${tmp} -jar ${pic} SortSam I= ${merged}/${base}.bam O= ${sort_dir}/${base}.sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${tmp}
done
```

### Remove duplicates
Flags: same as above

```
#! /bin/bash

files=(${sort_dir}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${sort_dir}/${base}.sort.bam O= ${rmd_dir}/${base}.rmd.sort.bam M= ${rmd_dir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
```

### Remove low quality reads
Flags:
  -q 20 = min quality score
  -F 0x0004 = remove unmapped reads
  - 
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

samtools mpileup -B -Q 0 -f ${ref_genome} ${final_bam}/*.bam > ${mpileup_dir}/${project_name}.mpileup
```

### Sync Files
```
#! /bin/bash

java -ea -Xmx7g -jar ${sync} --input ${mpileup_dir}/${project_name}.mpileup --output ${mpileup_dir}/${project_name}.sync --fastq-type sanger --min-qual 20 --threads 2
```

