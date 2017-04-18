# Running Novoalign mapper

Starting with sequence files that have already been inspected with md5sum and fastqc and have been trimmed using trimmomatic (See other files for running trimmomatic)

Novoalign link tutorial: http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/basic-short-read-mapping/

### Need to create directory for project (project_dir) and to house mapping outputs

ex. mkdir novoalign (project), and mkdir novo_dir (mapping outputs)

### Find path to call Novoalign

ex. /usr/local/novoalign


### Make a scripts directory to house scripts to run 

ex. mkdir novo_scripts

### Novoindex reference

Need index dir: mkdir novo_index

Example Script:
```
#! /bin/bash

#Create variable for location of reference genome
ref_genome=/home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta

#Variables for project, novoalign, and output directory
project_dir=/home/paul/episodicData/novoalign
novoalign=/usr/local/novoalign
novo_index=${project_dir}/novo_index

#Index the reference
${novoalign}/novoindex ${novo_index}/dmel-all-chromosome-r5.57_2.nix  ${ref_genome}
```

### Run Novoalign
Note: Compressed read files are not supported in unlicensed versions.

Need to unzip (in trimmomatic output dir)
```
gunzip *.gz
```

    > Novoalign Flags
        - d -- Full pathname of indexed reference sequence from novoindex
        - f -- Files containing the read sequences to be aligned  
        - o -- Specifies output report format and options (SAM)  
        - i ###,## -- Sets fragment orientation and approximate fragment length for proper pairs.
             ex. -i 250 50  Defaults to paired end Illumina or Mate Pair ABI with 250bp insert and 50bp standard deviation
                - Using 400,100 based on initial mapping with novoalign first run through


The script
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
${novoalign}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 400,100 -o SAM > ${novo_dir}/${base}_novo.sam

done
```

Takes a long time: Only uses 1 thread (100% computer)

From novoalign reference manual: -c 99 Sets the number of threads to be used. On licensed versions it defaults 
to the number of CPUs as reported by sysinfo(). On free version the option is disabled

To avoid this problem, run scripts for each in parallel

1) make the script to make multiple scripts

Make dir for all output scripts: mkdir split_mappingScripts

Run in Scripts dir ( not split_mappingScripts)

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
echo "${novoalign}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 400,100 -o SAM > ${novo_dir}/${base}_novo.sam" > ./split_mappingScripts/${base}.sh

done
```

2) Create script to call all and run in parallel (us & which puts job in background)
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
echo "${map_scripts}/${base}.sh &" > ${scripts}/novo_parallel_map.sh

done
```
Change permissions and run novo_parallel_map.sh to run all files at once (can change input parametes to run subsets on different days)

Rezip files in trim_dir (saves space)
```
gzip *.fastq
```

### sam-bam files: 
-- need bam directory (mkdir novo_bam)

    >Flags:
        - b -- output is .bam
        - S -- input is .sam
        - q 20 -- quality mapping score of 20 (standard throughout all experiments)
        

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


## Next Steps:

### Merge 
--  mkdir novo_merge
> no flags
> be sure to merge the base generation seperatly (two sequence runs)

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
-- mkdir novo_pic and mkdir novo_tmp (was helpful, explained in flags)

    > Flags; 
        - Xmx2g -- 2 Gb of memory allocated
        - Djava.io.tmpdir=${tmp} -- using my temporary direcory due to errors in space allocation to avoid errors while running (not necessary)
        - I -- input
        - O -- output
        - VALIDATION_STRINGENCY=SILENT -- stops Picard from reporting every issue that would ultimately be displayed
        - SO=coordinate -- sort order based on coordinate
    
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


files=(${novo_merge}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
java -Xmx2g -Djava.io.tmpdir=${novo_tmp} -jar ${pic} SortSam I= ${novo_merge}/${base}.bam O= ${novo_pic}/${base}_novo_sort.bam VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${novo_tmp}
done
```

### Remove Duplicates 
-- mkdir novo_rmd

    > Flags;
        - Xmx2g -- ""
        - MarkDuplicates -- ""
        - I -- ""
        - O -- ""
        - M -- creates an output file of statistics of duplicates found
        - VALIDATION_STRINGENCY=SILENT -- ""
        - REMOVE_DUPLICATES= true -- get rid of any found duplicated regions

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
-- mkdir novo_final

    > Flags;
        -q 20 -- ""
        -F 0x0004 -- remove any unmapped reads (hexidecimal value for unmapped = 0x0004)
        -b -- ""

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

### Create mpileup
-- mkdir novo_mpileup

    > Flags;
        -B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment
        -Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)
        -f -- path to reference sequence
        

- Needs the reference genome: indexed version or not?

```
#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_final=${project_dir}/novo_final

#Path to output directory
novo_mpileup=${project_dir}/novo_mpileup

#Variable for reference genome (non-indexed for novoAlign)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz

samtools mpileup -B -Q 0 -f ${ref_genome} ${novo_final}/*.bam > ${novo_mpileup}/${project_name}.mpileup
```

### Create .sync file 
--use mpileup dir

    > Flags;
        -Xmx7g -- ""
        --input -- ""
        --output -- ""
        --fastq-type -- needed for base encoding
        --min-qual 20 -- already set to 20 before, but since custome script, use 20 as safer assumption
        --threads 2 -- ""
        
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

### Testing out GATK
-- mkdir novo_GATK
-- unzip reference?
-- need to do as a loop with base again?

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


files=(${novo_final}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -Xmx8g -jar ${gatk} -I ${novo_final}/${base}.bam -R ${ref_genome} -T RealignerTargetCreator -o ${novo_GATK}/${base}.intervals

java -Xmx8g -jar ${gatk} -I ${novo_final}/${base}.bam -R ${ref_genome} -T IndelRealigner -targetIntervals ${novo_GATK}/${base}.intervals -o ${novo_GATK}/${base}_realigned.bam

done
```
Error with reference
need .dict (java -jar CreateSequenceDictionary.jar R= Homo_sapiens_assembly18.fasta O= Homo_sapiens_assembly18.dict)
also need .fai (have but rerun? -- samtools faidx Homo_sapiens_assembly18.fasta )


Needs a readgroup (dummy can work)
java -jar picard.jar AddOrReplaceReadGroups I=input.bam O=output.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

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

java -jar /usr/local/picard-tools-1.131/picard.jar AddOrReplaceReadGroups I=${novo_final}/${base}.bam O=${novo_final}/${base}_RG.bam RGID=L001_L002 RGLB=library1 RGPL=illumina RGPU=None RGSM=${base}



- samtools workflow with GATK:               http://www.htslib.org/workflow/
- CRISP:                     https://github.com/vibansal/crisp
