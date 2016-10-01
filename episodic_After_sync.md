#Scripts to use with sync files

Have two sync files for different mapping (BWA and Bowtie)
- both have had sections removed (all with Het, all with U and all with dmel_mitochondrion_genome

BWA:
```
episodic_data_main.sync
```

Bowtie:
```
episodic_data_bowtie_main.sync
```

### CMH Test
For loop script
- need to change the populations based on comparisons wanted
- should be able to repeat populations (pop 13 for base comparisons)
- change the 15's to 13's for sure
- 

BWA:
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

perl ${cmh_test} --min-count 3 --min-coverage 5 --max-coverage 100 --population ${population} --input ${mpileup_dir}/${project_name}_main.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-40

done
```

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
mkdir ${mpileup_dir}/${population}
pop_dir=${mpileup_dir}/${population}

perl ${cmh_test} --min-count 3 --min-coverage 5 --max-coverage 100 --population ${population} --input ${mpileup_dir}/${project_name}_main.sync --output ${pop_dir}/${project_name}_${population}.cmh.txt

perl ${cmh_gwas} --input ${pop_dir}/${project_name}_${population}.cmh.txt --output ${pop_dir}/${project_name}_${population}.cmh.gwas --min-pvalue 1.0e-40

done
```

