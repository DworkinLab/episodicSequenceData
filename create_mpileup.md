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

```
#! /bin/bash

project_dir = /home/paul/episodicData

raw_dir = ${project_dir}/raw_dir

trimmomatic = /usr/local/trimmomatic
# location of trimmomatic on machine
adapt_path = /usr/local/trimmomatic/adapters
# path to adapter sequences
# ?? need to change the apater type!!!
bwa_path = /usr/local/bwa/0.7.8

# ?? ref_genome needs to be difined:

```



### Should not need to change anything below here

Create all working Directories
```
#! /bin/bash

cd {project_dir}

mkdir ${project_dir}/trim_dir

mkdir ${project_dir}/index_dir

mkdir ${project_dir}/bwa_dir


mkdir ${project_dir}/sam_dir


mkdir ${project_dir}/bam_dir

```

Defining all directories (cp to start of all scripts?)
```
project_dir = /home/paul/episodicData
raw_dir = ${project_dir}/raw_dir
trim_dir = ${project_dir}/trim_dir
index_dir = ${project_dir}/index_dir
bwa_dir = ${project_dir}/bwa_dir
sam_dir = ${project_dir}/sam_dir
bam_dir = ${project_dir}/bam_dir 

```




