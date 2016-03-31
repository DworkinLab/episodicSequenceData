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

Define all directories 
```
#! bin/bash

trim_dir = ${project_dir}/trim_dir
index_dir = ${project_dir}/index_dir
bwa_dir = ${project_dir}/bwa_dir
sam_dir = ${project_dir}/sam_dir
bam_dir = ${project_dir}/bam_dir 

```




