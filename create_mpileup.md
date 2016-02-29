# Full run through of all sequence data - Trimming to mpileup file format
## Need a thourough log
### Change parameters here at top, rest should fall in place
Should run md5sum and fastqc seperatly (before running quality control)
  - set it up so create project name and create a raw directory (raw_dir) and automatically create other directories
  - need known path to project name (i.e /home/paul/episodicData)
  - need to move paths to other directories on machine (i.e bwa or trim) at the top

project_dir = 
raw_dir = ${project_dir}/raw_dir
mkdir ${project_dir}/trim_dir
trim_dir = ${project_dir}/trim_dir
trimmomatic = 
  - location of trimmomatic on machine
adapt_path = 
  - path to adapter sequences
files =
  - files to input
index_dir = 
  - location of index directory (to put reference genome into)
bwa_dir = 
  - 
ref_genome = 
  - genome from index_dir for reference (will have the extension)
sam_dir = 
bam_dir
