# Analysis with custom R scripts
________________________________________________________________________________________________________________

Layout of work done with R on generated .sync files

## Before R: Setting Up .sync Files

The .sync files are generated using Popoolation2 depending on the data cleaning done to generate .bam files 

The .sync files used can have the unused regions removed for R analysis: use grep
```
#Remove Heterochromatin regions ('Het'), the Unmapped Regions ('U'), and the Mitochondria genome ('dmel_mitochondrion_genome'):
# -v == does not contain (keeps those without the 'Het' etc.)
# ${sync} == path to the location of sync files
# ${project_name} == name of the .sync file that uniquely identifies it 

grep -v 'Het' ${sync}/${project_name}.sync > ${sync}/${project_name}_less_het.sync
grep -v 'U' ${sync}/${project_name}_less_het.sync > ${sync}/${project_name}_removed_U_Het.sync
grep -v 'dmel_mitochondrion_genome' ${sync}/${project_name}_removed_U_Het.sync > ${sync}/${project_name}_main.sync

#Remove intermediates (-f == force remove) to have main.sync files
rm -f ${sync}/${project_name}_less_het.sync
rm -f ${sync}/${project_name}_removed_U_Het.sync
```

Left with a .sync file with only the 3R, 3L, 2R, 2L, X and 4 chromosomal arms

Split these more: split into individual files based on chromosome using grep
```
# ${sync} == path to the location of sync files
# ${project_name} == name of the .sync file that uniquely identifies it 
# '&' runs each in parallel: all will run together not one at a time
# For 4th chromo: using ^4 which only keeps rows that begin with 4 (layout of .sync begins with chromosome)

grep '3R' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_3R.sync &
grep '2R' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_2R.sync &
grep '3L' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_3L.sync &
grep '2L' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_2L.sync &
grep '^4' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_4.sync &
grep 'X' ${sync}/${project_name}_main.sync > ${sync}/${project_name}_X.sync
```

.sync files for each different chromsome

Depending on the size, these may be able to run individually (likely possible with small 4th) but R only works if using smaller data sets: may need to break up further



## Running In R

### Packages:

Packages Needed through analysis: be sure to have installed on either local or remote machine based on place of analysis

```
#install.packages('dplyr')
library(dplyr)

#install.packages('tidyr')
library(tidyr)

#install.packages('ggplot2')
library(ggplot2)
```

### Converting Sync file to counts

Similar to some 



