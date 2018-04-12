##########################
## basic submit command ##
## Audrey Hendricks     ##
## April 3 2018         ##
##########################

#!/bin/bash
#submit with  qsub -q all.q basic_submit.sh
#qstat -u yourusername

#$ -cwd
#$ -o logfilename.log
#$ -e errorfilename.err
#$ -S /bin/bash

Rscript Rcodetosubmitname.R

