#!/bin/bash
#SBATCH -J duplicator
#SBATCH -N 1
#SBATCH -t 148:00:00
#SBATCH -n 8
#SBATCH -A snic001-12-216

echo "env_duplicator_make.sbatch"

#
# 
#

/software/testing/R/2.15.1/bin/R --no-save < ~/fantom5/R/duplicator/env_duplicator_make.R

# End of script
