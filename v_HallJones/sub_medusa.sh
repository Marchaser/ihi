#!/bin/bash
#
# Execute the job from the current directory.
#$ -cwd
#
# Merge the standard error with the standard output.
#$ -j y
#
# Maximum runtime/wallclock. Please change it if your job requires more 
# than the system default, 24 hours.
#$ -l h_rt=24:00:00
#
# Maximum virtual memory per CPU core. Please change it if your job 
# requires more than the system default, 3G.
#$ -l h_vmem=3G
#
#$ -N base
#
 
# Set up proper modules.
. /etc/profile
module load intel/2013sp1
module load matlab/R2014a
module list
 
# Run your batch/serial job.
matlab -nojvm -nodisplay -r MAIN > base.txt 
 
exit 0

