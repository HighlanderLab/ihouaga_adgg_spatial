#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N job2.sh              
#$ -cwd                  
#$ -l h_rt=00:05:00 
#$ -l h_vmem=1G
# Initialise the environment modules
# Load R
module load R/4.0
# Run the program
R --no-save --no-restore -f adgg_descriptive_version1.0.0.R