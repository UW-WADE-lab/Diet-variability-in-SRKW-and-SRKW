#!/bin/bash
#SBATCH --partition=node # Queue selection
#SBATCH --job-name=Oorc_2016-2021_fecal_dada2_QAQC_tax95_halibut # Job name
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amy.vancise@noaa.gov # Where to send mail\
#SBATCH --ntasks=1 # Run a single task\
#SBATCH --cpus-per-task=16  # Number of CPU cores per task\
#SBATCH --mem=100000 # Job memory request\
#SBATCH --time=48:00:00 # Time limit hrs:min:sec\
#SBATCH --output=Oorc_2016-2021_fecal_dada2_QAQC_tax95_halibut.log # Standard output/error\
#export OMP_NUM_THREADS=8

module load R

#run dada2 script
Rscript dada2QAQC.R
