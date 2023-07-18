#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:45:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=18   # 36 processor core(s) per node 
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --job-name="lokib"
#SBATCH --output="lokib.out"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load matlab/R2022a
export omp_num_threads=8
matlab -nodisplay -nosplash -nodesktop -r "clear all; run_pulse(2000); exit;"
