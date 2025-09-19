#!/bin/bash

#SBATCH --time=06:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out 
module load r 
R CMD BATCH NMA_simulation_both_with_prev_row_2.R
exit

