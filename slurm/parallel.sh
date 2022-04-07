#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00
#SBATCH --mem=2G
#SBATCH --partition=slow

srun /home/cwlui/cmpt431-project/linear_eqn_solver_parallel --strategy 1 --nThreads 8 --printSolution 0 --granularity 20 --inputFile /home/cwlui/cmpt431-project/inputs/generated.txt