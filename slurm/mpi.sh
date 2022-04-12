#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=05:00
#SBATCH --partition=slow

srun /home/cwlui/cmpt431-project/linear_eqn_solver_mpi 