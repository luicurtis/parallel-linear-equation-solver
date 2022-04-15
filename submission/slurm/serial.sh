#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00
#SBATCH --mem=2G
#SBATCH --partition=slow

srun ./linear_eqn_solver_serial --printSolution 0 --inputFile ./inputs/generated.txt
