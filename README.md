# How to Use and Generate Input Files
To genereate input files:
1. run `make`
2. execute `./generate`
3. A text file, `generated.txt` will be created and can be used as an input file to all programs

Note: If the option `--inputFile` is not specified, the default file that is used is `/inputs/generated.txt`. There are already 3 input files generated and located in `/inputs`

# How to Run Implementations
Run `make` to create binaries of all implementations.

## Running on Local Machine
### Serial Version
Input Arguments:

- `--inputFile`: specify the path to the desired input file
- `--printSolution`: Value of 0 or 1. 0 to have no solution printing and 1 to print the resulting solution vector values 

Example:

`./linear_eqn_solver_serial --inputFile inputs/generated.txt --printSolution 0`

### Parallel Version
Input Arguments:
- `--inputFile`: specify the path to the desired input file
- `--printSolution`: Value of 0 or 1. 0 to have no solution printing and 1 to print the resulting solution vector values 
- `--nThreads`: Number of threads. Default value of 1.
- `--strategy`: Value of 1, 2, or 3. Default value of 1.
    - 1 - Static Mapping
    - 2 - Dynamic Mapping
    - 3 - Equal Mapping
- `--granularity`: Granularity for dynamic mapping. Only affects `--strategy 2`. Default value of 1.

Example:

`./linear_eqn_solver_parallel --inputFile inputs/generated.txt --printSolution 0 --nThreads 4 --strategy 1 --granularity 1`

### Distributed Version:
Input Arguments:
- `--inputFile`: specify the path to the desired input file

Example:

`mpirun -n 1 ./linear_eqn_solver_mpi --inputFile inputs/generated.txt`

## Running on Cluster
### Serial Version:
Input Arguments: Same as above. Change values within the shell script

Within the project directory run: `sbatch slurm/serial.sh`

### Parallel Version:
Input Arguments: Same as above. Change values within the shell script

Within the project directory run: `sbatch slurm/parallel.sh`

### Distributed Version:
Input Arguments: Same as above. Change values within the shell script

Within the project directory run: `sbatch slurm/mpi.sh`




