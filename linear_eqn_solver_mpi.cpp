#include <bits/stdc++.h>
#include <stdio.h>
#include <mpi.h>

#include <atomic>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>

#include "core/utils.h"

#define size 1000

// function to calculate the values of the unknowns
void backSub(double** mat, double (&x)[size], double* back_sub_time_taken) {
  timer t;

  t.start();
  // calculate variables from bottom row to top row
  for (int i = size - 1; i >= 0; i--) {
    x[i] = mat[i][size];

    // Initialize j to i+1 since matrix is upper triangular
    for (int j = i + 1; j < size; j++) {
      // subtract all the lhs values except the coefficient of the variable
      // whose value is being calculated
      x[i] -= mat[i][j] * x[j];
    }

    // divide the RHS by the coefficient of the variable being calculated
    x[i] = x[i] / mat[i][i];
  }

  *back_sub_time_taken = t.stop();
}

// function for elementary operation of swapping two rows
void swap_row(double** mat, int i, int j) {
  for (int col = 0; col <= size; col++) {
    double temp = mat[i][col];
    mat[i][col] = mat[j][col];
    mat[j][col] = temp;
  }
}

void print(double** mat) {
  for (int i = 0; i < size; i++, printf("\n"))
    for (int j = 0; j <= size; j++) printf("%lf ", mat[i][j]);

  printf("\n");
}

void forwardElimination(double** mat, int world_rank, int world_size) {
  for (int pivot = 0; pivot < size; pivot++) {
    // root needs to make sure pivot is non zero
    if (mat[pivot][pivot] == 0 && world_rank == world_size-1) {
      // swap with non zero row
      int swapIndex = -1;
      for (int i = pivot + 1; i < size; i++) {
        if (mat[i][pivot] != 0) {
          swapIndex = i;
        }
      }
      if (swapIndex == -1) {
        MPI_Abort(MPI_COMM_WORLD, 1);  // Matrix is singular
      }
      swap_row(mat, pivot, swapIndex);
    }

    // Determine which rows to work on
    int min_rows_per_process = (size - pivot - 1) / world_size;
    int start = (min_rows_per_process * world_rank) + pivot + 1;
    int end = (world_rank == world_size - 1)
                  ? size-1
                  : (min_rows_per_process * (world_rank + 1)) + pivot;

    // root broadcasts the pivot row and the section of matrix to work on for each thread
    MPI_Bcast(&(mat[pivot][pivot]), size+1-pivot, MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    int sendCount = min_rows_per_process*(size+1);
    if (world_rank == world_size-1) {
      MPI_Scatter(&(mat[pivot+1][0]), sendCount, MPI_DOUBLE, MPI_IN_PLACE, sendCount, 
            MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    } else {
      MPI_Scatter(&(mat[pivot+1][0]), sendCount, MPI_DOUBLE, &(mat[start][0]), sendCount, 
            MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    }

    for (int currRow = start; currRow <= end; currRow++) {
      // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to 0
      double f = mat[currRow][pivot] / mat[pivot][pivot];
      for (int j = pivot + 1; j <= size; j++) {
        mat[currRow][j] -= mat[pivot][j] * f;
      }
      mat[currRow][pivot] = 0;
    }
    // gather finished sections of matrix at root
    if (world_rank == world_size-1) {
      MPI_Gather(MPI_IN_PLACE, sendCount, MPI_DOUBLE, &(mat[pivot+1][0]), sendCount, 
          MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    } else {
      MPI_Gather(&(mat[start][0]), sendCount, MPI_DOUBLE, NULL, 0, 
          MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    }
  }
}

void gaussian_elimination_mpi(double** mat, int world_rank, int world_size) {
  timer t1;
  double time_taken = 0.0;

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  forwardElimination(mat, world_rank, world_size);

  if (world_rank == world_size-1) {
    double x[size];  // An array to store solution
    double back_sub_time_taken = 0.0;
    backSub(mat, x, &back_sub_time_taken);
    // Print Statistics
    printf("\nSolution for the system:\n");
    for (int i = 0; i < size; i++) printf("x[%d]: %lf\n", i, round(x[i]));
    // verify solution
    for (int i = 0; i < 1000; i++) {
      assert(round(x[i]) == double(i * 2 - 100));
    }
    std::cout << "Solution Validated\n";
    time_taken = t1.stop();
    std::cout << "Total Time taken (in seconds) : " << time_taken << "\n";
  }
  // -------------------------------------------------------------------


}

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "linear_eqn_solver_parallel",
      "Solver for a linear equation using parallel thread execution");
  options.add_options(
      "",
      {
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "inputs/generated.txt")},
      });

  auto cl_options = options.parse(argc, argv);
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

  MPI_Init(NULL, NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  //note that the root process is the last one (world_size-1)
  if (world_rank == world_size-1) {
    std::cout << "Creating Empty 1000 x 1000 Matrix\n";
  }

  // allocate input matrix
  double** mat = new double*[size];
  mat[0] = new double[size * (size + 1)];
  for (int i = 1; i < size; ++i) {
    mat[i] = mat[0] + i * (size + 1);
  }

  if (world_rank == world_size-1) {
    std::cout << "Reading Input File : " << input_file_path << "\n";
    // read input matrix from file
    std::ifstream myfile(input_file_path, std::ios_base::in);
    for (int r = 0; r < size; r++) {
      for (int c = 0; c < size + 1; c++) {
        myfile >> mat[r][c];
      }
    }
    std::cout << "1000 x 1000 Matrix Read Into Memory\n";
  }

  gaussian_elimination_mpi(mat, world_rank, world_size);
  delete[] mat;
  MPI_Finalize();
  return 0;
}