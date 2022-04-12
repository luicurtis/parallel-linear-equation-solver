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
#include <vector>

#define size 1000

using namespace std;

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
  double bTime = 0.0;
  double sTime = 0.0;
  double gTime = 0.0;
  timer commTimer;

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
    commTimer.start();
    MPI_Bcast(&(mat[pivot][pivot]), size+1-pivot, MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    bTime += commTimer.stop();
    int sendCount = min_rows_per_process*(size+1);
    commTimer.start();
    if (world_rank == world_size-1) {
      MPI_Scatter(&(mat[pivot+1][0]), sendCount, MPI_DOUBLE, MPI_IN_PLACE, sendCount, 
            MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    } else {
      MPI_Scatter(&(mat[pivot+1][0]), sendCount, MPI_DOUBLE, &(mat[start][0]), sendCount, 
            MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    }
    sTime += commTimer.stop();

    for (int currRow = start; currRow <= end; currRow++) {
      // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to 0
      double f = mat[currRow][pivot] / mat[pivot][pivot];
      for (int j = pivot + 1; j <= size; j++) {
        mat[currRow][j] -= mat[pivot][j] * f;
      }
      mat[currRow][pivot] = 0;
    }
    // gather finished sections of matrix at root
    commTimer.start();
    if (world_rank == world_size-1) {
      MPI_Gather(MPI_IN_PLACE, sendCount, MPI_DOUBLE, &(mat[pivot+1][0]), sendCount, 
          MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    } else {
      MPI_Gather(&(mat[start][0]), sendCount, MPI_DOUBLE, NULL, 0, 
          MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
    }
    gTime += commTimer.stop();
  }

  printf("Thread:%i, bTime:%f, sTime:%f, gTime:%f\n",world_rank,bTime,sTime,gTime);
}


void forwardEliminationStatic(double** mat, int world_rank, int world_size) {
  double bTime = 0.0;
  double cTime = 0.0;
  timer commTimer;

  // Determine which rows to work on
  int min_rows_per_process = size / world_size;
  int start = (min_rows_per_process * world_rank);
  int end = (world_rank == world_size - 1)
                ? size-1
                : (min_rows_per_process * (world_rank + 1))-1;

  // create communicators
  commTimer.start();
  int containsPivotRank = 0;
  vector<MPI_Comm> commVector;
  for (int i = 0; i < world_size-1; i++) {
    int color = (world_rank >= containsPivotRank) ? 1 : MPI_UNDEFINED;;
    MPI_Comm newComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &newComm);
    commVector.push_back(newComm);
  }
  cTime = commTimer.stop();

  for (int pivot = 0; pivot <= end; pivot++) {
    if (pivot >= start && pivot <= end) {
      // make sure pivot is non-zero 
      if (mat[pivot][pivot] == 0) {
        // swap with non-zero row
        int swapIndex = -1;
        for (int i = pivot + 1; i <= end; i++) {
          if (mat[i][pivot] != 0) {
            swapIndex = i;
          }
        }
        if (swapIndex == -1) {
          printf("Failed because matrix is singular!\n");
          MPI_Abort(MPI_COMM_WORLD, 1);  // Matrix is singular
        }
        swap_row(mat, pivot, swapIndex);
      }
    }

    // thread that contains the pivot broadcast the pivot 
    int currRank = min(world_size-1, pivot/min_rows_per_process);
    if (currRank != containsPivotRank) {
      containsPivotRank = currRank;
    }
    if (containsPivotRank != world_size-1) {
      MPI_Request req;
      commTimer.start();
      if (world_rank == containsPivotRank) {
        MPI_Ibcast(&(mat[pivot][pivot]), size+1-pivot, MPI_DOUBLE, containsPivotRank, commVector.at(containsPivotRank), &req);
      } else {
        MPI_Ibcast(&(mat[pivot][pivot]), size+1-pivot, MPI_DOUBLE, containsPivotRank, commVector.at(containsPivotRank), &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
      }
      bTime += commTimer.stop();
    }

    for (int currRow = max(start,pivot+1); currRow <= end; currRow++) {
      // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to 0
      double f = mat[currRow][pivot] / mat[pivot][pivot];
      for (int j = pivot + 1; j <= size; j++) {
        mat[currRow][j] -= mat[pivot][j] * f;
      }
      mat[currRow][pivot] = 0;
    }
  }

  printf("Thread:%i, bTime:%f, cTime:%f\n",world_rank,bTime,cTime);
}

void gaussian_elimination_mpi(double** mat, int world_rank, int world_size) {
  timer t1;
  double time_taken = 0.0;

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  forwardEliminationStatic(mat, world_rank, world_size);

  if (world_rank == world_size-1) {
    //print(mat);
    double x[size];  // An array to store solution
    double back_sub_time_taken = 0.0;
    backSub(mat, x, &back_sub_time_taken);
    // Print Statistics
    time_taken = t1.stop();
    cout << "Total Time taken (in seconds) : " << time_taken << "\n";

    printf("\nSolution for the system:\n");
    //for (int i = 0; i < size; i++) printf("x[%d]: %lf\n", i, round(x[i]));
    // verify solution
    for (int i = 0; i < 1000; i++) {
      assert(round(x[i]) == double(i * 2 - 100));
    }
    cout << "Solution Validated\n";
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
           cxxopts::value<string>()->default_value(
               "inputs/generated.txt")},
      });

  auto cl_options = options.parse(argc, argv);
  string input_file_path = cl_options["inputFile"].as<string>();

  MPI_Init(NULL, NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  //note that the root process is the last one (world_size-1)
  if (world_rank == world_size-1) {
    cout << "Creating Empty 1000 x 1000 Matrix\n";
  }

  // allocate input matrix
  double** mat = new double*[size];
  mat[0] = new double[size * (size + 1)];
  for (int i = 1; i < size; ++i) {
    mat[i] = mat[0] + i * (size + 1);
  }
  // read input matrix from file
  ifstream myfile(input_file_path, ios_base::in);
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size + 1; c++) {
      myfile >> mat[r][c];
    }
  }

  gaussian_elimination_mpi(mat, world_rank, world_size);
  delete[] mat;
  MPI_Finalize();
  return 0;
}