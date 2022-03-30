#include <bits/stdc++.h>
#include <stdio.h>

#include <atomic>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>

#include "core/utils.h"

#define size 1000

// function to reduce matrix to r.e.f.  Returns a value to
// indicate whether matrix is singular or not
void forwardElim(double** mat, int tid, int start, int end, double* time_taken,
                 CustomBarrier* barrier, std::atomic<bool>& singular_flag);

// function to calculate the values of the unknowns
void backSub(double** mat);

// function to get matrix content
void gaussian_elimination_parallel(double** mat, uint n_threads) {
  std::vector<std::thread> threads(n_threads);

  // Dividing up rows for n threads
  std::vector<int> start_row(n_threads, 0);
  std::vector<int> end_row(n_threads, 0);
  int min_rows_per_process = size / n_threads;
  int excess_rows = size % n_threads;
  int cur_row = 0;
  for (int i = 0; i < n_threads; i++) {
    start_row[i] = cur_row;
    if (i < excess_rows) {
      end_row[i] = start_row[i] + min_rows_per_process;
    } else {
      end_row[i] = start_row[i] + min_rows_per_process - 1;
    }
    cur_row = end_row[i] + 1;
  }

  std::vector<double> local_time_taken(n_threads, 0.0);
  CustomBarrier barrier(n_threads);
  std::atomic<bool> singular_flag(false);
  timer t1;
  double time_taken = 0.0;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  for (int i = 0; i < n_threads; i++) {
    threads.push_back(std::thread(forwardElim, mat, i, start_row[i], end_row[i],
                                  &local_time_taken[i], &barrier,
                                  std::ref(singular_flag)));
  }

  for (std::thread& t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  if (singular_flag == true) {
    printf("Singular Matrix.\n");
    /* if the RHS of equation corresponding to
       zero row  is 0, * system has infinitely
       many solutions, else inconsistent*/
    if (mat[singular_flag][size])
      printf("Inconsistent System.");
    else
      printf(
          "May have infinitely many "
          "solutions.");

    return;
  }

  backSub(mat);

  time_taken = t1.stop();
  // -------------------------------------------------------------------

  // Print Statistics
  std::cout << "thread_id, starting row, ending row, time_taken" << std::endl;
  for (uint i = 0; i < n_threads; i++) {
    std::cout << i << ", " << start_row[i] << ", " << end_row[i] << ", "
              << local_time_taken[i] << std::endl;
  }
  std::cout << "Total Time taken (in seconds) : " << time_taken << "\n";
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

// function to reduce matrix to r.e.f.
void forwardElim(double** mat, int tid, int start, int end, double* time_taken,
                 CustomBarrier* barrier, std::atomic<bool>& singular_flag) {
  timer t;
  t.start();

  for (int pivot = 0; pivot < size; pivot++) {
    // T0 needs to make sure pivot is non zero
    if (mat[pivot][pivot] == 0 && tid == 0) {
      // swap with non zero row
      int swapIndex = -1;
      for (int i = pivot + 1; i < size; i++) {
        if (mat[i][pivot] != 0) {
          swapIndex = i;
        }
      }
      if (swapIndex == -1) {
        singular_flag = true;
        barrier->wait();
        break;  // Matrix is singular
      }

      swap_row(mat, pivot, swapIndex);
    }

    // If pivot is 0, all other threads need to wait until T0 swaps the rows
    barrier->wait();
    if (singular_flag == true) break;

    if (pivot + 1 <= end) {
      int currRow = (start > pivot) ? start : pivot + 1;
      for (currRow; currRow <= end; currRow++) {
        // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to
        // 0
        double f = mat[currRow][pivot] / mat[pivot][pivot];
        // TODO: Make this calculation local and then replace the entire row
        for (int j = pivot + 1; j <= size; j++) {
          mat[currRow][j] -= mat[pivot][j] * f;
        }
        mat[currRow][pivot] = 0;
      }
    }
    barrier->wait();
  }

  *time_taken = t.stop();
}

// function to calculate the values of the unknowns
void backSub(double** mat) {
  double x[size];  // An array to store solution

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

  printf("\nSolution for the system:\n");
  for (int i = 0; i < size; i++) printf("%lf\n", round(x[i]));
}

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "linear_eqn_solver_parallel",
      "Solver for a linear equation using parallel thread execution");
  options.add_options(
      "",
      {
          {"nThreads", "Number of Threads",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "inputs/generated.txt")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

  // Check edge cases on inputs
  if (n_threads <= 0) {
    throw std::invalid_argument(
        "The commandline arguments: --nThreads must be at least 1\n");
  }

  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout << "Creating Empty 1000 x 1000 Matrix\n";

  // allocate input matrix
  double** mat = new double*[size];
  mat[0] = new double[size * (size + 1)];
  for (int i = 1; i < size; ++i) {
    mat[i] = mat[0] + i * (size + 1);
  }

  std::cout << "Reading Input File : " << input_file_path << "\n";
  // read input matrix from file
  std::ifstream myfile(input_file_path, std::ios_base::in);
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size + 1; c++) {
      myfile >> mat[r][c];
    }
  }

  std::cout << "1000 x 1000 Matrix Filled In\n";

  gaussian_elimination_parallel(mat, n_threads);
  delete[] mat;
  return 0;
}