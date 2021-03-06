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

class DynamicMapping {
 public:
  uint k;
  uint n_rows;
  uint num_of_threads;
  uint pivot;
  std::atomic<uint> threads_done;
  std::atomic<uint> next_row;

  DynamicMapping()
      : k(1),
        n_rows(0),
        num_of_threads(1),
        pivot(1),
        threads_done(0),
        next_row(0) {}

  DynamicMapping(uint k, uint n, uint n_threads) {
    this->k = k;
    n_rows = n;
    num_of_threads = n_threads;
    threads_done = 0;
    next_row = 1;
    pivot = 1;
  }

  int getNextRowToBeProcessed() {
    uint cur_next = next_row.fetch_add(k);
    if (cur_next >= n_rows) {
      uint cur_threads = threads_done.fetch_add(1);
      if (cur_threads + 1 == num_of_threads) {
        // last thread should reset mapping state
        threads_done = 0;
        pivot++;
        next_row = pivot;
      }
      return -1;
    } else {
      return cur_next;
    }
  }
};

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

// function to reduce matrix to r.e.f.  Returns a value to
// indicate whether matrix is singular or not
// function to reduce matrix to r.e.f.
void forwardElimStatic(double** mat, int tid, int start, int end,
                       double* time_taken, CustomBarrier* barrier,
                       std::atomic<bool>& singular_flag, double* barrier1_time,
                       double* barrier2_time) {
  timer t;
  timer b1;
  timer b2;
  double b1_time = 0.0;
  double b2_time = 0.0;

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
    b1.start();
    barrier->wait();
    b1_time += b1.stop();

    if (singular_flag == true) break;

    if (pivot + 1 <= end) {
      int currRow = (start > pivot) ? start : pivot + 1;
      for (currRow; currRow <= end; currRow++) {
        // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to
        // 0
        double f = mat[currRow][pivot] / mat[pivot][pivot];
        for (int j = pivot + 1; j <= size; j++) {
          mat[currRow][j] -= mat[pivot][j] * f;
        }
        mat[currRow][pivot] = 0;
      }
    }
    b2.start();
    barrier->wait();
    b2_time += b2.stop();
  }

  *time_taken = t.stop();
  *barrier1_time = b1_time;
  *barrier2_time = b2_time;
}

void forwardElimDynamic(double** mat, int tid, double* time_taken,
                        CustomBarrier* barrier,
                        std::atomic<bool>& singular_flag, DynamicMapping* dm,
                        uint k, int* n_row_processed, double* barrier1_time,
                        double* barrier2_time) {
  timer t;
  timer b1;
  timer b2;
  double b1_time = 0.0;
  double b2_time = 0.0;
  int row_processed = 0;

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
    b1.start();
    barrier->wait();
    b1_time += b1.stop();
    if (singular_flag == true) break;

    while (true) {
      int currRow = dm->getNextRowToBeProcessed();
      // std::cout<<currRow<<std::endl;
      if (currRow == -1) break;
      for (uint j = 0; j < k; j++) {
        // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to
        // 0
        double f = mat[currRow][pivot] / mat[pivot][pivot];
        for (int j = pivot + 1; j <= size; j++) {
          mat[currRow][j] -= mat[pivot][j] * f;
        }
        mat[currRow][pivot] = 0;

        row_processed++;
        currRow++;
        if (currRow >= size) break;
      }
    }
    b2.start();
    barrier->wait();
    b2_time += b2.stop();
  }

  *time_taken = t.stop();
  *barrier1_time = b1_time;
  *barrier2_time = b2_time;
  *n_row_processed = row_processed;
}

void forwardElimEqual(double** mat, int tid, double* time_taken,
                      CustomBarrier* barrier, std::atomic<bool>& singular_flag,
                      int* n_row_processed, uint n_threads,
                      double* barrier1_time, double* barrier2_time) {
  timer t;
  timer b1;
  timer b2;
  double b1_time = 0.0;
  double b2_time = 0.0;
  int row_processed = 0;
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
    b1.start();
    barrier->wait();
    b1_time += b1.stop();
    if (singular_flag == true) break;

    // Determine how many rows to compute
    int min_rows_per_process = (size - pivot - 1) / n_threads;
    int start = (min_rows_per_process * tid) + pivot + 1;
    int end = (tid == n_threads - 1)
                  ? size-1
                  : (min_rows_per_process * (tid + 1)) + pivot;

    for (int currRow = start; currRow <= end; currRow++) {
      // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to
      // 0
      double f = mat[currRow][pivot] / mat[pivot][pivot];
      for (int j = pivot + 1; j <= size; j++) {
        mat[currRow][j] -= mat[pivot][j] * f;
      }
      mat[currRow][pivot] = 0;
      row_processed++;
    }
    b2.start();
    barrier->wait();
    b2_time += b2.stop();
  }

  *time_taken = t.stop();
  *barrier1_time = b1_time;
  *barrier2_time = b2_time;
  *n_row_processed = row_processed;
}

void gaussian_elimination_parallel_static(double** mat, uint n_threads,
                                          uint printSolution) {
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
  double x[size];  // An array to store solution
  double back_sub_time_taken = 0.0;
  std::vector<double> barrier1_time(n_threads, 0.0);
  std::vector<double> barrier2_time(n_threads, 0.0);

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  for (int i = 0; i < n_threads; i++) {
    threads.push_back(std::thread(forwardElimStatic, mat, i, start_row[i],
                                  end_row[i], &local_time_taken[i], &barrier,
                                  std::ref(singular_flag), &barrier1_time[i],
                                  &barrier2_time[i]));
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
      printf("May have infinitely many solutions.");
    return;
  }

  backSub(mat, x, &back_sub_time_taken);

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  // Print Statistics
  if (printSolution) {
    printf("\nSolution for the system:\n");
    for (int i = 0; i < size; i++) printf("x[%d]: %lf\n", i, round(x[i]));
  }

  // verify solution
  for (int i = 0; i < size; i++) {
    double sum = 0;
    for (int j = 0; j < size; j++) {
        sum += mat[i][j] * x[j];
    }
    assert(abs(mat[i][size] - sum) / mat[i][size] <= 0.005);
  }

  std::cout << "Solution Validated\n";
  std::cout << "Back Sub Time taken (in seconds) : " << back_sub_time_taken
            << "\n";
  std::cout << "thread_id, starting row, ending row, barrier1_time, "
               "barrier2_time, time_taken"
            << std::endl;
  for (uint i = 0; i < n_threads; i++) {
    std::cout << i << ", " << start_row[i] << ", " << end_row[i] << ", "
              << std::setprecision(TIME_PRECISION) << barrier1_time[i] << ", "
              << barrier2_time[i] << ", " << local_time_taken[i] << std::endl;
  }
  std::cout << "Total Time taken (in seconds) : " << time_taken << "\n";
}

void gaussian_elimination_parallel_dynamic(double** mat, uint n_threads, uint k,
                                           uint printSolution) {
  std::vector<std::thread> threads(n_threads);
  std::vector<double> local_time_taken(n_threads, 0.0);
  std::vector<int> rows_processed(n_threads, 0);
  CustomBarrier barrier(n_threads);
  DynamicMapping dm(k, size, n_threads);
  std::atomic<bool> singular_flag(false);
  timer t1;
  double time_taken = 0.0;
  double x[size];  // An array to store solution
  double back_sub_time_taken = 0.0;
  std::vector<double> barrier1_time(n_threads, 0.0);
  std::vector<double> barrier2_time(n_threads, 0.0);

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  for (int i = 0; i < n_threads; i++) {
    threads.push_back(
        std::thread(forwardElimDynamic, mat, i, &local_time_taken[i], &barrier,
                    std::ref(singular_flag), &dm, k, &rows_processed[i],
                    &barrier1_time[i], &barrier2_time[i]));
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
      printf("May have infinitely many solutions.");
    return;
  }

  backSub(mat, x, &back_sub_time_taken);

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  // Print Statistics
  if (printSolution) {
    printf("\nSolution for the system:\n");
    for (int i = 0; i < size; i++) printf("x[%d]: %lf\n", i, round(x[i]));
  }

  // verify solution
  for (int i = 0; i < size; i++) {
    double sum = 0;
    for (int j = 0; j < size; j++) {
        sum += mat[i][j] * x[j];
    }
    assert(abs(mat[i][size] - sum) / mat[i][size] <= 0.005);
  }

  std::cout << "Solution Validated\n";
  std::cout << "Back Sub Time taken (in seconds) : " << back_sub_time_taken
            << "\n";
  std::cout << "thread_id, num_rows, barrier1_time, barrier2_time, time_taken"
            << std::endl;
  for (uint i = 0; i < n_threads; i++) {
    std::cout << i << ", " << rows_processed[i] << ", "
              << std::setprecision(TIME_PRECISION) << barrier1_time[i] << ", "
              << barrier2_time[i] << ", " << local_time_taken[i] << std::endl;
  }
  std::cout << "Total Time taken (in seconds) : " << time_taken << "\n";
}

void gaussian_elimination_parallel_equal(double** mat, uint n_threads,
                                         uint printSolution) {
  std::vector<std::thread> threads(n_threads);
  std::vector<double> local_time_taken(n_threads, 0.0);
  std::vector<int> rows_processed(n_threads, 0);
  CustomBarrier barrier(n_threads);
  std::atomic<bool> singular_flag(false);
  timer t1;
  double time_taken = 0.0;
  double x[size];  // An array to store solution
  double back_sub_time_taken = 0.0;
  std::vector<double> barrier1_time(n_threads, 0.0);
  std::vector<double> barrier2_time(n_threads, 0.0);

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  for (int i = 0; i < n_threads; i++) {
    threads.push_back(
        std::thread(forwardElimEqual, mat, i, &local_time_taken[i], &barrier,
                    std::ref(singular_flag), &rows_processed[i], n_threads,
                    &barrier1_time[i], &barrier2_time[i]));
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
      printf("May have infinitely many solutions.");
    return;
  }

  backSub(mat, x, &back_sub_time_taken);

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  // Print Statistics
  if (printSolution) {
    printf("\nSolution for the system:\n");
    for (int i = 0; i < size; i++) printf("x[%d]: %lf\n", i, round(x[i]));
  }

  // verify solution
  for (int i = 0; i < size; i++) {
    double sum = 0;
    for (int j = 0; j < size; j++) {
        sum += mat[i][j] * x[j];
    }
    assert(abs(mat[i][size] - sum) / mat[i][size] <= 0.005);
  }

  std::cout << "Solution Validated\n";
  std::cout << "Back Sub Time taken (in seconds) : " << back_sub_time_taken
            << "\n";
  std::cout << "thread_id, num_rows, barrier1_time, barrier2_time, time_taken"
            << std::endl;
  for (uint i = 0; i < n_threads; i++) {
    std::cout << i << ", " << rows_processed[i] << ", "
              << std::setprecision(TIME_PRECISION) << barrier1_time[i] << ", "
              << barrier2_time[i] << ", " << local_time_taken[i] << std::endl;
  }
  std::cout << "Total Time taken (in seconds) : " << time_taken << "\n";
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
          {"strategy", "Task mapping strategy",
           cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
          {"granularity", "Row mapping Granularity",
           cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
          {"printSolution", "Toggle for solution printing",
           cxxopts::value<uint>()->default_value("1")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  uint strategy = cl_options["strategy"].as<uint>();
  uint k = cl_options["granularity"].as<uint>();
  uint printSolution = cl_options["printSolution"].as<uint>();

  // Check edge cases on inputs
  if (n_threads <= 0) {
    throw std::invalid_argument(
        "The commandline arguments: --nThreads must be at least 1\n");
  }

  if (strategy <= 0 || strategy > 3) {
    throw std::invalid_argument(
        "The commandline argument: --strategy only accepts values 1 and 2 \n"
        "    1 represents static mapping execution\n"
        "    2 represents dynamic mapping\n"
        "    3 represents equal mapping for every pivot");
  }

  if (k <= 0) {
    throw std::invalid_argument(
        "The commandline argument: --granularity must be a positive integer "
        "value\n");
  }

  if (printSolution > 1) {
    throw std::invalid_argument(
        "The commandline argument: --printSolution can only be 0 (printing "
        "off) or 1 (printing on)");
  }

  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout << "Strategy : " << strategy << std::endl;
  std::cout << "Granularity : " << k << std::endl;
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

  std::cout << "1000 x 1000 Matrix Read Into Memory\n";

  switch (strategy) {
    case 1:
      gaussian_elimination_parallel_static(mat, n_threads, printSolution);
      break;
    case 2:
      gaussian_elimination_parallel_dynamic(mat, n_threads, k, printSolution);
      break;
    case 3:
      gaussian_elimination_parallel_equal(mat, n_threads, printSolution);
      break;
    default:
      gaussian_elimination_parallel_static(mat, n_threads, printSolution);
  }
  delete[] mat;
  return 0;
}
