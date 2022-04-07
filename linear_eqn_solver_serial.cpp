#include <bits/stdc++.h>
#include <stdio.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include "core/utils.h"

#define size 1000

// function to reduce matrix to r.e.f.  Returns a value to
// indicate whether matrix is singular or not
int forwardElim(double** mat);

// function to calculate the values of the unknowns
void backSub(double** mat, double (&x)[size], double *back_sub_time_taken);

// function to get matrix content
void gaussianElimination(double** mat, uint printSolution) {
  timer t1;
  double time_taken = 0.0;
  double x[size]; // An array to store solution
  double back_sub_time_taken = 0.0;

  // -------------------------------------------------------------------
  t1.start();
  int singular_flag = forwardElim(mat);
  if (singular_flag != -1) {
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
  backSub(mat, x, &back_sub_time_taken);

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  // verify solution
  for (int i = 0; i < 1000; i++) {
    assert(round(x[i]) == double(i*2 - 100));
  }

  // Print Statistics
  if (printSolution) {
    printf("\nSolution for the system:\n");
    for (int i = 0; i < size; i++) printf("x[%d]: %lf\n", i, round(x[i]));
  }

  std::cout << "Solution Validated\n";
  std::cout << "Back Sub Time taken (in seconds) : " << back_sub_time_taken << "\n";
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
int forwardElim(double** mat) {
  for (int pivot = 0; pivot < size; pivot++) {
    // make sure pivot is non zero
    if (mat[pivot][pivot] == 0) {
      // swap with non zero row
      int swapIndex = -1;
      for (int i = pivot + 1; i < size; i++) {
        if (mat[i][pivot] != 0) {
          swapIndex = i;
        }
      }
      if (swapIndex == -1) return pivot;  // Matrix is singular

      swap_row(mat, pivot, swapIndex);
    }

    for (int currRow = pivot + 1; currRow < size; currRow++) {
      // subtract f * pivot row from currRow to reduce mat[currRow][pivot] to
      // 0
      double f = mat[currRow][pivot] / mat[pivot][pivot];
      for (int j = pivot + 1; j <= size; j++) {
        mat[currRow][j] -= mat[pivot][j] * f;
      }
      mat[currRow][pivot] = 0;
    }
  }
  return -1;
}

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

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "linear_eqn_solver_serial",
      "Solver for a linear equation using serial execution");
  options.add_options("", {
                              {"inputFile", "Input graph file path",
                               cxxopts::value<std::string>()->default_value(
                                   "inputs/generated.txt")},
                              {"printSolution", "Toggle for solution printing",
                               cxxopts::value<uint>()->default_value("1")},
                          });

  auto cl_options = options.parse(argc, argv);
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  uint printSolution = cl_options["printSolution"].as<uint>();

  if (printSolution > 1) {
    throw std::invalid_argument(
      "The commandline argument: --printSolution can only be 0 (printing off) or 1 (printing on)");
  }

  std::cout << "Number of Threads : 1" << std::endl;
  std::cout << "Creating Empty 1000 x 1000 Matrix\n";
  // allocate input matrix
  double** mat = new double*[size];
  mat[0] = new double[size * (size + 1)];
  for (int i = 1; i < size; ++i) {
    mat[i] = mat[0] + i * (size + 1);
  }

  std::cout << "Reading Input File : " << input_file_path<< "\n";
  // read input matrix from file
  // std::ifstream myfile("data.txt", std::ios_base::in);
  std::ifstream myfile(input_file_path, std::ios_base::in);
  for (int r = 0; r < size; r++) {
    for (int c = 0; c < size + 1; c++) {
      myfile >> mat[r][c];
    }
  }

  std::cout << "1000 x 1000 Matrix Filled In\n";

  gaussianElimination(mat, printSolution);
  delete[] mat;
  return 0;
}