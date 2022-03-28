#include<bits/stdc++.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
  
#define size 1000
  
// function to reduce matrix to r.e.f.  Returns a value to 
// indicate whether matrix is singular or not
int forwardElim(double** mat);
  
// function to calculate the values of the unknowns
void backSub(double** mat);
  
// function to get matrix content
void gaussianElimination(double** mat) {
    int singular_flag = forwardElim(mat);
    if (singular_flag != -1) {
        printf("Singular Matrix.\n");
        /* if the RHS of equation corresponding to
           zero row  is 0, * system has infinitely
           many solutions, else inconsistent*/
        if (mat[singular_flag][size])
            printf("Inconsistent System.");
        else
            printf("May have infinitely many "
                   "solutions.");
  
        return;
    }
    backSub(mat);
}
  
// function for elementary operation of swapping two rows
void swap_row(double** mat, int i, int j) {
    for (int col=0; col<=size; col++) {
        double temp = mat[i][col];
        mat[i][col] = mat[j][col];
        mat[j][col] = temp;
    }
}
  
void print(double** mat) {
    for (int i=0; i<size; i++, printf("\n"))
        for (int j=0; j<=size; j++)
            printf("%lf ", mat[i][j]);
  
    printf("\n");
}
  
// function to reduce matrix to r.e.f.
int forwardElim(double** mat) {
    for (int piviot=0; piviot<size; piviot++) {
        // make sure pivot is non zero
        if (mat[piviot][piviot] == 0) {
            //swap with non zero row
            int swapIndex = -1;
            for (int i = piviot+1; i < size; i++) {
                if (mat[i][piviot] != 0) {
                    swapIndex = i;
                }
            }
            if (swapIndex == -1)
                return piviot; // Matrix is singular
                
            swap_row(mat, piviot, swapIndex);
        }
  
        for (int currRow=piviot+1; currRow<size; currRow++) {
            //subtract f * piviot row from currRow to reduce mat[currRow][piviot] to 0
            double f = mat[currRow][piviot]/mat[piviot][piviot];
            for (int j=piviot+1; j<=size; j++) {
                mat[currRow][j] -= mat[piviot][j]*f;
            }
            mat[currRow][piviot] = 0;
        }
    }
    return -1;
}
  
// function to calculate the values of the unknowns
void backSub(double** mat) {
    double x[size];  // An array to store solution
  
    //calculate variables from bottom row to top row
    for (int i = size-1; i >= 0; i--) {
        x[i] = mat[i][size];
  
        // Initialize j to i+1 since matrix is upper triangular
        for (int j=i+1; j<size; j++) {
            // subtract all the lhs values except the coefficient of the variable whose value is being calculated
            x[i] -= mat[i][j]*x[j];
        }
  
        // divide the RHS by the coefficient of the variable being calculated
        x[i] = x[i]/mat[i][i];
    }
  
    printf("\nSolution for the system:\n");
    for (int i=0; i<size; i++)
        printf("%lf\n", round(x[i]));
}
  

int main() {
    // allocate input matrix
    double** mat = new double*[size];
    mat[0] = new double[size * (size+1)];
    for (int i = 1; i < size; ++i) {
        mat[i] = mat[0] + i * (size+1);
    }
    // read input matrix from file
    //std::ifstream myfile("data.txt", std::ios_base::in);
    std::ifstream myfile("generated.txt", std::ios_base::in);
    for (int r = 0; r < size; r++) {
        for(int c = 0; c < size+1; c++) {
            myfile >> mat[r][c];
        }
    }
    gaussianElimination(mat);
    delete [] mat;
    return 0;
}