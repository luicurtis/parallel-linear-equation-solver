#include<bits/stdc++.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
  
#define size 3
  

int main() {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-100000.0, 110000.0);
    int* variables = new int[1000];
    for (int i = 0; i < 1000; i++) {
        variables[i] = i*2 - 100;
    }
    ofstream myfile;
    myfile.open ("generated.txt");
    for (int row = 0; row < 1000; row++) {
        double sum = 0;
        for (int col = 0; col < 1000; col++) {
            double coefficient = dist(mt);
            myfile << coefficient << " ";
            sum += coefficient * variables[col];
        }
        myfile << sum << " ";
    }

    return 0;
}