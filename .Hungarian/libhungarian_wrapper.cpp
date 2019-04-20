#include <iostream>
#include "Hungarian.h"

HungarianAlgorithm HA;
vector<vector<double>> Array2Vector(double* costArray, int rows, int cols); 
extern "C" void Linear_OptAssign(double* costArray, int rows, int cols, int* numRows, int* numCols); 

vector<vector<double>> Array2Vector(double* costArray, int rows, int cols) {
    // Transpose if rows > cols to meet the requirements of the algorithm
    bool isTransposed = false;
    if (rows > cols) { int tmp = cols; cols = rows; rows = tmp; isTransposed = true; }

    vector<vector<double>> costMatrix(rows);
    for (int i = 0; i < rows; i++) costMatrix[i].resize(cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) { 
            if (!isTransposed) costMatrix[i][j] = *(costArray + i * cols + j); 
            else costMatrix[i][j] = *(costArray + j * rows + i); 
        }
    }

    return costMatrix;
}

void Linear_OptAssign(double* costArray, int rows, int cols, int* numRows, int* numCols) {
    // Transpose if rows > cols to meet the requirements of the algorithm
    bool isTransposed = false;
    if (rows > cols) isTransposed = true;

    vector<vector<double>> costMatrix = Array2Vector(costArray, rows, cols);
    vector<int> assignment;
    HA.Solve(costMatrix, assignment); 

    if (!isTransposed) for (int i = 0; i < rows; i++) { numRows[i] = i; numCols[i] = assignment[i]; }
    else for (int i = 0; i < cols; i++) { numRows[i] = assignment[i]; numCols[i] = i; }

}


int main() {
    double testArray[] = { 10, 19, 8, 15, 0, 10, 18, 7, 17, 0, 13, 16, 9, 14, 0, 12, 19, 8, 18, 0 };
    int numRows[4] = {0};
    int numCols[4] = {0};

    Linear_OptAssign(testArray, 4, 5, numRows, numCols);
    return 0;
}