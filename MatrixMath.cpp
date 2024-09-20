#include "MatrixMath.h"
#include <vector>

//Matrix-Vector Multiplication
void MatrixMath::VectorMult(double ans[4], double T[4][4], double P[4])
{
    for(int i = 0; i < 3; i++){ //For loop over the T matrix rows
        for(int j = 0; j < 3; j++){//For loop over the P elements
            ans[i] += T[i][j] * P[j];
        }
    }
}

//Matrix-Matrix Multiplication
void MatrixMath::MatrixMult(double ans[4][4],double T_L[4][4],double T_R[4][4])
{
    for(int i = 0; i < 4; i++){ //Iterates through rows of T_L matrix
        for(int j = 0; j < 4; j++){ //Iterates through columns of T_R matrix
            for(int k = 0; k < 4; k++){ //Performs dot product of row i of T_L with column j of T_R
                ans[i][j] += T_L[i][k] * T_R[k][j];
            }
        }
    }
}
