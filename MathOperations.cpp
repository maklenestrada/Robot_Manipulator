#include "MathOperations.h"
#include <vector>
#include <cmath>

//Matrix-Vector Multiplication
void MathOperations::VectorMult(double ans[4], double T[4][4], double P[4])
{
    for(int i = 0; i < 3; i++){ //For loop over the T matrix rows
        for(int j = 0; j < 3; j++){//For loop over the P elements
            ans[i] += T[i][j] * P[j];
        }
    }
}

//Matrix-Matrix Multiplication
void MathOperations::MatrixMult(double ans[4][4], double T_L[4][4], double T_R[4][4])
{
    for(int i = 0; i < 4; i++){ //Iterates through rows of T_L matrix
        for(int j = 0; j < 4; j++){ //Iterates through columns of T_R matrix
            for(int k = 0; k < 4; k++){ //Performs dot product of row i of T_L with column j of T_R
                ans[i][j] += T_L[i][k] * T_R[k][j];
            }
        }
    }
}

//Inverse of T Matrix
void MathOperations::InverseTMatrix(double ans[4][4],double T[4][4])
{
    // Upper left 3x3 of result is transpose of upper left 3x3 of tran
    for (int i=0 ; i<3 ; ++i)
        for (int j =0 ; j<3 ; ++j)
            ans[i][j] = T[j][i] ;
    // Set the values for the last column of the result
    ans[3][0] = ans[3][1] = ans[3][2] = 0.0 ;
    ans[3][3] = 1.0 ;
    // Initialize the values of the last column of the result.
    ans[0][3] = ans[1][3] = ans[2][3] = 0.0 ;
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            ans[i][3] -= ans[i][j] * T[j][3] ;
}

//Computes cross product
void MathOperations::CrossProduct(double ans[3],double X[3],double Y[3])
{
    ans[0] = X[1] * Y[2] - X[2] * Y[1];
    ans[1] = X[2] * Y[0] - X[0] * Y[2];
    ans[2] = X[0] * Y[1] - X[1] * Y[0];
}

//Computes dot product
void MathOperations::DotProduct(double* ans,double X[3],double Y[3])
{
    if (ans != nullptr){
        *ans = X[0] * Y[0] + X[1] * Y[1] + X[2] * Y[2];
    }

}

//Calculates norm
void MathOperations::Norm(double* ans,double vec[3])
{
    if (ans != nullptr){
        *ans = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    }
}

//Multiply vector by a scalor
void MathOperations::MultVecByScalar(double ans[3], double vec[3], double c)
{
    for(int i = 0; i <3; i++){
        ans[i] = vec[i] * c;
    }
}

//Add two vectors together
void MathOperations::AddVectors(double ans[3], double x[3],double y[3])
{
    for(int i = 0; i < 3; i++){
        ans[i] = x[i] + y[i];
    }
}

//Check if value is near desired tolerance
int MathOperations::ValueNear(double val, double goal, double tol)
{
    if ((val > goal - tol) && (val < goal + tol))
        return 1;
    else
        return 0;
}