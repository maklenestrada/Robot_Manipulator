#ifndef ROBOT_MANIPULATOR__MATRIXMATH_H
#define ROBOT_MANIPULATOR__MATRIXMATH_H
class MatrixMath {
public:
    //Constructor
    MatrixMath() {}

    //Function for matrix-vector multiplication (T (4x4) times P (4x1))
    void VectorMult(double ans[4], double T[4][4], double P[4]);

    //Function for matrix-matrix multiplication (T (4x4) times T (4x4))
    void MatrixMult(double ans[4][4],double T_L[4][4],double T_R[4][4]);

    //Function to calculate the inverse of a T (4x4) Matrix
    void InverseTMatrix(double ans[4][4],double T[4][4]);
};
#endif //ROBOT_MANIPULATOR__MATRIXMATH_H
