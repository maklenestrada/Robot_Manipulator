#ifndef ROBOT_MANIPULATOR__MATHOPERATIONS_H
#define ROBOT_MANIPULATOR__MATHOPERATIONS_H
class MathOperations {
public:
    //Constructor
    MathOperations() {}

    //Function for matrix-vector multiplication (T (4x4) times P (4x1))
    void VectorMult(double ans[4], double T[4][4], double P[4]);

    //Function for matrix-matrix multiplication (T (4x4) times T (4x4))
    void MatrixMult(double ans[4][4],double T_L[4][4],double T_R[4][4]);

    //Function to calculate the inverse of a T (4x4) Matrix
    void InverseTMatrix(double ans[4][4],double T[4][4]);

    //Function for cross product
    void CrossProduct(double ans[3],double X[3],double Y[3]);

    //Function for dot product
    void DotProduct(double* ans,double X[3],double Y[3]);

    //Function to calculate norm of vector
    void Norm(double* ans,double vec[3]);

    //Function to multiply vector times scalar
    void MultVecByScalar(double ans[3], double vec[3], double c);

    //Function to check if value is near a desired tolerance
    int ValueNear(double val, double goal, double tol);
};
#endif //ROBOT_MANIPULATOR__MATHOPERATIONS_H
