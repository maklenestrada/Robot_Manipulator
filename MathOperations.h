#ifndef ROBOT_MANIPULATOR__MATHOPERATIONS_H
#define ROBOT_MANIPULATOR__MATHOPERATIONS_H
class MathOperations {
public:
    //Constructor
    MathOperations() {}

    //Function to root polynomial
    int Poly_Solve(double root_r[], double root_c[], int d, double xcof[]);

    //Function to solve two bi-quadratic equations
    void SolveBiQuadratic(double a1,double b1,double d1,double e1,double f1,double g1,double h1,double i1,
                           double j1,double a2,double b2,double d2,double e2,double f2,double g2,double h2,
                           double i2,double j2,double* C8,double* C7,double* C6,double* C5,double* C4,
                           double* C3,double* C2,double* C1,double* C0);

    //Function to take the determinant of a 3x3 matrix
    double MatrixDet_R3(double Mat[3][3]);

    //Function to take the inverse of a 3x3 matrix
    void MatrixInv_R3(double inv[3][3], double Mat[3][3]);

    //Function to multiply 3x3 matrix by 3x1 vector
    void MultMatVec_R3(double ans[3], double Mat[3][3], double Vec[3]);

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

    //Function to add to vectors
    void AddVectors(double ans[3], double x[3],double y[3]);

    //Function to check if value is near a desired tolerance
    int ValueNear(double val, double goal, double tol);
};
#endif //ROBOT_MANIPULATOR__MATHOPERATIONS_H
