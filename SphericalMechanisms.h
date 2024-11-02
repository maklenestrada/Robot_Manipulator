#ifndef ROBOT_MANIPULATOR_SPHERICALMECHANISMS_H
#define ROBOT_MANIPULATOR_SPHERICALMECHANISMS_H
#include <iostream>

class SphericalMechanisms{
    public:
    //Constructor
    SphericalMechanisms() {}

    //Function that uses trig approach to solve eq in the form Ac + Bs + D = 0
    int SolveTrig(double A, double B, double D, double *ang_a, double *ang_b);

    //Function to take the inner roots of two bi-quadratic equations and give corresponding outside root value
    void CalcOutsideRootFromInner(double a1,double b1,double d1,double e1,double f1,double g1,
                                  double h1,double i1,double j1,double a2,double b2,double d2,
                                  double e2,double f2,double g2,double h2,double i2,double j2,
                                  int d,double x_r[],double x_i[],double* x2_real,double* x2_imag);

    //Function to find outside root using Sylvester's solution to two bi-quadratic equations
    void CalcSylvesterOutsideRoot(double a1,double b1,double d1,double e1,double f1,double g1,
                                  double h1,double i1,double j1,double a2,double b2,double d2,
                                  double e2,double f2,double g2,double h2,double i2,double j2,
                                  double x_r,double x_i,double* x2_r,double* x2_i);
};
#endif //ROBOT_MANIPULATOR_SPHERICALMECHANISMS_H
