#ifndef ROBOT_MANIPULATOR_SPHERICALMECHANISMS_H
#define ROBOT_MANIPULATOR_SPHERICALMECHANISMS_H
#include <iostream>

class SphericalMechanisms{
    public:
    //Constructor
    SphericalMechanisms() {}

    //Function that uses trig approach to solve eq in the form Ac + Bs + D = 0
    int SolveTrig(double A, double B, double D, double *ang_a, double *ang_b);

};
#endif //ROBOT_MANIPULATOR_SPHERICALMECHANISMS_H
