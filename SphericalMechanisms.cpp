#include "SphericalMechanisms.h"
#include "MathOperations.h"
#include "MathOperationsTemp.h"
#include <math.h>
#include <complex>

#define D2R M_PI/180.0
#define R2D 180/M_PI

using namespace std;
//Function that uses trig approach to solve eq in the form Ac + Bs + D = 0
int SphericalMechanisms::SolveTrig(double A, double B, double D, double* ang_a, double* ang_b)
{
     double normAB;
     normAB = sqrt(A*A + B*B);

     double c_gamma, s_gamma;
     c_gamma = A/normAB;
     s_gamma = B/normAB;

     double gamma;
     gamma = atan2(s_gamma,c_gamma)*R2D;

     double th_min_gamm1 = acos(-D/normAB)*R2D ;
     double th_min_gamm2 = 360 - th_min_gamm1;

     *ang_a = th_min_gamm1 + gamma;
     //Check to see that (th_min_gamm + gamm) is not great than 360, if so subtract by 360
     if(th_min_gamm2 + gamma >= 360){
         *ang_b = th_min_gamm2 + gamma - 360;
     }else {
         *ang_b = th_min_gamm2 + gamma;
     }


}

//Function to take the inner roots of two bi-quadratic equations and give corresponding outside root value
void SphericalMechanisms::CalcOutsideRootFromInner(double a1,double b1,double d1,double e1,double f1,double g1,
                                                        double h1,double i1,double j1,double a2,double b2,double d2,
                                                        double e2,double f2,double g2,double h2,double i2,double j2,
                                                        int d,double x_r[],double x_i[],double* x2_real,double* x2_imag)
{
    double x2_r;
    double x2_i;

    for(int i = 0; i < d; i++)
    {
        CalcSylvesterOutsideRoot(a1, b1, d1, e1, f1, g1, h1, i1, j1,
                                 a2, b2, d2, e2, f2, g2, h2, i2, j2, x_r[i], x_i[i],
                                 &x2_r, &x2_i);

        x2_real[i] = x2_r;
        x2_imag[i] = x2_i;
    }
}

//Function to find outside root using Sylvester's solution to two bi-quadratic equations
void SphericalMechanisms::CalcSylvesterOutsideRoot(double a1,double b1,double d1,double e1,double f1,double g1,
                                                   double h1,double i1,double j1,double a2,double b2,double d2,
                                                   double e2,double f2,double g2,double h2,double i2,double j2,
                                                   double x_r,double x_i,double* x2_r,double* x2_i)
{
    MathOperations MathOps;

    //Combine x to have real and imag part
    complex<double> x(x_r, x_i);

    // Calculate L1, M1, N1, L2, M2, N2 as complex numbers
    complex<double> L1 = a1 * (x * x) + b1 * x + d1;
    complex<double> M1 = e1 * (x * x) + f1 * x + g1;
    complex<double> N1 = h1 * (x * x) + i1 * x + j1;
    complex<double> L2 = a2 * (x * x) + b2 * x + d2;
    complex<double> M2 = e2 * (x * x) + f2 * x + g2;
    complex<double> N2 = h2 * (x * x) + i2 * x + j2;

    //Solve for vector [x2^3,x2^2,x2] = (Mat)^-1 * Vec
    complex<double> Mat[3][3] = {
            {0, L1, M1},
            {0, L2, M2},
            {L1, M1, N1}
    };
    complex<double> Vec[3] = {-N1, -N2, 0};
    complex<double> Inv[3][3];
    complex<double> ans[3];
    MatrixInv(Inv, Mat);
    MatrixVecMult(ans, Inv, Vec);

    // Extract the result for x2 as a complex number
    complex<double> x2 = ans[2]; //The third comp of the ans vec gives the x2 value corresponding to x1 root

    //Extract the real and imaginary parts of x2
    *x2_r = x2.real();
    *x2_i = x2.imag();
}