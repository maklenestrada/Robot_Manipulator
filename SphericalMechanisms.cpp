#include "SphericalMechanisms.h"
#include "MathOperations.h"
#include <math.h>

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
//        cout << x_r[i] << endl;
//        cout << x_i[i] << endl;
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
//    //Solve using real part of root
//    double L1_r = a1 * (x_r * x_r) + b1 * x_r + d1;
//    double M1_r = e1 * (x_r * x_r) + f1 * x_r + g1;
//    double N1_r = h1 * (x_r * x_r) + i1 * x_r + j1;
//    double L2_r = a2 * (x_r * x_r) + b2 * x_r + d2;
//    double M2_r = e2 * (x_r * x_r) + f2 * x_r + g2;
//    double N2_r = h2 * (x_r * x_r) + i2 * x_r + j2;
//
//    //Solve using imag part of root
//    double L1_i = a1 * (x_i * x_i) + b1 * x_i + d1;
//    double M1_i = e1 * (x_i * x_i) + f1 * x_i + g1;
//    double N1_i = h1 * (x_i * x_i) + i1 * x_i + j1;
//    double L2_i = a2 * (x_i * x_i) + b2 * x_i + d2;
//    double M2_i = e2 * (x_i * x_i) + f2 * x_i + g2;
//    double N2_i = h2 * (x_i * x_i) + i2 * x_i + j2;

    MathOperations MathOps;
    //Solve for vector [x2^3,x2^2,x2] = (Mat)^-1 * Vec
    //Calculating real root of x2
    double Mat_r[3][3] = {
            {0,  L1_r, M1_r},
            {0,  L2_r, M2_r},
            {L1_r, M1_r, N1_r}};
    double Vec_r[3] = {-N1_r, -N2_r, 0};
    double Inv_r[3][3] = {0};
    double ans_r[3] = {0};
    MathOps.MatrixInv_R3(Inv_r, Mat_r);
    MathOps.MultMatVec_R3(ans_r, Mat_r, Vec_r); //The third comp of this vec gives the x2_r value corresponding to x1_r root
    *x2_r = ans_r[2];
    cout << ans_r[2] << endl;
    //Calculating imag root of x2
    double Mat_i[3][3] = {
            {0,  L1_i, M1_i},
            {0,  L2_i, M2_i},
            {L1_i, M1_i, N1_i}};
    double Vec_i[3] = {-N1_i, -N2_i, 0};
    double Inv_i[3][3] = {0};
    double ans_i[3] = {0};
    MathOps.MatrixInv_R3(Inv_i, Mat_i);
    MathOps.MultMatVec_R3(ans_i, Mat_i, Vec_i);//The third comp of this vec gives the x2_i value corresponding to x1_i root
    *x2_i = ans_i[2];

}