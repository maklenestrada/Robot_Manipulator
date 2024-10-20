#include "SphericalMechanisms.h"
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
