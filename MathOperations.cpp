#include "MathOperations.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;
//For function that roots polynomial
#define Sqrt(x)  sqrt((double)(x))
#define Fabs(x)  fabs((double)(x))

/* This routine will evaluate the roots of a polynomial of
      degree "d" ("d" must be less than or equal to 36).
   "root_r" and "root_c" are the real and complex parts of the
      'd' solutions to the original equation.
   "xcof" is an array of coefficients, ordered from smallest
      to largest power.

   xcof[16] x^16 + xcof[15] x^15 + ... + xcof[1] x + xcof[0] = 0 */
int MathOperations::Poly_Solve(double root_r[], double root_c[], int d, double xcof[])
{
    double coef[37], dis, X, Y, Z[37], X0, XX[40],YY[40],
            U, V, dUx, dUy, den , dX, dY, dXY, XY, C, B[40] ;
    int i, k, deg, cnt ;
    int lst, lflip, ltry ;

    if (d > 36)
        return (0) ;

    for (i=0 ; i<=d ; ++i)
        coef[i] = xcof[i] ;

    deg = d ;

    while (coef[deg] == 0.0)
        deg-- ;  //The leading coefficient was zero.

    if (deg <1)
        return (-1) ;  //The polynomial must be at least of degree 1.

    cnt = 0 ;  //cnt keeps track of the number of roots found

    if (deg == 1)
        goto solve_linear ;

    if (deg == 2)
        goto solve_quad ;

    /**************************/
    /*  Set initial values.   */
    /**************************/
    L30:
    lst = 0 ;   //lst counts the number of different starting values
    lflip = 0 ; /*lflip determines whether the inverse polynomial is
                 being considered */

    X = 0.00608 ;
    Y = 0.2939 ;

    L35:
    X0 = X ;
    X  = -5.0*Y ;
    Y  =  2.0*X0 ;

    ltry = 0 ;  //ltry counts the # of interations for a starting value

    lst++ ;

    L38:
    XX[0] = 1.0 ;
    YY[0] = 0.0 ;

    for (i=1 ; i<=deg ; ++i)
        //Evaluate x^16, x^15, etc where x is complex
    {XX[i] = X * XX[i-1] - Y * YY[i-1] ;
        YY[i] = X * YY[i-1] + Y * XX[i-1] ;
    }
    U = coef[0] ;
    V = 0.0 ;

    for (i=1 ; i<=deg ; ++i) //Evaluate the polynomial.
    {U += coef[i] * XX[i] ;
        V += coef[i] * YY[i] ;
    }

    dUx = 0.0 ;
    dUy = 0.0 ;

    for (i=1 ; i<=deg ; ++i)
    {dUx += i*coef[i] * XX[i-1] ;
        dUy -= i*coef[i] * YY[i-1] ;
    }
    den = dUx*dUx + dUy*dUy ;

    dX =  (V*dUy - U*dUx)/den ;
    dY = -(U*dUy + V*dUx)/den ;

    X += dX ;   //Next try for root.
    Y += dY ;

    if (Fabs(X) < 40.0)
    {dXY = Sqrt(dX*dX + dY*dY) ;
        XY =  Sqrt(X*X + Y*Y) ;

        if (Fabs(dXY/XY) > 0.0000000002)   //was 0.0000001
        {ltry++ ;
            if (ltry<400)     //was 300
                goto L38 ;
            else
                goto flip_poly ;
        }
        else
            goto reduce_poly ;
    }
    flip_poly:
    lflip++ ;
    ltry = 0 ;

    for (k=0 ; k<=deg ; ++k)
        Z[k] = coef[deg-k] ;

    for (k=0 ; k<=deg ; ++k)
        coef[k] = Z[k] ;

    if (lflip ==1)
    {X = 0.189 ;
        Y = -0.132 ;
        goto L38 ;
    }

    if (lflip ==2)
        if (lst < 4 )
            goto L35 ;
    return (-300) ;  /*A solution was not found for 300 iterations
                         for 4 starting values. */

    /**************************/
    reduce_poly:

    if (Fabs(Y) < 0.000006)  //was 0.0000005
        Y = 0.0 ;
    cnt++ ;

    if (lflip ==1)
    {for (k=0 ; k<=deg ; ++k)  //flip it back
            Z[k] = coef[deg-k] ;
        for (k=0 ; k<=deg ; ++k)
            coef[k] = Z[k] ;

        den = X*X + Y*Y ; //The root to the orig. eqn is 1/(X+iY)
        root_r[cnt-1] = X = X/den ;
        root_c[cnt-1] = Y = Y/den ;
    }

    else
    {root_r[cnt-1] = X ;
        root_c[cnt-1] = Y ;
    }

    if (Y==0.0)
    {//Reduce the equation by one degree.

        C = X ;
        B[deg] = 0.0 ;
        for (k=deg-1 ; k>=0 ; --k)
            B[k] = coef[k+1] + C * B[k+1] ;  /*115*/

        deg-- ;  //Reduce the degree of the polynomial by 1

        for (k=0 ; k<=deg ; ++k)
            coef[k] = B[k] ;

        if (deg ==2)
            goto solve_quad ;
        else if (deg ==1)
            goto solve_linear ;

        else
            goto L30 ;
    }

    else
    {//Reduce the equation by the complex conjugates.
        cnt++ ;
        root_r[cnt-1] =  X ;
        root_c[cnt-1] = -Y ;

        B[deg-2] = coef[deg] ;
        B[deg-3] = coef[deg-1] + 2.0* X * B[deg-2] ;

        for (k=deg-4 ; k>=0 ; --k)
        {B[k] = coef[k+2]  - (X*X+Y*Y) * B[k+2] + 2.0 * X * B[k+1] ;
        }
        deg -= 2 ;

        for (k=0 ; k<=deg ; ++k)
            coef[k] = B[k] ;

        if (deg==2)
            goto solve_quad ;
        if (deg==1)
            goto solve_linear ;
        else
            goto L30 ;
    }


    /**************************/
    solve_quad:
    dis = coef[1]*coef[1] - 4.0*coef[2]*coef[0] ;

    X = -coef[1] / (2.0*coef[2]) ;

    if (dis>= 0.0)
    {Y = Sqrt(dis) / (2.0*coef[2]) ;
        root_r[cnt] = X+Y ;
        root_r[cnt+1] = X-Y ;
        root_c[cnt] = root_c[cnt+1] = 0.0 ;
    }

    else
    {Y = Sqrt(-dis)/ (2.0*coef[2]) ;
        root_r[cnt] =   root_r[cnt+1] = X ;
        root_c[cnt] = -(root_c[cnt+1] = Y) ;
    }
    return (1) ;

    solve_linear:
    root_r[cnt] = -coef[0] / coef[1] ;
    root_c[cnt] = 0.0 ;
    return (1) ;
}

//Solve two Bi-Quadratic Equations using Sylvester's Methods
void MathOperations::SolveBiQuadratic(double a1,double b1,double d1,double e1,double f1,double g1,double h1,double i1,
                                           double j1,double a2,double b2,double d2,double e2,double f2,double g2,double h2,
                                           double i2,double j2,double* C8,double* C7,double* C6,double* C5,double* C4,
                                           double* C3,double* C2,double* C1,double* C0)
{
    //These are the coefficients for the 8th order poly to solve the bi-quad equations
    *C8 = -a1 * a1 * h2 * h2 + 2 * a1 * a2 * h1 * h2 + a1 * e1 * e2 * h2 - a1 * e2 * e2 * h1
         - a2 * a2 * h1 * h1 - a2 * e1 * e1 * h2 + a2 * e1 * e2 * h1;

    *C7 = -2 * a1 * a1 * h2 * i2 + 2 * a1 * a2 * h1 * i2 + 2 * a1 * a2 * h2 * i1
         - 2 * a1 * b1 * h2 * h2 + 2 * a1 * b2 * h1 * h2 + a1 * e1 * e2 * i2 + a1 * e1 * f2 * h2
         - a1 * e2 * e2 * i1 + a1 * e2 * f1 * h2 - 2 * a1 * e2 * f2 * h1 - 2 * a2 * a2 * h1 * i1
         + 2 * a2 * b1 * h1 * h2 - 2 * a2 * b2 * h1 * h1 - a2 * e1 * e1 * i2 + a2 * e1 * e2 * i1
         - 2 * a2 * e1 * f1 * h2 + a2 * e1 * f2 * h1 + a2 * e2 * f1 * h1 + b1 * e1 * e2 * h2
         - b1 * e2 * e2 * h1 - b2 * e1 * e1 * h2 + b2 * e1 * e2 * h1;

    *C6 = -2 * a1 * a1 * h2 * j2 - a1 * a1 * i2 * i2 + 2 * a1 * a2 * h1 * j2
         + 2 * a1 * a2 * h2 * j1 + 2 * a1 * a2 * i1 * i2 - 4 * a1 * b1 * h2 * i2
         + 2 * a1 * b2 * h1 * i2 + 2 * a1 * b2 * h2 * i1 - 2 * a1 * d1 * h2 * h2
         + 2 * a1 * d2 * h1 * h2 + a1 * e1 * e2 * j2 + a1 * e1 * f2 * i2 + a1 * e1 * g2 * h2
         - a1 * e2 * e2 * j1 + a1 * e2 * f1 * i2 - 2 * a1 * e2 * f2 * i1 + a1 * e2 * g1 * h2
         - 2 * a1 * e2 * g2 * h1 + a1 * f1 * f2 * h2 - a1 * f2 * f2 * h1 - 2 * a2 * a2 * h1 * j1
         - a2 * a2 * i1 * i1 + 2 * a2 * b1 * h1 * i2 + 2 * a2 * b1 * h2 * i1
         - 4 * a2 * b2 * h1 * i1 + 2 * a2 * d1 * h1 * h2 - 2 * a2 * d2 * h1 * h1
         - a2 * e1 * e1 * j2 + a2 * e1 * e2 * j1 - 2 * a2 * e1 * f1 * i2 + a2 * e1 * f2 * i1
         - 2 * a2 * e1 * g1 * h2 + a2 * e1 * g2 * h1 + a2 * e2 * f1 * i1 + a2 * e2 * g1 * h1
         - a2 * f1 * f1 * h2 + a2 * f1 * f2 * h1 - b1 * b1 * h2 * h2 + 2 * b1 * b2 * h1 * h2
         + b1 * e1 * e2 * i2 + b1 * e1 * f2 * h2 - b1 * e2 * e2 * i1 + b1 * e2 * f1 * h2
         - 2 * b1 * e2 * f2 * h1 - b2 * b2 * h1 * h1 - b2 * e1 * e1 * i2 + b2 * e1 * e2 * i1
         - 2 * b2 * e1 * f1 * h2 + b2 * e1 * f2 * h1 + b2 * e2 * f1 * h1 + d1 * e1 * e2 * h2
         - d1 * e2 * e2 * h1 - d2 * e1 * e1 * h2 + d2 * e1 * e2 * h1;

    *C5 = -2 * a1 * a1 * i2 * j2 + 2 * a1 * a2 * i1 * j2 + 2 * a1 * a2 * i2 * j1
         - 4 * a1 * b1 * h2 * j2 - 2 * a1 * b1 * i2 * i2 + 2 * a1 * b2 * h1 * j2
         + 2 * a1 * b2 * h2 * j1 + 2 * a1 * b2 * i1 * i2 - 4 * a1 * d1 * h2 * i2
         + 2 * a1 * d2 * h1 * i2 + 2 * a1 * d2 * h2 * i1 + a1 * e1 * f2 * j2 + a1 * e1 * g2 * i2
         + a1 * e2 * f1 * j2 - 2 * a1 * e2 * f2 * j1 + a1 * e2 * g1 * i2 - 2 * a1 * e2 * g2 * i1
         + a1 * f1 * f2 * i2 + a1 * f1 * g2 * h2 - a1 * f2 * f2 * i1 + a1 * f2 * g1 * h2
         - 2 * a1 * f2 * g2 * h1 - 2 * a2 * a2 * i1 * j1 + 2 * a2 * b1 * h1 * j2
         + 2 * a2 * b1 * h2 * j1 + 2 * a2 * b1 * i1 * i2 - 4 * a2 * b2 * h1 * j1
         - 2 * a2 * b2 * i1 * i1 + 2 * a2 * d1 * h1 * i2 + 2 * a2 * d1 * h2 * i1
         - 4 * a2 * d2 * h1 * i1 - 2 * a2 * e1 * f1 * j2 + a2 * e1 * f2 * j1
         - 2 * a2 * e1 * g1 * i2 + a2 * e1 * g2 * i1 + a2 * e2 * f1 * j1 + a2 * e2 * g1 * i1
         - a2 * f1 * f1 * i2 + a2 * f1 * f2 * i1 - 2 * a2 * f1 * g1 * h2 + a2 * f1 * g2 * h1
         + a2 * f2 * g1 * h1 - 2 * b1 * b1 * h2 * i2 + 2 * b1 * b2 * h1 * i2
         + 2 * b1 * b2 * h2 * i1 - 2 * b1 * d1 * h2 * h2 + 2 * b1 * d2 * h1 * h2
         + b1 * e1 * e2 * j2 + b1 * e1 * f2 * i2 + b1 * e1 * g2 * h2 - b1 * e2 * e2 * j1
         + b1 * e2 * f1 * i2 - 2 * b1 * e2 * f2 * i1 + b1 * e2 * g1 * h2 - 2 * b1 * e2 * g2 * h1
         + b1 * f1 * f2 * h2 - b1 * f2 * f2 * h1 - 2 * b2 * b2 * h1 * i1 + 2 * b2 * d1 * h1 * h2
         - 2 * b2 * d2 * h1 * h1 - b2 * e1 * e1 * j2 + b2 * e1 * e2 * j1 - 2 * b2 * e1 * f1 * i2
         + b2 * e1 * f2 * i1 - 2 * b2 * e1 * g1 * h2 + b2 * e1 * g2 * h1 + b2 * e2 * f1 * i1
         + b2 * e2 * g1 * h1 - b2 * f1 * f1 * h2 + b2 * f1 * f2 * h1 + d1 * e1 * e2 * i2
         + d1 * e1 * f2 * h2 - d1 * e2 * e2 * i1 + d1 * e2 * f1 * h2 - 2 * d1 * e2 * f2 * h1
         - d2 * e1 * e1 * i2 + d2 * e1 * e2 * i1 - 2 * d2 * e1 * f1 * h2 + d2 * e1 * f2 * h1
         + d2 * e2 * f1 * h1;

    *C4 = -a1 * a1 * j2 * j2 + 2 * a1 * a2 * j1 * j2 - 4 * a1 * b1 * i2 * j2
         + 2 * a1 * b2 * i1 * j2 + 2 * a1 * b2 * i2 * j1 - 4 * a1 * d1 * h2 * j2
         - 2 * a1 * d1 * i2 * i2 + 2 * a1 * d2 * h1 * j2 + 2 * a1 * d2 * h2 * j1
         + 2 * a1 * d2 * i1 * i2 + a1 * e1 * g2 * j2 + a1 * e2 * g1 * j2 - 2 * a1 * e2 * g2 * j1
         + a1 * f1 * f2 * j2 + a1 * f1 * g2 * i2 - a1 * f2 * f2 * j1 + a1 * f2 * g1 * i2
         - 2 * a1 * f2 * g2 * i1 + a1 * g1 * g2 * h2 - a1 * g2 * g2 * h1 - a2 * a2 * j1 * j1
         + 2 * a2 * b1 * i1 * j2 + 2 * a2 * b1 * i2 * j1 - 4 * a2 * b2 * i1 * j1
         + 2 * a2 * d1 * h1 * j2 + 2 * a2 * d1 * h2 * j1 + 2 * a2 * d1 * i1 * i2
         - 4 * a2 * d2 * h1 * j1 - 2 * a2 * d2 * i1 * i1 - 2 * a2 * e1 * g1 * j2
         + a2 * e1 * g2 * j1 + a2 * e2 * g1 * j1 - a2 * f1 * f1 * j2 + a2 * f1 * f2 * j1
         - 2 * a2 * f1 * g1 * i2 + a2 * f1 * g2 * i1 + a2 * f2 * g1 * i1 - a2 * g1 * g1 * h2
         + a2 * g1 * g2 * h1 - 2 * b1 * b1 * h2 * j2 - b1 * b1 * i2 * i2 + 2 * b1 * b2 * h1 * j2
         + 2 * b1 * b2 * h2 * j1 + 2 * b1 * b2 * i1 * i2 - 4 * b1 * d1 * h2 * i2
         + 2 * b1 * d2 * h1 * i2 + 2 * b1 * d2 * h2 * i1 + b1 * e1 * f2 * j2 + b1 * e1 * g2 * i2
         + b1 * e2 * f1 * j2 - 2 * b1 * e2 * f2 * j1 + b1 * e2 * g1 * i2 - 2 * b1 * e2 * g2 * i1
         + b1 * f1 * f2 * i2 + b1 * f1 * g2 * h2 - b1 * f2 * f2 * i1 + b1 * f2 * g1 * h2
         - 2 * b1 * f2 * g2 * h1 - 2 * b2 * b2 * h1 * j1 - b2 * b2 * i1 * i1
         + 2 * b2 * d1 * h1 * i2 + 2 * b2 * d1 * h2 * i1 - 4 * b2 * d2 * h1 * i1
         - 2 * b2 * e1 * f1 * j2 + b2 * e1 * f2 * j1 - 2 * b2 * e1 * g1 * i2 + b2 * e1 * g2 * i1
         + b2 * e2 * f1 * j1 + b2 * e2 * g1 * i1 - b2 * f1 * f1 * i2 + b2 * f1 * f2 * i1
         - 2 * b2 * f1 * g1 * h2 + b2 * f1 * g2 * h1 + b2 * f2 * g1 * h1 - d1 * d1 * h2 * h2
         + 2 * d1 * d2 * h1 * h2 + d1 * e1 * e2 * j2 + d1 * e1 * f2 * i2 + d1 * e1 * g2 * h2
         - d1 * e2 * e2 * j1 + d1 * e2 * f1 * i2 - 2 * d1 * e2 * f2 * i1 + d1 * e2 * g1 * h2
         - 2 * d1 * e2 * g2 * h1 + d1 * f1 * f2 * h2 - d1 * f2 * f2 * h1 - d2 * d2 * h1 * h1
         - d2 * e1 * e1 * j2 + d2 * e1 * e2 * j1 - 2 * d2 * e1 * f1 * i2 + d2 * e1 * f2 * i1
         - 2 * d2 * e1 * g1 * h2 + d2 * e1 * g2 * h1 + d2 * e2 * f1 * i1 + d2 * e2 * g1 * h1
         - d2 * f1 * f1 * h2 + d2 * f1 * f2 * h1;

    *C3 = -2 * a1 * b1 * j2 * j2 + 2 * a1 * b2 * j1 * j2 - 4 * a1 * d1 * i2 * j2
         + 2 * a1 * d2 * i1 * j2 + 2 * a1 * d2 * i2 * j1 + a1 * f1 * g2 * j2 + a1 * f2 * g1 * j2
         - 2 * a1 * f2 * g2 * j1 + a1 * g1 * g2 * i2 - a1 * g2 * g2 * i1 + 2 * a2 * b1 * j1 * j2
         - 2 * a2 * b2 * j1 * j1 + 2 * a2 * d1 * i1 * j2 + 2 * a2 * d1 * i2 * j1
         - 4 * a2 * d2 * i1 * j1 - 2 * a2 * f1 * g1 * j2 + a2 * f1 * g2 * j1 + a2 * f2 * g1 * j1
         - a2 * g1 * g1 * i2 + a2 * g1 * g2 * i1 - 2 * b1 * b1 * i2 * j2 + 2 * b1 * b2 * i1 * j2
         + 2 * b1 * b2 * i2 * j1 - 4 * b1 * d1 * h2 * j2 - 2 * b1 * d1 * i2 * i2
         + 2 * b1 * d2 * h1 * j2 + 2 * b1 * d2 * h2 * j1 + 2 * b1 * d2 * i1 * i2
         + b1 * e1 * g2 * j2 + b1 * e2 * g1 * j2 - 2 * b1 * e2 * g2 * j1 + b1 * f1 * f2 * j2
         + b1 * f1 * g2 * i2 - b1 * f2 * f2 * j1 + b1 * f2 * g1 * i2 - 2 * b1 * f2 * g2 * i1
         + b1 * g1 * g2 * h2 - b1 * g2 * g2 * h1 - 2 * b2 * b2 * i1 * j1 + 2 * b2 * d1 * h1 * j2
         + 2 * b2 * d1 * h2 * j1 + 2 * b2 * d1 * i1 * i2 - 4 * b2 * d2 * h1 * j1
         - 2 * b2 * d2 * i1 * i1 - 2 * b2 * e1 * g1 * j2 + b2 * e1 * g2 * j1 + b2 * e2 * g1 * j1
         - b2 * f1 * f1 * j2 + b2 * f1 * f2 * j1 - 2 * b2 * f1 * g1 * i2 + b2 * f1 * g2 * i1
         + b2 * f2 * g1 * i1 - b2 * g1 * g1 * h2 + b2 * g1 * g2 * h1 - 2 * d1 * d1 * h2 * i2
         + 2 * d1 * d2 * h1 * i2 + 2 * d1 * d2 * h2 * i1 + d1 * e1 * f2 * j2 + d1 * e1 * g2 * i2
         + d1 * e2 * f1 * j2 - 2 * d1 * e2 * f2 * j1 + d1 * e2 * g1 * i2 - 2 * d1 * e2 * g2 * i1
         + d1 * f1 * f2 * i2 + d1 * f1 * g2 * h2 - d1 * f2 * f2 * i1 + d1 * f2 * g1 * h2
         - 2 * d1 * f2 * g2 * h1 - 2 * d2 * d2 * h1 * i1 - 2 * d2 * e1 * f1 * j2
         + d2 * e1 * f2 * j1 - 2 * d2 * e1 * g1 * i2 + d2 * e1 * g2 * i1 + d2 * e2 * f1 * j1
         + d2 * e2 * g1 * i1 - d2 * f1 * f1 * i2 + d2 * f1 * f2 * i1 - 2 * d2 * f1 * g1 * h2
         + d2 * f1 * g2 * h1 + d2 * f2 * g1 * h1;

    *C2 = -2 * a1 * d1 * j2 * j2 + 2 * a1 * d2 * j1 * j2 + a1 * g1 * g2 * j2
         - a1 * g2 * g2 * j1 + 2 * a2 * d1 * j1 * j2 - 2 * a2 * d2 * j1 * j1
         - a2 * g1 * g1 * j2 + a2 * g1 * g2 * j1 - b1 * b1 * j2 * j2 + 2 * b1 * b2 * j1 * j2
         - 4 * b1 * d1 * i2 * j2 + 2 * b1 * d2 * i1 * j2 + 2 * b1 * d2 * i2 * j1
         + b1 * f1 * g2 * j2 + b1 * f2 * g1 * j2 - 2 * b1 * f2 * g2 * j1 + b1 * g1 * g2 * i2
         - b1 * g2 * g2 * i1 - b2 * b2 * j1 * j1 + 2 * b2 * d1 * i1 * j2 + 2 * b2 * d1 * i2 * j1
         - 4 * b2 * d2 * i1 * j1 - 2 * b2 * f1 * g1 * j2 + b2 * f1 * g2 * j1 + b2 * f2 * g1 * j1
         - b2 * g1 * g1 * i2 + b2 * g1 * g2 * i1 - 2 * d1 * d1 * h2 * j2
         - d1 * d1 * i2 * i2 + 2 * d1 * d2 * h1 * j2 + 2 * d1 * d2 * h2 * j1
         + 2 * d1 * d2 * i1 * i2 + d1 * e1 * g2 * j2 + d1 * e2 * g1 * j2 - 2 * d1 * e2 * g2 * j1
         + d1 * f1 * f2 * j2 + d1 * f1 * g2 * i2 - d1 * f2 * f2 * j1 + d1 * f2 * g1 * i2
         - 2 * d1 * f2 * g2 * i1 + d1 * g1 * g2 * h2 - d1 * g2 * g2 * h1 - 2 * d2 * d2 * h1 * j1
         - d2 * d2 * i1 * i1 - 2 * d2 * e1 * g1 * j2 + d2 * e1 * g2 * j1 + d2 * e2 * g1 * j1
         - d2 * f1 * f1 * j2 + d2 * f1 * f2 * j1 - 2 * d2 * f1 * g1 * i2 + d2 * f1 * g2 * i1
         + d2 * f2 * g1 * i1 - d2 * g1 * g1 * h2 + d2 * g1 * g2 * h1;

    *C1 = -2 * b1 * d1 * j2 * j2 + 2 * b1 * d2 * j1 * j2 + b1 * g1 * g2 * j2 - b1 * g2 * g2 * j1
         + 2 * b2 * d1 * j1 * j2 - 2 * b2 * d2 * j1 * j1 - b2 * g1 * g1 * j2 + b2 * g1 * g2 * j1
         - 2 * d1 * d1 * i2 * j2 + 2 * d1 * d2 * i1 * j2 + 2 * d1 * d2 * i2 * j1
         + d1 * f1 * g2 * j2 + d1 * f2 * g1 * j2 - 2 * d1 * f2 * g2 * j1 + d1 * g1 * g2 * i2
         - d1 * g2 * g2 * i1 - 2 * d2 * d2 * i1 * j1 - 2 * d2 * f1 * g1 * j2 + d2 * f1 * g2 * j1
         + d2 * f2 * g1 * j1 - d2 * g1 * g1 * i2 + d2 * g1 * g2 * i1;

    *C0 = -d1 * d1 * j2 * j2 + 2 * d1 * d2 * j1 * j2 + d1 * g1 * g2 * j2
         - d1 * g2 * g2 * j1 - d2 * d2 * j1 * j1 - d2 * g1 * g1 * j2 + d2 * g1 * g2 * j1;

}

//Function to take the determinant of a 3x3 matrix
double MathOperations::MatrixDet_R3(double Mat[3][3])
{
    return Mat[0][0] * (Mat[1][1] * Mat[2][2] - Mat[1][2] * Mat[2][1]) -
           Mat[0][1] * (Mat[1][0] * Mat[2][2] - Mat[1][2] * Mat[2][0]) +
           Mat[0][2] * (Mat[1][0] * Mat[2][1] - Mat[1][1] * Mat[2][0]);
}

//Function to calculate the inverse of  3x3 matrix
void MathOperations::MatrixInv_R3(double Inv[3][3], double Mat[3][3])
{
    double Det = MatrixDet_R3(Mat);  // Check if the matrix is linearly dependent
    if (Det == 0) {
        throw runtime_error("Matrix is singular and cannot be inverted.");
    }
    Inv[0][0] =  (Mat[1][1] * Mat[2][2] - Mat[1][2] * Mat[2][1]) / Det;
    Inv[0][1] = -(Mat[0][1] * Mat[2][2] - Mat[0][2] * Mat[2][1]) / Det;
    Inv[0][2] =  (Mat[0][1] * Mat[1][2] - Mat[0][2] * Mat[1][1]) / Det;

    Inv[1][0] = -(Mat[1][0] * Mat[2][2] - Mat[1][2] * Mat[2][0]) / Det;
    Inv[1][1] =  (Mat[0][0] * Mat[2][2] - Mat[0][2] * Mat[2][0]) / Det;
    Inv[1][2] = -(Mat[0][0] * Mat[1][2] - Mat[0][2] * Mat[1][0]) / Det;

    Inv[2][0] =  (Mat[1][0] * Mat[2][1] - Mat[1][1] * Mat[2][0]) / Det;
    Inv[2][1] = -(Mat[0][0] * Mat[2][1] - Mat[0][1] * Mat[2][0]) / Det;
    Inv[2][2] =  (Mat[0][0] * Mat[1][1] - Mat[0][1] * Mat[1][0]) / Det;
}

//Function to multiply 3x3 matrix by 3x1 vector
void MathOperations::MultMatVec_R3( double ans[3], double Mat[3][3], double Vec[3])
{
    for (int i = 0; i < 3; ++i) {
        ans[i] = 0; // Initialize each element of the result to zero
        for (int j = 0; j < 3; ++j) {
            ans[i] += Mat[i][j] * Vec[j]; // Matrix-vector multiplication
        }
    }
}

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

// Function to solve a 2x2 linear system
// Ax = B
// x = Ainv*B
void MathOperations::solveLinearSystem(double A[2][2], double B[2], double* x1, double* x2) {

    //Calculate Inverse of A
    double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    double Ainv[2][2];
    Ainv[0][0] = (1/det) * A[1][1];
    Ainv[0][1] = (1/det) * (-A[0][1]);
    Ainv[1][0] = (1/det) * (-A[1][0]);
    Ainv[1][1] = (1/det) * A[0][0];

    //Multiply Ainv times B (i.e. x = Ainv*B)

    double x[2] = {0};

    for(int i = 0; i < 2; i++){ //For loop over the T matrix rows
        for(int j = 0; j < 2; j++){//For loop over the P elements
            x[i] += Ainv[i][j] * B[j];
        }
    }

    *x1 = x[0];
    *x2 = x[1];

}