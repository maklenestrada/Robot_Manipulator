#include "ReverseAnalysis.h"
#include <iostream>
#include "RobotKinematics.h"
#include "MathOperations.h"
#include "SphericalMechanisms.h"
#include <math.h>
#include <cmath>

using namespace std;
#define D2R M_PI/180.0
#define R2D 180/M_PI
void ReverseAnalysis::ReverseAnalysis_GEP60(double P_tool_6[3],double P_tool_F[3],double S6_F[3],double a67_F[3],double S6,
                                            double* valid,double th[6][8]) {
    //Defining GE P60 Robot Constant Mechanism Parameters
    //Twist Angles
    double alpha_12,alpha_23,alpha_34,alpha_45,alpha_56;
    alpha_12 = 270*D2R;
    alpha_23 = 0;
    alpha_34 = 0;
    alpha_45 = 270*D2R;
    alpha_56 = 90*D2R;
    //Joint Offsets
    double S2,S3,S4,S5;
    S2 = 0; //cm
    S3 = 0; //cm
    S4 = 9.8; //cm
    S5 = 14.5; //cm
    //Link Length (cm)
    double a12, a23, a34, a45, a56;
    a12 = 0;
    a23 = 70;
    a34 = 90;
    a45 = 0;
    a56 = 0;
    double a67, alpha_67;
    a67 = 0;
    alpha_67 = 90*D2R;

    RobotKinematics Robot_GEP60(alpha_12,alpha_23,alpha_34,alpha_45,alpha_56,alpha_67,S2,S3,S4,S5,a12,a23,a34,a45,a56,a67);

    //Closed Loop Analysis
    double a71 = 0;
    double S7 = 0;
    double S1 = 0;
    double alpha_71 = 0;
    double th7 = 0;
    double gamma1 = 0;
    Robot_GEP60.Closed_Loop(P_tool_6, P_tool_F, S6_F,a67_F,
                    &a71, &S7, &S1, &alpha_71, &th7, &gamma1);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    Using the closed_Loop function I get values for theta 7 to use for the analysis
     The next step is to solve for theta 1
     */

    for(int i = 0; i < 8; i++)
    {
        th[0][i] = th7 * D2R;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Equation 11.55
    // c1 [-S6 * Y7 + S7 * s71] + s1 [ - S6 * X7 - a71] + S4 = 0
    // Ac1 + Bs1 + D = 0
    SphericalMechanisms SphMech;

    double Y7 = - (sin(alpha_71)*cos(alpha_67) + cos(alpha_71)*sin(alpha_67)*cos(th7));
    double X7 = sin(alpha_67)*sin(th7);
    double A = -S6*Y7 + S7*sin(alpha_71);
    double B =  - S6*X7 - a71;
    double D = S4;
    double th1A, th1B;

    SphMech.SolveTrig(A,B,D, &th1A, &th1B);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    To solve for th1 the vector loop equation is dotted with S2, and set 14 of the direction cosine tables is used to
    get an equation that contains th1 as the only unkown
    //Equation 11.55
    // c1 [-S6 * Y7 + S7 * s71] + s1 [ - S6 * X7 - a71] + S4 = 0
     The next step is to solve for theta 5
     */

    for(int i = 0; i < 4; i++)
    {
        th[1][i] = th1A * D2R;
    }
    for(int i = 4; i < 8; i++)
    {
        th[1][i] = th1B * D2R;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Z17 = Z5
    // Z17 = c5

    double Xb1 = 0;
    double Yb1 = 0;
    double Zb1 = 0;
    double Z17 = 0;
    double c5[8];
    for(int i = 0; i < 8; i++){
        Xb1 = sin(alpha_12) * sin(th[1][i]);
        Yb1 = - (sin(alpha_71)*cos(alpha_12) + cos(alpha_71)*sin(alpha_12)*cos(th[1][i]));
        Zb1 = cos(alpha_71)*cos(alpha_12) - sin(alpha_71)*sin(alpha_12)*cos(th[1][i]);
        Z17 = sin(alpha_67) * (Xb1 * sin(th[0][i]) + Yb1 * cos(th7)) + cos(alpha_67) * Zb1;
        if(i == 0 || i == 1 || i == 4 || i == 5){
            th[2][i] = acos(Z17);
        } else {
            th[2][i] = 2*M_PI - acos(Z17);
        }

    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    To solve for th5 I used acos(Z17) then found the other angle for cos by 2PI - th
     Equations Used:
     (aij = alpha_ij)
    Xb1 = sin(a12) * sin(th1);
    Yb1 = - ( sin(a71) * cos(a12) + cos(a71) * sin(a12) * cos(th1));
    Zb1 = cos(a71) * cos(a12) - sin(a71) * sin(a12) * cos(th1)
    Z17 = sin(a67) * (Xb1 * sin(th7) + Yb1 * cos(th7)) + cos(a67) * Zb1
    The next step is to solve for theta 6
     */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double X17 = 0;
    double Y17 = 0;
    double c6 = 0;
    double s6 = 0;
    Xb1 = 0;
    Yb1 = 0;
    Zb1 = 0;

    for(int i = 0; i < 8; i++)
    {
        Xb1 = sin(alpha_12) * sin(th[1][i]);
        Yb1 = - (sin(alpha_71)*cos(alpha_12) + cos(alpha_71)*sin(alpha_12)*cos(th[1][i]));
        Zb1 = cos(alpha_71)*cos(alpha_12) - sin(alpha_71)*sin(alpha_12)*cos(th[1][i]);
        X17 = Xb1*cos(th[0][i]) - Yb1*sin(th[0][i]);
        Y17 = cos(alpha_67) * (Xb1*sin(th[0][i]) + Yb1*cos(th[0][i])) - sin(alpha_67)*Zb1;

        c6 = X17/sin(th[2][i]);
        s6 = Y17/sin(th[2][i]);

        th[3][i] = atan2(s6,c6);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    To solve for th6 I X17 = X56 & Y17 = -Xstar56 is used to get c6 and s6
     Equations Used:
     (aij = alpha_ij)
    Xb1 = sin(a12) * sin(th1);
    Yb1 = - ( sin(a71) * cos(a12) + cos(a71) * sin(a12) * cos(th1));
    Zb1 = cos(a71) * cos(a12) - sin(a71) * sin(a12) * cos(th1)
    X17 = Xb1 * cos(th7) - Yb1 * sin(th7)
    Y17 = cos(a67) * (Xb1 * sin(th7) + Yb1 * cos(th7)) - sin(a67) * Zb1
    c6 = X17/sin(th5)
    s6 = Y17/sin(th5)
    th6 = atan2(s71,c71)
    The next step is to solve for theta 3
     */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    double X6 = 0;double Y6 = 0;double Z6 = 0;
    double X67 = 0; double Y67 = 0; double Z67 = 0;
    double X671 = 0; double Y671 = 0;
    double X71 = 0; double Y71 = 0;
    double X1 = 0; double Y1 = 0;
    X7 = 0; Y7 = 0;double Z7 = 0;
    double K1 = 0; double K2 = 0;
    double c3 = 0;
    for(int i = 0; i < 8; i++)
    {
        X6 = sin(alpha_56)*sin(th[3][i]);
        Y6 = -(sin(alpha_67)*cos(alpha_56) + cos(alpha_67)*sin(alpha_56)*cos(th[3][i]));
        Z6 = cos(alpha_67)*cos(alpha_56) - sin(alpha_67)*sin(alpha_56)*cos(th[3][i]);

        X67 = X6*cos(th[0][i]) - Y6*sin(th[0][i]);
        Y67 = cos(alpha_71) * (X6*sin(th[0][i]) + Y6*cos(th[0][i])) - sin(alpha_71)*Z6;
        Z67 = sin(alpha_71) * (X6*sin(th[0][i]) + Y6*cos(th[0][i])) +cos(alpha_71)*Z6;

        X671 = X67*cos(th[1][i]) - Y67*sin(th[1][i]);
        Y671 = cos(alpha_12) * (X67*sin(th[1][i]) + Y67*cos(th[1][i])) - sin(alpha_12)*Z67;

        X7 = sin(alpha_67) * sin(th[0][i]);
        Y7 = - (sin(alpha_71) * cos(alpha_67) + cos(alpha_71) * sin(alpha_67) * cos(th[0][i]));
        Z7 = cos(alpha_71)*cos(alpha_67) - sin(alpha_71)*sin(alpha_67)*cos(th[0][i]);

        X71 = X7*cos(th[1][i]) - Y7*sin(th[1][i]);
        Y71 = cos(alpha_12) * (X7*sin(th[1][i]) + Y7*cos(th[1][i])) - sin(alpha_12)*Z7;

        X1 = sin(alpha_71)*sin(th[1][i]);
        Y1 = -(sin(alpha_12)*cos(alpha_71) + cos(alpha_12)*sin(alpha_71)*cos(th[1][i]));

        K1 = -S5*X671 - S6*X71 - S7*X1 -a71*cos(th[1][i]);
        K2 = -S1 - S5*Y671 - S6*Y71 - S7*Y1;

        c3 = (K1*K1 + K2*K2 - a23*a23 - a34*a34)/(2 * a23 * a34);
        if(i == 0 || i == 2 || i == 4 || i == 6){
            th[4][i] = acos(c3);
        } else {
            th[4][i] = 2*M_PI - acos(c3);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    To solve for th3 I used these equations
     (aij = alpha_ij)
    K1 = S5*X671 - S6*X71 - S7*X1 -a71*cos(th1)
    K2 = -S1 - S5*Y671 - S6*Y71 - S7*Y1
    cos(th3) = (K1*K1 + K2*K2 - a23*a23 - a34*a34)/(2 * a23 * a34)
    X671 = X67*cos(th1) - Y67*sin(th1)
    Y671 = cos(a12) * (X67*sin(th1) + Y67*cos(th1)) - sin(a12)*Z67

    X67 = X6*cos(th7) - Y6*sin(th7)
    Y67 = cos(a71) * (X6*sin(th7) + Y6*cos(th7)) - sin(a71)*Z6
    Z67 = sin(a71) * (X6*sin(th7) + Y6*cos(th7)) +cos(a71)*Z6

    X6 = sin(a56)*sin(th6)
    Y6 = -(sin(a67)*cos(a56) + cos(a67)*sin(a56)*cos(th6))
    Z6 = cos(a67)*cos(a56) - sin(a67)*sin(a56)*cos(th6)

    X71 = X7*cos(th1) - Y7*sin(th1)
    Y71 = cos(a12) * (X7*sin(th1) + Y7*cos(th1)) - sin(a12)*Z7

    X1 = sin(a71)*sin(th1)
    Y1 = -(sin(a12)*cos(a71) cos(a12)*sin(a71)*cos(th1))

    X7 = sin(a67) * sin(th7)
    Y7 = - (sin(a71) * cos(a67) + cos(a71) * sin(a67) * cos(th7))
    Z7 = cos(a71)*cos(a67) - sin(a71)*sin(a67)*cos(th7)

    The next step is to solve for theta 2
     */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MathOperations MathOps;

    double Amat[2][2] = {0};
    double b[2] = {0};
    double x[2] = {0};
    X6 = 0;Y6 = 0;Z6 = 0;
    X67 = 0;Y67 = 0;Z67 = 0;
    X671 = 0;Y671 = 0;
    X7 = 0;Y7 = 0;Z7 = 0;
    X71 = 0;Y71 = 0;
    X1 = 0; Y1 = 0;
    K1 = 0;K2 = 0;
    double c2 = 0;double s2 = 0;

    for(int i = 0; i < 8; i++) {
        //Set up A matrix
        Amat[0][0] = a23 + a34 * cos(th[4][i]);
        Amat[0][1] = -a34 * sin(th[4][i]);
        Amat[1][0] = -a34 * sin(th[4][i]);
        Amat[1][1] = -a23 - a34 * cos(th[4][i]);

        //Solve for K1 and K2
        X6 = sin(alpha_56)*sin(th[3][i]);
        Y6 = -(sin(alpha_67)*cos(alpha_56) + cos(alpha_67)*sin(alpha_56)*cos(th[3][i]));
        Z6 = cos(alpha_67)*cos(alpha_56) - sin(alpha_67)*sin(alpha_56)*cos(th[3][i]);

        X67 = X6*cos(th[0][i]) - Y6*sin(th[0][i]);
        Y67 = cos(alpha_71) * (X6*sin(th[0][i]) + Y6*cos(th[0][i])) - sin(alpha_71)*Z6;
        Z67 = sin(alpha_71) * (X6*sin(th[0][i]) + Y6*cos(th[0][i])) +cos(alpha_71)*Z6;

        X671 = X67*cos(th[1][i]) - Y67*sin(th[1][i]);
        Y671 = cos(alpha_12) * (X67*sin(th[1][i]) + Y67*cos(th[1][i])) - sin(alpha_12)*Z67;

        X7 = sin(alpha_67) * sin(th[0][i]);
        Y7 = - (sin(alpha_71) * cos(alpha_67) + cos(alpha_71) * sin(alpha_67) * cos(th[0][i]));
        Z7 = cos(alpha_71)*cos(alpha_67) - sin(alpha_71)*sin(alpha_67)*cos(th[0][i]);

        X71 = X7*cos(th[1][i]) - Y7*sin(th[1][i]);
        Y71 = cos(alpha_12) * (X7*sin(th[1][i]) + Y7*cos(th[1][i])) - sin(alpha_12)*Z7;

        X1 = sin(alpha_71)*sin(th[1][i]);
        Y1 = -(sin(alpha_12)*cos(alpha_71) + cos(alpha_12)*sin(alpha_71)*cos(th[1][i]));

        K1 = -S5*X671 - S6*X71 - S7*X1 - a71*cos(th[1][i]);
        K2 = -S1 - S5*Y671 - S6*Y71 - S7*Y1;

        b[0] = K1;
        b[1] = K2;

        //Solve system of equations for c2 and s2
        MathOps.solveLinearSystem(Amat, b, &c2, &s2);

        th[5][i] = atan2(s2,c2);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    To solve for th2 I solved two linear equations
     (aij = alpha_ij)
     c2 * (a23 + a34*c3) + s2 * (-a34s3) = K1
     c2 * (-a34*S3) + s2 * (-a23 - a34*c3) = K2

    The next step is to solve for theta 4
     */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    double X67123 = 0; double Y67123 = 0;
    double X6712 = 0; double Y6712 = 0; double Z6712 = 0;
    X671 = 0; Y671 = 0; double Z671 = 0;
    X67 = 0; Y67 = 0; Z67 = 0;
    X6 = 0; Y6 = 0; Z6 = 0;
    double c4 = 0; double s4 = 0;

    for(int i = 0; i < 8; i++)
    {
        X6 = sin(alpha_56)*sin(th[3][i]);
        Y6 = -(sin(alpha_67)*cos(alpha_56) + cos(alpha_67)*sin(alpha_56)*cos(th[3][i]));
        Z6 = cos(alpha_67)*cos(alpha_56) - sin(alpha_67)*sin(alpha_56)*cos(th[3][i]);

        X67 = X6*cos(th[0][i]) - Y6*sin(th[0][i]);
        Y67 = cos(alpha_71) * (X6*sin(th[0][i]) + Y6*cos(th[0][i])) - sin(alpha_71)*Z6;
        Z67 = sin(alpha_71) * (X6*sin(th[0][i]) + Y6*cos(th[0][i])) +cos(alpha_71)*Z6;

        X671 = X67*cos(th[1][i]) - Y67*sin(th[1][i]);
        Y671 = cos(alpha_12) * (X67*sin(th[1][i]) + Y67*cos(th[1][i])) - sin(alpha_12)*Z67;
        Z671 = sin(alpha_12) * (X67*sin(th[1][i]) + Y67*cos(th[1][i])) + cos(alpha_12)*Z67;

        X6712 = X671*cos(th[5][i]) - Y671*sin(th[5][i]);
        Y6712 = cos(alpha_23) * (X671*sin(th[5][i]) + Y671*cos(th[5][i])) - sin(alpha_23)*Z671;
        Z6712 = sin(alpha_23) * (X671*sin(th[5][i]) + Y671*cos(th[5][i])) + cos(alpha_23)*Z671;

        X67123 = X6712*cos(th[4][i]) - Y6712*sin(th[4][i]);
        Y67123 = cos(alpha_34) * (X6712*sin(th[4][i]) + Y6712*cos(th[4][i])) - sin(alpha_34)*Z6712;

        c4 = -Y67123;
        s4 = -X67123;

        th[6][i] = atan2(s4,c4);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    To solve for th4 I used these equations
     (aij = alpha_ij)
    c4 = -Y67123
    s4 = -X67123
    (aij = alpha_ij)
    X67123 = X6712*cos(th3) - Y6712*sin(th3)
    Y67123 = cos(a34) * (X6712*sin(th3) + Y6712*cos(th3)) - sin(a34)*Z6712

    X6712 = X671*cos(th2) - Y671*sin(th2)
    Y6712 = cos(a23) * (X671*sin(th2) + Y671*cos(th2)) - sin(a23)*Z671
    Z6712 = sin(a23) * (X671*sin(th2) + Y671*cos(th2)) + cos(a23)*Z671

    X671 = X67*cos(th1) - Y67*sin(th1)
    Y671 = cos(a12) * (X67*sin(th1) + Y67*cos(th1)) - sin(a12)*Z67
    Z671 = sin(a12) * (X67*sin(th1) + Y67*cos(th1)) + cos(a12)*Z67

    X67 = X6*cos(th7) - Y6*sin(th7)
    Y67 = cos(a71) * (X6*sin(th7) + Y6*cos(th7)) - sin(a71)*Z6
    Z67 = sin(a71) * (X6*sin(th7) + Y6*cos(th7)) +cos(a71)*Z6

    X6 = sin(a56)*sin(th6)
    Y6 = -(sin(a67)*cos(a56) + cos(a67)*sin(a56)*cos(th6))
    Z6 = cos(a67)*cos(a56) - sin(a67)*sin(a56)*cos(th6)
     */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    for (int i = 0; i < 7; ++i) {
//        for (int j = 0; j < 8; ++j) {
//            // Print each value with 4 decimal places
//            cout <<  th[i][j]*R2D << "\t";
//        }
//        cout << "\n"; // Move to the next row after each row of elements
//    }

}



