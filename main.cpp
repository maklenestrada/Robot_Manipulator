#include <iostream>
#include <iomanip>
#include "RobotKinematics.h"
#include "MathOperations.h"
#include "SphericalMechanisms.h"
#include "ReverseAnalysis.h"
#include <math.h>

#define D2R M_PI/180.0
#define R2D 180/M_PI

using namespace std;

//Homework 8 Reverse Analysis
int main(){
    //Defining GE P60 Robot Constant Mechanism Parameters
    //Twist Angles
    double alpha_12,alpha_23,alpha_34,alpha_45,alpha_56,alpha_67;
    alpha_12 = 270*D2R;
    alpha_23 = 0;
    alpha_34 = 0;
    alpha_45 = 270*D2R;
    alpha_56 = 90*D2R;
    alpha_67 = 90*D2R;
    //Joint Offsets
    double S2,S3,S4,S5;
    S2 = 0; //cm
    S3 = 0; //cm
    S4 = 9.8; //cm
    S5 = 14.5; //cm
    //Link Length (cm)
    double a12, a23, a34, a45, a56, a67;
    a12 = 0;
    a23 = 70;
    a34 = 90;
    a45 = 0;
    a56 = 0;
    a67 = 0;

    RobotKinematics Robot_GEP60(alpha_12,alpha_23,alpha_34,alpha_45,alpha_56,alpha_67,S2,S3,S4,S5,a12,a23,a34,a45,a56,a67);

    // Units in cm.
    double P_tool_6[3] = {2.0, 3.0, 5.0};
    double P_tool_F[3] = {80.0, 80.0, 18.0};
    double S6_F[3] = {-1.0, 1.0, 1.0};
    double a67_F[3] = {1.0, 2.0,  -1.0};
    double S6 = 15.24;
    double valid;  // Number of valid solutions.
    double th[7][8];
    ReverseAnalysis RevAnalysis;
    RevAnalysis.ReverseAnalysis_GEP60(P_tool_6, P_tool_F, S6_F, a67_F, S6,
                    &valid, th);
    double T_6toF[4][4] = {0};

    for (int i = 0; i < 8; ++i) {
        cout << "**********************" << endl;
        cout << "Solution #" << i + 1 << endl;
        cout << "phi1[" << i + 1 << "] = " << fixed << setprecision(3)
             << th[1][i] * R2D <<  "\u00B0" << endl;
        cout << "th2 [" << i + 1 << "] = " << fixed << setprecision(3)
             << th[5][i] * R2D <<  "\u00B0" << endl;
        cout << "th3 [" << i + 1 << "] = " << fixed << setprecision(3)
             << th[5][i] * R2D <<  "\u00B0"  << endl;
        cout << "th4 [" << i + 1 << "] = " << fixed << setprecision(3)
             << th[6][i] * R2D <<  "\u00B0" << endl;
        cout << "th5 [" << i + 1 << "] = " << fixed << setprecision(3)
             << th[2][i] * R2D <<  "\u00B0" << endl;
        cout << "th6 [" << i + 1 << "] = " << fixed << setprecision(3)
             << th[3][i] * R2D<<  "\u00B0"  << endl;
        cout << endl;


    Robot_GEP60.Forward(T_6toF, th[1][i],th[5][i],th[4][i],th[6][i],th[2][i],th[3][i],S6);

    MathOperations MatrixOp;
    //Getting the Point in Fixed Frame
    double P1_F[4] = {0};
    double P1_6[4] = {2.0, 3.0, 5.0, 1.0};
    MatrixOp.VectorMult(P1_F,T_6toF,P1_6);

        cout << "P_tool_F   = " << fixed << setprecision(2)
             << P1_F[0] << "  "
             << P1_F[1] << "  "
             << P1_F[2] << endl;
        cout << "v6f  = " << fixed << setprecision(2)
             << T_6toF[0][2] << "  "
             << T_6toF[1][2] << "  "
             << T_6toF[2][2] << endl;
        cout << "v67f = " << fixed << setprecision(2)
             << T_6toF[0][0] << "  "
             << T_6toF[1][0] << "  "
             << T_6toF[2][0] << endl;
        cout << "**********************" << endl;
        cout << endl;
    }
    return 0;
}




//int main() {
//    //Defining GE P60 Robot Constant Mechanism Parameters
//    //Twist Angles
//    double alpha_12,alpha_23,alpha_34,alpha_45,alpha_56;
//    alpha_12 = 270*D2R;
//    alpha_23 = 0;
//    alpha_34 = 0;
//    alpha_45 = 270*D2R;
//    alpha_56 = 90*D2R;
//    //Joint Offsets
//    double S2,S3,S4,S5;
//    S2 = 0; //cm
//    S3 = 0; //cm
//    S4 = 9.8; //cm
//    S5 = 14.5; //cm
//    //Link Length (cm)
//    double a12, a23, a34, a45, a56;
//    a12 = 0;
//    a23 = 70;
//    a34 = 90;
//    a45 = 0;
//    a56 = 0;
//
//    //Defining Variable Parameters (HW4)
//    double a67, alpha_67;
//    a67 = 0;
//    alpha_67 = 90*D2R;
//
//    //New For HW4
//    RobotKinematics Robot_GEP60(alpha_12,alpha_23,alpha_34,alpha_45,alpha_56,alpha_67,S2,S3,S4,S5,a12,a23,a34,a45,a56,a67);
//    //////////////////////////////////////////////////////////////////////////////////////////////////////
//    //Defining Variable Parameters (HW3)
//    double phi1, th2, th3, th4, th5, th6, S6;
//    phi1 = 50*D2R;
//    th2  = 120*D2R;
//    th3  = 295*D2R;
//    th4  = 30*D2R;
//    th5  = 190*D2R;
//    th6  = 100*D2R;
//    S6   = 15.24; // cm
//
//    //Create a T matrix and point P
//    double T_6toF[4][4] = {0};
//    double P1_F[4] = {0};
//    double P1_6[4] = {3.2, 4.1, 5.5, 1.0};
//
//    //Perform Forward Analysis
//    Robot_GEP60.Forward(T_6toF, phi1,th2,th3,th4,th5,th6,S6);
//
//    //Matrix Math Object
//    MathOperations MatrixOp;
//    //Getting the Point in Fixed Frame
//    MatrixOp.VectorMult(P1_F,T_6toF,P1_6);
//
//    //Output the Results
//    //For HW3
//    cout << "Tool Point in Fixed = " << P1_F[0] << ","
//         << P1_F[1] << "," << P1_F[2] << " cm" << endl ;
//    cout << "S6 in Fixed = " << T_6toF[0][2] << ","
//         << T_6toF[1][ 2] << "," << T_6toF[2][ 2] << endl ;
//    cout << "a67 in Fixed = " << T_6toF[0][0] << ","
//         << T_6toF[1][ 0] << "," << T_6toF[2][ 0] << endl ;
//    //////////////////////////////////////////////////////////////////////////////////////////////////////
////    //Inputs for Closed Loop Analysis (HW4)
////    //Values from Example 5.7 in book these are in inches
////    //a71,alpha_71,S7,th7,S1,gamma1 values to test case in PG 51
////    // Test General Case
////    double P_tool_6[3] = {5,3,7};
////    double P_tool_F[3] = {25,23,24};
////    double S6_F[3] = {0.177,0.884,-0.433};
////    double a67_F[3] = {-0.153,0.459,0.875};
////
////
////    //Closed Loop Analysis
////    double a71 = 0;
////    double S7 = 0;
////    double S1 = 0;
////    double alpha_71 = 0;
////    double th7 = 0;
////    double gamma1 = 0;
////    Robot_GEP60.Closed_Loop(P_tool_6, P_tool_F, S6_F,a67_F,
////                    &a71, &S7, &S1, &alpha_71, &th7, &gamma1);
////
////    //Testing Special Case 1: S1 and S7 are parallel
////    double P_tool_6_S1[3] = {5,3,2};
////    double P_tool_F_S1[3] = {30,11,24};
////    double S6_F_S1[3] = {1,0,0};
////    double a67_F_S1[3] = {0,1,0};
////    double a71_S1 = 0;
////    double S7_S1 = 0;
////    double S1_S1 = 0;
////    double alpha_71_S1 = 0;
////    double th7_S1 = 0;
////    double gamma1_S1 = 0;
////    Robot_GEP60.Closed_Loop(P_tool_6_S1, P_tool_F_S1, S6_F_S1,a67_F_S1,
////                            &a71_S1, &S7_S1, &S1_S1, &alpha_71_S1, &th7_S1, &gamma1_S1);
////
////    //Testing Special Case 2: S1 and S7 are collinear
////    double P_tool_6_S2[3] = {0,0,0};
////    double P_tool_F_S2[3] = {0,0,24};
////    double S6_F_S2[3] = {1,0,0};
////    double a67_F_S2[3] = {0,1,0};
////    double a71_S2 = 0;
////    double S7_S2 = 0;
////    double S1_S2 = 0;
////    double alpha_71_S2 = 0;
////    double th7_S2 = 0;
////    double gamma1_S2 = 0;
////    Robot_GEP60.Closed_Loop(P_tool_6_S2, P_tool_F_S2, S6_F_S2,a67_F_S2,
////                            &a71_S2, &S7_S2, &S1_S2, &alpha_71_S2, &th7_S2, &gamma1_S2);
////
////    //////////////////////////////////////////////////////////////////////////////////////////////////////
////    //Testing Values for practice exam 2022
////    double P_tool_5[3] = {0,0,0};
////    double P_tool_F[3] = {15,-12,10};
////    double S5_F[3] = {.177,.884,-.433};
////    double a56_F[3] = {-.153,.459,.875};
////
////    double a61 = 0;
////    double S6 = 0;
////    double S1 = 0;
////    double alpha_61 = 0;
////    double th6 = 0;
////    double gamma1 = 0;
////    Robot_GEP60.Closed_Loop(P_tool_5, P_tool_F, S5_F,a56_F,
////                            &a61, &S6, &S1, &alpha_61, &th6, &gamma1);
////    cout << "Practice Exam 2022" << endl;
////    cout << "Link Length a61 = " << a61 << endl;
////    cout << "Twist Angle Alpha 61 = " << alpha_61 <<  "\u00B0" << endl; // \u00B0 gives degree symbol
////    cout << "Joint Offset S6 = " << S6 << endl;
////    cout << "Joint Angle Theta 6 = " << th6 <<  "\u00B0" << endl;
////    cout << "Joint Offset S1  = " << S1<< endl;
////    cout << "Gamma 1 = " << gamma1 <<  "\u00B0" << endl;
////    //////////////////////////////////////////////////////////////////////////////////////////////////////
////
////    //Homework 5, Spherical Closed-Loop Mechanisms
////    //Spherical Mechanisms Object
////    SphericalMechanisms SphMech;
////
////    double A =  -0.1481;
////    double B = 0.7244;
////    double D =   0.6023;
////    double ang_a = 0;
////    double ang_b = 0;
////    SphMech.SolveTrig(A, B, D, &ang_a, &ang_b);
////
////   // cout << ang_a << " , " << ang_b << endl;
////
////    //Homework 7, Solve two bi-quadratic equations
////    //Using Values from HW7 Test Case
////
////    MathOperations MathOps;
////    double a1 = 4;
////    double b1 = 32;
////    double d1 = -11.2;
////    double e1 = -3.3;
////    double f1 = 4.6;
////    double g1 = 11.4;
////    double h1 = 5.5;
////    double i1 = -1.2;
////    double j1 = -3.3;
////
////    double a2 = 2;
////    double b2 = -11;
////    double d2 = -4.2;
////    double e2 = 1.3;
////    double f2 = 2.6;
////    double g2 = 7.4;
////    double h2 = -1.5;
////    double i2 = -2.2;
////    double j2 = -5.3;
////
////    //8th order Polynomial Coefficients
////    double C8,C7,C6,C5,C4,C3,C2,C1,C0;
////
////    //Input the constants a1-j2 and calculate coefficients C8-C0
////    MathOps.SolveBiQuadratic(a1,b1,d1,e1,f1,g1,h1,i1,j1,a2,b2,d2,e2,f2,g2,h2,i2,j2,
////                                &C8,&C7,&C6,&C5,&C4,&C3,&C2,&C1,&C0);
////    cout << "C8: " << C7<< endl;
////    //Take the Root of 8th order poly
////    double xcof[] = {C0,C1,C2,C3,C4,C5,C6,C7,C8};
////    int d = 8;  //degree of polynomial
////    double root_r[36] = {0};
////    double root_c[36] = {0};
////    MathOps.Poly_Solve(root_r, root_c, d, xcof);
////    for(int i = 0; i < 36; ++i) {
////        cout << root_c[i] << " ";
////    }
////    cout << endl;
////    //Find the associated outside roots
////    double x2_real[36] = {0};
////    double x2_imag[36] = {0};
////    SphMech.CalcOutsideRootFromInner(a1, b1, d1, e1, f1, g1, h1, i1, j1, a2, b2, d2, e2, f2, g2, h2, i2, j2, d, root_r, root_c, x2_real, x2_imag);
////
////
//////    //Testing to see if math functions work
//////    double test[3][3] = {{1,2,-1},{2,1,2},{-1,2,1}};
//////    double inv[3][3] = {0};
//////    cout << MathOps.MatrixDet_R3(test) << endl;
//////    MathOps.MatrixInv_R3(inv,test);
//////    for (int i = 0; i < 3; ++i) {          // Loop over rows
//////        for (int j = 0; j < 3; ++j) {      // Loop over columns
//////            cout << inv[i][j] << " ";
//////        }
//////        cout << endl;                       // Newline after each row
//////    }
////    cout << "Bi-Quadratic Equations" << endl;
////    cout << "Eq 1: (4 x1^2 + 32 x1 + -11.2) x2^2 + (-3.3 x1^2 + 4.6 x1 + 11.4) x2 + (5.5 x1^2 + -1.2 x1 + -3.3) = 0" << endl;
////    cout << "Eq 2: (2 x1^2 + -11 x1 + -4.2) x2^2 + (1.3 x1^2 + 2.6 x1 + 7.4) x2 + (-1.5 x1^2 + -2.2 x1 + -5.3) = 0" << endl;
////    cout << "ans 0:" << endl;
////    cout << "   x1 = " << root_r[0] << " , " <<  root_c[0] << " i" << endl;
////    cout << "   x2 = " << x2_real[0] << " , " <<  x2_imag[0] << " i" << endl;
////    cout << "ans 1:" << endl;
////    cout << "   x2 = " << root_r[1] << " , " <<  root_c[1] << " i" << endl;
////    cout << "   x2 = " << x2_real[1] << " , " <<  x2_imag[1] << " i" << endl;
////    cout <<"ans 2" << endl;
////    cout << "   x1 = " << root_r[2] << " , " <<  root_c[2] << " i" << endl;
////    cout << "   x2 = " << x2_real[2] << " , " <<  x2_imag[2] << " i" << endl;
////    cout << "ans 3" << endl;
////    cout << "   x1 = " << root_r[3] << " , " <<  root_c[3] << " i" << endl;
////    cout << "   x2 = " << x2_real[3] << " , " <<  x2_imag[3] << " i" << endl;
////    cout << "ans 4" << endl;
////    cout << "   x1 = " << root_r[4] << " , " <<  root_c[4] << " i" << endl;
////    cout << "   x2 = " << x2_real[4] << " , " <<  x2_imag[4] << " i" << endl;
////    cout << "ans 5"<< endl;
////    cout << "   x1 = " << root_r[5] << " , " <<  root_c[5] << " i" << endl;
////    cout << "   x2 = " << x2_real[5] << " , " <<  x2_imag[5] << " i" << endl;
////    cout << "ans 6" << endl;
////    cout << "   x1 = " << root_r[6] << " , " <<  root_c[6] << " i" << endl;
////    cout << "   x2 = " << x2_real[6] << " , " <<  x2_imag[6] << " i" << endl;
////    cout << "ans 7" << endl;
////    cout << "   x1 = " << root_r[7] << " , " <<  root_c[7] << " i" << endl;
////    cout << "   x2 = " << x2_real[7] << " , " <<  x2_imag[7] << " i" << endl;
////
////    //Output the Results
////    //For HW3
//////    cout << "Tool Point in Fixed = " << P1_F[0] << ","
//////         << P1_F[1] << "," << P1_F[2] << " cm" << endl ;
//////    cout << "S6 in Fixed = " << T_6toF[0][2] << ","
//////         << T_6toF[1][ 2] << "," << T_6toF[2][ 2] << endl ;
//////    cout << "a67 in Fixed = " << T_6toF[0][0] << ","
//////         << T_6toF[1][ 0] << "," << T_6toF[2][ 0] << endl ;
////
//////    //For HW4 verifying values from example 5.7 (General Case)
//////    cout << "Testing General Case" << endl;
//////    cout << "P_tool_6, P_tool_F,S6_F,a67_F values are from example 5.7 in book" << endl;
//////    cout << "Link Length a71 = " << a71 << endl;
//////    cout << "Twist Angle Alpha 71 = " << alpha_71 <<  "\u00B0" << endl; // \u00B0 gives degree symbol
//////    cout << "Joint Offset S7 = " << S7 << endl;
//////    cout << "Joint Angle Theta 7 = " << th7 <<  "\u00B0" << endl;
//////    cout << "Joint Offset S1  = " << S1<< endl;
//////    cout << "Gamma 1 = " << gamma1 <<  "\u00B0" << endl;
//////    cout << endl;
//////
////////    //Testing Special Case 1
//////    cout << "Testing Special Case 1: S1 and S7 are parallel" << endl;
//////    cout << "Input Values: P_tool_6 = {5,3,2}, P_tool_F = {30,11,24}, S6_F = {1,0,0}, a67_F = {0,1,0}" << endl;
//////    cout << "Link Length a71 = " << a71_S1 << endl;
//////    cout << "Twist Angle Alpha 71 = " << alpha_71_S1 <<  "\u00B0" << endl; // \u00B0 gives degree symbol
//////    cout << "Joint Offset S7 = " << S7_S1 << endl;
//////    cout << "Joint Angle Theta 7 = " << th7_S1 <<  "\u00B0" << endl;
//////    cout << "Joint Offset S1  = " << S1_S1 << endl;
//////    cout << "Gamma 1 = " << gamma1_S1 <<  "\u00B0" << endl;
//////    cout << endl;
//////    //Testing Special Case 2
//////    cout << "Testing Special Case 2: S1 and S7 are collinear" << endl;
//////    cout << "Input Values: P_tool_6 = {0,0,0}, P_tool_F = {0,0,24}, S6_F = {1,0,0}, a67_F = {0,1,0}" << endl;
//////    cout << "Link Length a71 = " << a71_S2 << endl;
//////    cout << "Twist Angle Alpha 71 = " << alpha_71_S2 <<  "\u00B0" << endl; // \u00B0 gives degree symbol
//////    cout << "Joint Offset S7 = " << S7_S2 << endl;
//////    cout << "Joint Angle Theta 7 = " << th7_S2 <<  "\u00B0" << endl;
//////    cout << "Joint Offset S1  = " << S1_S2 << endl;
//////    cout << "Gamma 1 = " << gamma1_S2 <<  "\u00B0" << endl;
//////    return 0;
////
//}
