#include <iostream>
#include "RobotKinematics.h"
#include "MathOperations.h"
#include <math.h>

#define D2R M_PI/180.0

using namespace std;

int main() {
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

    //Defining Variable Parameters (HW4)
    double a67, alpha_67;
    a67 = 0;
    alpha_67 = 90*D2R;

    //New For HW4
    RobotKinematics Robot_GEP60(alpha_12,alpha_23,alpha_34,alpha_45,alpha_56,alpha_67,S2,S3,S4,S5,a12,a23,a34,a45,a56,a67);

    //Defining Variable Parameters (HW3)
    double phi1, th2, th3, th4, th5, th6, S6;
    phi1 = 50*D2R;
    th2  = 120*D2R;
    th3  = 295*D2R;
    th4  = 30*D2R;
    th5  = 190*D2R;
    th6  = 100*D2R;
    S6   = 15.24; // cm

    //Create a T matrix and point P
    double T_6toF[4][4] = {0};
    double P1_F[4] = {0};
    double P1_6[4] = {3.2, 4.1, 5.5, 1.0};

    //Perform Forward Analysis
    Robot_GEP60.Forward(T_6toF, phi1,th2,th3,th4,th5,th6,S6);

    //Matrix Math Object
    MathOperations MatrixOp;
    //Getting the Point in Fixed Frame
    MatrixOp.VectorMult(P1_F,T_6toF,P1_6);

    //Inputs for Closed Loop Analysis (HW4)
    //Values from Example 5.7 in book these are in inches
    //a71,alpha_71,S7,th7,S1,gamma1 values to test case in PG 51
    // Test General Case
    double P_tool_6[3] = {5,3,7};
    double P_tool_F[3] = {25,23,24};
    double S6_F[3] = {0.177,0.884,-0.433};
    double a67_F[3] = {-0.153,0.459,0.875};

    //Closed Loop Analysis
    double a71 = 0;
    double S7 = 0;
    double S1 = 0;
    double alpha_71 = 0;
    double th7 = 0;
    double gamma1 = 0;
    Robot_GEP60.Closed_Loop(P_tool_6, P_tool_F, S6_F,a67_F,
                    &a71, &S7, &S1, &alpha_71, &th7, &gamma1);

    //Testing Special Case 1: S1 and S7 are parallel
    double P_tool_6_S1[3] = {5,3,2};
    double P_tool_F_S1[3] = {30,11,24};
    double S6_F_S1[3] = {1,0,0};
    double a67_F_S1[3] = {0,1,0};
    double a71_S1 = 0;
    double S7_S1 = 0;
    double S1_S1 = 0;
    double alpha_71_S1 = 0;
    double th7_S1 = 0;
    double gamma1_S1 = 0;
    Robot_GEP60.Closed_Loop(P_tool_6_S1, P_tool_F_S1, S6_F_S1,a67_F_S1,
                            &a71_S1, &S7_S1, &S1_S1, &alpha_71_S1, &th7_S1, &gamma1_S1);

    //Testing Special Case 2: S1 and S7 are collinear
    double P_tool_6_S2[3] = {0,0,0};
    double P_tool_F_S2[3] = {0,0,24};
    double S6_F_S2[3] = {1,0,0};
    double a67_F_S2[3] = {0,1,0};
    double a71_S2 = 0;
    double S7_S2 = 0;
    double S1_S2 = 0;
    double alpha_71_S2 = 0;
    double th7_S2 = 0;
    double gamma1_S2 = 0;
    Robot_GEP60.Closed_Loop(P_tool_6_S2, P_tool_F_S2, S6_F_S2,a67_F_S2,
                            &a71_S2, &S7_S2, &S1_S2, &alpha_71_S2, &th7_S2, &gamma1_S2);


    //Output the Results
    //For HW3
//    cout << "Tool Point in Fixed = " << P1_F[0] << ","
//         << P1_F[1] << "," << P1_F[2] << " cm" << endl ;
//    cout << "S6 in Fixed = " << T_6toF[0][2] << ","
//         << T_6toF[1][ 2] << "," << T_6toF[2][ 2] << endl ;
//    cout << "a67 in Fixed = " << T_6toF[0][0] << ","
//         << T_6toF[1][ 0] << "," << T_6toF[2][ 0] << endl ;

    //For HW4 verifying values from example 5.7 (General Case)
    cout << "Testing General Case" << endl;
    cout << "P_tool_6, P_tool_F,S6_F,a67_F values are from example 5.7 in book" << endl;
    cout << "Link Length a71 = " << a71 << endl;
    cout << "Twist Angle Alpha 71 = " << alpha_71 <<  "\u00B0" << endl; // \u00B0 gives degree symbol
    cout << "Joint Offset S7 = " << S7 << endl;
    cout << "Joint Angle Theta 7 = " << th7 <<  "\u00B0" << endl;
    cout << "Joint Offset S1  = " << S1<< endl;
    cout << "Gamma 1 = " << gamma1 <<  "\u00B0" << endl;
    cout << endl;

//    //Testing Special Case 1
    cout << "Testing Special Case 1: S1 and S7 are parallel" << endl;
    cout << "Input Values: P_tool_6 = {5,3,2}, P_tool_F = {30,11,24}, S6_F = {1,0,0}, a67_F = {0,1,0}" << endl;
    cout << "Link Length a71 = " << a71_S1 << endl;
    cout << "Twist Angle Alpha 71 = " << alpha_71_S1 <<  "\u00B0" << endl; // \u00B0 gives degree symbol
    cout << "Joint Offset S7 = " << S7_S1 << endl;
    cout << "Joint Angle Theta 7 = " << th7_S1 <<  "\u00B0" << endl;
    cout << "Joint Offset S1  = " << S1_S1 << endl;
    cout << "Gamma 1 = " << gamma1_S1 <<  "\u00B0" << endl;
    cout << endl;
    //Testing Special Case 2
    cout << "Testing Special Case 2: S1 and S7 are collinear" << endl;
    cout << "Input Values: P_tool_6 = {0,0,0}, P_tool_F = {0,0,24}, S6_F = {1,0,0}, a67_F = {0,1,0}" << endl;
    cout << "Link Length a71 = " << a71_S2 << endl;
    cout << "Twist Angle Alpha 71 = " << alpha_71_S2 <<  "\u00B0" << endl; // \u00B0 gives degree symbol
    cout << "Joint Offset S7 = " << S7_S2 << endl;
    cout << "Joint Angle Theta 7 = " << th7_S2 <<  "\u00B0" << endl;
    cout << "Joint Offset S1  = " << S1_S2 << endl;
    cout << "Gamma 1 = " << gamma1_S2 <<  "\u00B0" << endl;
    return 0;

}
