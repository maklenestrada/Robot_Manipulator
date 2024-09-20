#include <iostream>
#include "RobotKinematics.h"
#include "MatrixMath.h"
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

//    //Defining Variable Parameters (HW4)
//    double a67, alpha_67;
//    a67 = 0;
//    alpha_67 = 90*D2R;

    RobotKinematics Robot_GEP60(alpha_12,alpha_23,alpha_34,alpha_45,alpha_56,S2,S3,S4,S5,a12,a23,a34,a45,a56);

//    //New For HW4
//    RobotKinematics Robot_GEP60(alpha_12,alpha_23,alpha_34,alpha_45,alpha_56,alpha_67,S2,S3,S4,S5,a12,a23,a34,a45,a56,a67);

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
    MatrixMath MatrixOp;
    //Getting the Point in Fixed Frame
    MatrixOp.VectorMult(P1_F,T_6toF,P1_6);

    //Inputs for Closed Loop Analysis (HW4)
    //Pretty Sure these are from the Foward Analysis
    //I just set up the vectors rn but might have to get the values for these
    //from the T matrix
    double P_tool_6[3] = {0};
    double P_tool_F[3] = {0};
    double S6_F[3] = {0};
    double a67_F[3] = {0};

//    //Closed Loop Analysis
//    Robot_GEP60.Closed_Loop(phi1,th2,th3,th4,th5,th6,S6,
//                              P_tool_6[3],P_tool_F[3],S6_F[3],a67_F[3]);

    //Output the Results
    cout << "Tool Point in Fixed = " << P1_F[0] << ","
         << P1_F[1] << "," << P1_F[2] << " cm" << endl ;
    cout << "S6 in Fixed = " << T_6toF[0][2] << ","
         << T_6toF[1][ 2] << "," << T_6toF[2][ 2] << endl ;
    cout << "a67 in Fixed = " << T_6toF[0][0] << ","
         << T_6toF[1][ 0] << "," << T_6toF[2][ 0] << endl ;

    return 0;

}
