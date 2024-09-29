#include "RobotKinematics.h"
#include "MathOperations.h"
#include <math.h>
#include <vector>
#include <cmath>

#define D2R M_PI/180.0
#define R2D 1/D2R
using namespace std;

////Constructor Implementation
//RobotKinematics::RobotKinematics(double alpha_12, double alpha_23, double alpha_34, double alpha_45, double alpha_56,
//                                 double S2, double S3, double S4, double S5, double a12, double a23, double a34,
//                                 double a45, double a56) : alpha_12(alpha_12), alpha_23(alpha_23), alpha_34(alpha_34),
//                                                           alpha_45(alpha_45), alpha_56(alpha_56), S2(S2), S3(S3), S4(S4), S5(S5), a12(a12),
//                                                           a23(a23), a34(a34), a45(a45), a56(a56) {}

//New Constructor Implementation (HW4)
RobotKinematics::RobotKinematics(double alpha_12, double alpha_23, double alpha_34, double alpha_45, double alpha_56,double alpha_67,
                                 double S2, double S3, double S4, double S5, double a12, double a23, double a34,
                                 double a45, double a56, double a67) : alpha_12(alpha_12), alpha_23(alpha_23), alpha_34(alpha_34),
                                 alpha_45(alpha_45), alpha_56(alpha_56),alpha_67(alpha_67), S2(S2), S3(S3), S4(S4), S5(S5), a12(a12),
                                 a23(a23), a34(a34), a45(a45), a56(a56), a67(a67) {}

//Forward Analysis
void RobotKinematics::Forward(double T_6toF[4][4], double phi1,
                              double th2, double th3, double th4,
                              double th5, double th6, double S6)
{
    //Defining T_Fto1
    double T_1toF[4][4];

    //Precompute sin and cos values
    double c1 = cos(phi1);
    double s1 = sin(phi1);

    //Populate the 4x4 Matrix
    // Row 1
    T_1toF[0][0] = c1;
    T_1toF[0][1] = -s1;
    T_1toF[0][2] = 0;
    T_1toF[0][3] = 0;
    // Row 2
    T_1toF[1][0] = s1;
    T_1toF[1][1] = c1;
    T_1toF[1][2] = 0;
    T_1toF[1][3] = 0;
    // Row 3
    T_1toF[2][0] = 0;
    T_1toF[2][1] = 0;
    T_1toF[2][2] = 1;
    T_1toF[2][3] = 0;
    // Row 4
    T_1toF[3][0] = 0;
    T_1toF[3][1] = 0;
    T_1toF[3][2] = 0;
    T_1toF[3][3] = 1;

    //Defining T_2to1, T_3to2, T_4to3,T_5to4,T_6to5
    //T_2to1
    double T_2to1[4][4] = {0};
    T_itoj(th2,alpha_12,a12,S2,T_2to1);
    //T_3to2
    double T_3to2[4][4] = {0};
    T_itoj(th3,alpha_23,a23,S3,T_3to2);
    //T_4to3
    double T_4to3[4][4] = {0};
    T_itoj(th4,alpha_34,a34,S4,T_4to3);
    //T_5to4
    double T_5to4[4][4] = {0};
    T_itoj(th5,alpha_45,a45,S5,T_5to4);
    //T_6to5
    double T_6to5[4][4] = {0};
    T_itoj(th6,alpha_56,a56,S6,T_6to5);

    //MathOperations Object
    MathOperations MatrixOp;

    //Multiply Matricies
    //T_5to4 * T_6to5
    double T_6to4[4][4] = {0};
    MatrixOp.MatrixMult(T_6to4,T_5to4,T_6to5);
    //T_4to3 * T_6to4
    double T_6to3[4][4] = {0};
    MatrixOp.MatrixMult(T_6to3,T_4to3,T_6to4);
    //T_3to2 * T_6to3
    double T_6to2[4][4] = {0};
    MatrixOp.MatrixMult(T_6to2,T_3to2,T_6to3);
    //T_2to1 * T_6to2
    double T_6to1[4][4] = {0};
    MatrixOp.MatrixMult(T_6to1,T_2to1,T_6to2);
    //T_1toF * T_6to1
    MatrixOp.MatrixMult(T_6toF,T_1toF,T_6to1);
}

//T_itoj Matrix
void RobotKinematics::T_itoj(double th_j,double alpha_ij,double a_ij, double S_j,double T_ij[4][4])
{
    //Precompute sin and cos
    double cj,sj, cij,sij;
    cj = cos(th_j);
    sj = sin(th_j);
    cij = cos(alpha_ij);
    sij = sin(alpha_ij);
    double Sj,aij;
    Sj = S_j;
    aij = a_ij;

    //Populate 4x4 T_ij Matrix
    //Row 1
    T_ij[0][0] = cj;
    T_ij[0][1] = -sj;
    T_ij[0][2] = 0;
    T_ij[0][3] = aij;
    //Row 2
    T_ij[1][0] = sj * cij;
    T_ij[1][1] = cj * cij;
    T_ij[1][2] = -sij;
    T_ij[1][3] = -sij * Sj;
    //Row 3
    T_ij[2][0] = sj * sij;
    T_ij[2][1] = cj * sij;
    T_ij[2][2] = cij;
    T_ij[2][3] = cij * Sj;
    //Row 4
    T_ij[3][0] = 0;
    T_ij[3][1] = 0;
    T_ij[3][2] = 0;
    T_ij[3][3] = 1;
}

//Closed Loop Analysis
void RobotKinematics::Closed_Loop(double phi1, double th2, double th3,
                 double th4, double th5, double th6,
                 double S6, double P_tool_6[3],
                 double P_tool_F[3], double S6_F[3], double a67_F[3])
{
    //MathOperations Object
    MathOperations MathOps;

    //Defining vector S7 (in fixed frame)
    //By definition must be perpendicular to a67 and S6
    double S7_F[3] = {0};
    MathOps.CrossProduct(S7_F,a67_F,S6_F);

    //Defining vector S1 (in fixed frame and is parallel to Z axis)
    double S1_F[3] = {0,0,1};

    //Calculating cos and sin values with twist angle alpha_71
    double c71;
    MathOps.DotProduct(&c71,S7_F,S1_F);

    //General Case
    if(!MathOps.ValueNear(abs(c71), 1.0, 0.0001))
    {
        //Calculate a71_F (EQ 5.10)
        double a71_F[3];
        double step1_a71[3],step2_a71,step3_a71;
        MathOps.CrossProduct(step1_a71,S7_F,S1_F);
        MathOps.Norm(&step2_a71,step1_a71);
        step3_a71 = 1/step2_a71;
        MathOps.MultVecByScalar(a71_F,step1_a71,step3_a71);

        //Calculate s71 (EQ 5.12)
        double s71;
        double step1_s71[3];
        MathOps.CrossProduct(step1_s71,S7_F,S1_F);
        MathOps.DotProduct(&s71,step1_s71,a71_F);

        //Calculating cos and sin values of joint angle theta 7 (th7) (EQ 5.13 & 5.14)
        double c7;
        MathOps.DotProduct(&c7,a67_F,a71_F);
        double s7;
        double step1_s7[3];
        MathOps.CrossProduct(step1_s7,a67_F,a71_F);
        MathOps.DotProduct(&s7,step1_s7,S7_F);

        //Calculating the cos and sin of gamma1 (EQ 5.15 & 5.16)
        double cg1;
        double X_F[3] = {1,0,0};
        MathOps.DotProduct(&cg1,a71_F,X_F);
        double sg1;
        double step1_sg1[3];
        MathOps.CrossProduct(step1_sg1,a71_F,X_F);
        MathOps.DotProduct(&sg1,step1_sg1,S1_F);

        //Calculating Joint Offset Value S7 (EQ 5.21)
        double S7;
        double step1_S7[3], step2_S7;
        double P6o_F[3] = {1,2,3}; // Okay Ithink this is the origin of the sixth frame but not sure go to office hours and ask what this is or ask mauro tomorrow
        MathOps.CrossProduct(step1_S7,S1_F,P6o_F);
        MathOps.DotProduct(&step2_S7,step1_S7,a71_F);
        S7 = step2_S7/s71;

        //Calculate Link Length a71 (EQ 5.22)
        double a71;
        double step1_a71l[3], step2_a71l; //l for length
        MathOps.CrossProduct(step1_a71l,P6o_F,S1_F);
        MathOps.DotProduct(&step2_a71l,step1_a71l,S7_F);
        a71 = step2_a71l/s71;

        //Calculate Joint Offset Value S1
        double S1;
        double step1_S1[3],step2_S1;
        MathOps.CrossProduct(step1_S1,P6o_F,S7_F);
        MathOps.DotProduct(&step2_S1,step1_S1,a71_F);
        S1 = step2_S1/s71;
    }
    else
    {
        //special case 1
        //
        //
        //
        //
//        if(a71 != 0)
//        {
//            //Regular special case1
//        }
//        else
//        {
//            //Special case 2
//        }

    }


}

//Getter Methods Implementation
double RobotKinematics::Getalpha_12() const { return alpha_12; }
double RobotKinematics::Getalpha_23() const { return alpha_23; }
double RobotKinematics::Getalpha_34() const { return alpha_34; }
double RobotKinematics::Getalpha_45() const { return alpha_45; }
double RobotKinematics::Getalpha_56() const { return alpha_56; }
double RobotKinematics::Getalpha_67() const { return alpha_67; }

double RobotKinematics::GetS2() const { return S2; }
double RobotKinematics::GetS3() const { return S3; }
double RobotKinematics::GetS4() const { return S4; }
double RobotKinematics::GetS5() const { return S5; }

double RobotKinematics::Geta12() const { return a12; }
double RobotKinematics::Geta23() const { return a23; }
double RobotKinematics::Geta34() const { return a34; }
double RobotKinematics::Geta45() const { return a45; }
double RobotKinematics::Geta56() const { return a56; }
double RobotKinematics::Geta67() const { return a67; }