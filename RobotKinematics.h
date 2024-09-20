#ifndef ROBOT_MANIPULATOR__ROBOTKINEMATICS_H
#define ROBOT_MANIPULATOR__ROBOTKINEMATICS_H

#include <iostream>

class RobotKinematics{
private:
    //Defining GE P60 Robot Constant Mechanism Parameters
    //Twist Angles
    double alpha_12,alpha_23,alpha_34,alpha_45,alpha_56;
    //Joint Offsets
    double S2,S3,S4,S5;
    //Link Length (cm)
    double a12, a23, a34, a45, a56;

public:
    //Constructor
    RobotKinematics(double alpha_12, double alpha_23, double alpha_34, double alpha_45, double alpha_56,
                    double S2, double S3, double S4, double S5,
                    double a12, double a23, double a34, double a45, double a56);

    //Function for the Forward Analysis
    void Forward(double T_6toF[4][4], double phi1,
                double th2, double th3, double th4,
                double th5, double th6, double S6);
    //Function to create the general T_itoj Matrix
    void T_itoj(double th_j,double alpha_ij,double a_ij, double S_j,double T_ij[4][4]);

    //Getter Methods to access private parameters
    double Getalpha_12() const;
    double Getalpha_23() const;
    double Getalpha_34() const;
    double Getalpha_45() const;
    double Getalpha_56() const;

    double GetS2() const;
    double GetS3() const;
    double GetS4() const;
    double GetS5() const;

    double Geta12() const;
    double Geta23() const;
    double Geta34() const;
    double Geta45() const;
    double Geta56() const;
};
#endif //ROBOT_MANIPULATOR__ROBOTKINEMATICS_H
