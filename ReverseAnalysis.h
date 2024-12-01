#ifndef ROBOT_MANIPULATOR_REVERSEANALYSIS_H
#define ROBOT_MANIPULATOR_REVERSEANALYSIS_H


class ReverseAnalysis {
public:
    //Constructor
    ReverseAnalysis() {}

    void ReverseAnalysis_GEP60(double P_tool_6[3],double P_tool_F[3],double S6_F[3],double a67_F[3],double S6,
            double* valid,double th[6][8]);



};


#endif //ROBOT_MANIPULATOR_REVERSEANALYSIS_H
