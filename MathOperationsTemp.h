#ifndef ROBOT_MANIPULATOR_MATHOPERATIONSTEMP_H
#define ROBOT_MANIPULATOR_MATHOPERATIONSTEMP_H

#include <complex>
#include <stdexcept>

// Template function for calculating the determinant of a 3x3 matrix
template <typename T>
T MatrixDet(const T Mat[3][3]) {
    return Mat[0][0] * (Mat[1][1] * Mat[2][2] - Mat[1][2] * Mat[2][1]) -
           Mat[0][1] * (Mat[1][0] * Mat[2][2] - Mat[1][2] * Mat[2][0]) +
           Mat[0][2] * (Mat[1][0] * Mat[2][1] - Mat[1][1] * Mat[2][0]);
}

// Template function for inverting a 3x3 matrix
template <typename T>
void MatrixInv(T Inv[3][3], const T Mat[3][3]) {
    T Det = MatrixDet(Mat);
    if (Det == T(0)) {
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
    }

    // Calculate the inverse using the cofactor matrix divided by the determinant
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

// Template function for multiplying a 3x3 matrix by a 3x1 vector
template <typename T>
void MatrixVecMult(T result[3], const T Mat[3][3], const T Vec[3]) {
    for (int i = 0; i < 3; ++i) {
        result[i] = T(0);  // Initialize the result to zero
        for (int j = 0; j < 3; ++j) {
            result[i] += Mat[i][j] * Vec[j];
        }
    }
}

#endif //ROBOT_MANIPULATOR_MATHOPERATIONSTEMP_H
