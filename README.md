# Robot_Manipulator

This project implements robot kinematics in C++. The code is structured into multiple source and header files to demonstrate modularity and proper code organization.

## Project Structure

├── MatrixMath.cpp        # Contains matrix operations like matrix-vector multiplication  
├── MatrixMath.h          # Header file for matrix operations  
├── RobotKinematics.cpp   # Implements the RobotKinematics class  
├── RobotKinematics.h     # Header file for the RobotKinematics class  
├── main.cpp              # Main program that utilizes both the RobotKinematics class and matrix functions  
├── README.md             # This README file  

## Files Overview

- **`MatrixMath.h`/`MatrixMath.cpp`**: Defines and implements matrix math functions such as `VectorMult()` for matrix-vector multiplication.
  
- **`RobotKinematics.h`/`RobotKinematics.cpp`**: Implements the `RobotKinematics` class, which handles robot kinematics calculations (e.g., forward kinematics).

- **`main.cpp`**: The main program file where the `RobotKinematics` class and matrix math functions are used together.

## Current Status (As of September 20, 2024)

- **Forward Kinematics Analysis**: The project currently implements forward kinematics analysis, allowing the calculation of the end effector's position based on given joint angles.
  
- **Work in Progress**: Currently working on implementing the **reverse kinematics analysis**, which will determine the required joint angles to achieve a desired end-effector position.


