# A regularized cutting plane approach for inverse MIP
This repository contains the code and instances associated with the following paper

"A regularized cutting plane algorithm for inverse mixed integer programming Problem" 

by Vishnu Vijayaraghavan, Kiavash Kianfar, Andrew Schaefer, submitted, 2019

## Implemented Algorithms 
    
      1. Cutting Plane Algorithm presented in "Wang, L., Cutting plane algorithms for the inverse mixed integer linear programming problem, Operations Research Letters 37: 114-116, 2009  
      
      2. Regularized Cutting Plane Algorithm developed in "Vijayaraghavan, V., Kianfar, K., Schaefer, A., A regularized cutting plane algorithm for inverse mixed integer programming Problem"  submitted, 2019

## Code

    Implementation in C++ with Optimization Problems solved using ILOG CPLEX 12.9 callable library
    on Windows Platform

    Files : 
       
       1. ProbMIP.cpp : Sets up the MIP/ Regularized MIP problems
                      
                        Implements functions to randomly generate Inverse MIP instances based on the Seed obtained
                        
                        Includes functions to solve the MIP instances and obtain the solutions

       2. GenLP.cpp   : Sets up the Cost Generator LP problem 
                         
                        Includes functions to add new cuts to the LP problem instance

                        Also includes functions to solve the LP instances and obtain the solutions  

       3. Main.cpp    : Includes the main() function that implements the algorithms

                        Takes the problem name as a command line argument while execution

                        The variable seed needs to be set to a positive integer value for the algorithm to generate random instances
                        of the Inverse MIP problem
 
## Instructions to Execute the Code

1. Include the files Main.cpp, ProbMIP.cpp, ProbMIP.h, GenLP.cpp, and GenLp.h

2. Pass the MIP problem to be used as a command line argument when calling the executable

3. The variable seed in Main.cpp to be set to a positive integer to generate a random instance,
   a random instances being a combination of randomly generated objective coefficient and required solution for the Inverse MIP problem

4. To execute the Cutting Plane Algorithm include the directive #define SOLVE_MIPCUT before the main() function in Main.cpp

5. To execute the Regularized Cutting Plane Algorithm include the directive #define SOLVE_MIPREG before the main() function in Main.cpp









