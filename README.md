[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-dark-v6.png#gh-dark-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/bilgesubilgin/Signaloid-Demo-Basic-GaussianElimination#gh-dark-mode-only)
[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-light-v6.png#gh-light-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/bilgesubilgin/Signaloid-Demo-Basic-GaussianElimination#gh-light-mode-only)


# Solve System of Linear Equations via Gaussian Elimination
## Overview
The code in this repository allows one to solve a specified system of linear equations Ax = y, where the entries of the coefficient matrix A and the solution vector y are distributions, via the Gaussian Elimination method implemented to run on the Signaloid Cloud Server. It outputs the solution only in the case when the system has a unique solution and it just prints out a description for other cases. To gauge the 'correctness' of the found unique solution, the code also plugs this solution into the provided system and finds a predicted value vector y' to be compared to the given value vector y.

## Inputs
The generic inputs are the (m x n) coefficient matrix A and the (m x 1) value vector y.
###For the provide example: 
- All inputs are assumed to be uniformly distributed with range size equal to 1.
- Expected values of the coefficient matrix A and the value vector y are:

         ⎛ 2.0   1.0   -1.0 ⎞                  ⎛   8.0 ⎞
  E[A] = ⎜-3.0  -1.0    2.0 ⎟   ,       E[y] = ⎜ -11.0 ⎟   .
         ⎝-2.0   1.0    2.0 ⎠                  ⎝  -3.0 ⎠

## Outputs
In the case of a unique solution the code outputs the (n x 1) solution vector x that satisfies Ax = y. Else, the method prints out the case for the solution, either no solutions or infinitely many.
It also outputs the value vector y' that is found after plugging the solution into the equations.

 When there is no uncertainty, the solution x to the above example system of equations is:
 
         ⎛  2.0 ⎞
  E[x] = ⎜  3.0 ⎟   .
         ⎝ -1.0 ⎠

## Repository Tree Structure
```
.
├── README.md
└── src
    └── solveViaGaussianElimination.c
```
