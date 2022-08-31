#include <math.h>
#include <stdio.h>
#include <uncertain.h>

/*
 *    Solve a linear system of equations using the Gaussian elimination method.
 *
 *    Overview:
 *    Consider a linear system of equations Ax = y, where A is the (m x n) coefficient matrix,
 *    x is a (n x 1) vector, and y is the (m x 1) value vector. Given the coefficient matrix A
 *    and the value vector y, this method outputs the solution vector x provided that the system
 *    has a unique solution. In case of no or inifinitely many solutions, it just prints out the case.
 *
 *    Inputs:
 *    - All inputs are assumed to be uniformly distributed with range size equal to 1.
 *    - Expected values of the coefficient matrix A and the value vector y are:
 *
 *            ⎛ 2.0   1.0   -1.0 ⎞                  ⎛   8.0 ⎞
 *     E[A] = ⎜-3.0  -1.0    2.0 ⎟   ,       E[y] = ⎜ -11.0 ⎟   .
 *            ⎝-2.0   1.0    2.0 ⎠                  ⎝  -3.0 ⎠
 *
 *    Outputs:
 *    In the case of a unique solution, the output is the solution vector x that satisfies Ax = y.
 *    Else, the method prints out the case for the solution, either no solutions or infinitely many.
 *    It also outputs the value vector y' that is found after plugging the solution into the equations.
 *
 *    When there is no uncertainty, the solution x to the above example system of equations is:
 *
 *            ⎛  2.0 ⎞
 *     E[x] = ⎜  3.0 ⎟   .
 *            ⎝ -1.0 ⎠
 */



// Intialize the distirbutions for the entries of the augmented matrix describing the linear system of equations.
static void
loadDistributions (int noEqs, int noVars, double matExpVals[][noVars], double * vecExpVals, double mat[][noVars], double * vec)
{
    // Fixed halfwidth for the uniform distributions of the system parameters.
    double halfWidth = 0.5, lowValue, highValue;
    for (int i = 0; i < noEqs; i++)
    {
        for (int j = 0; j < noVars; j++)
        {
            lowValue = matExpVals[i][j] - halfWidth;
            highValue = matExpVals[i][j] + halfWidth;
            mat[i][j] = libUncertainDoubleUniformDist(lowValue, highValue);
        }
        lowValue = vecExpVals[i] - halfWidth;
        highValue = vecExpVals[i] + halfWidth;
        vec[i] = libUncertainDoubleUniformDist(lowValue, highValue);
    }
}



// Find and return the argmax for the absolute values of entries of a given vector.
static int
argmaxAbs (int lenVector, double * vector)
{
    int i_max = 0;
    double maxValue = fabs(vector[0]);
    for (int i = 1; i < lenVector; i++)
    {
        if (fabs(vector[i]) > maxValue)
        {
            maxValue = fabs(vector[i]);
            i_max = i;
        }
    }
    return i_max;
}


// Matrix multiplication to test the found solution.
static void
matrixMultiply (int noRows, int noCols, double matrix[][noCols], double * inputVector, double * outputVector)
{
    double sum;
    for (int i = 0; i < noRows; i++)
    {
        sum = 0;
        for (int j = 0; j < noCols; j++)
        {
            sum = sum + matrix[i][j] * inputVector[j];
        }
        outputVector[i] = sum;
    }
}


// Reduce the given augmented matrix to the row echelon form.
static void
reduceToRowEchelon (int noRows, int noCols, double augmentedMatrix[][noCols])
{
    // Initialize the pivot row and column indices.
    int ii = 0, jj = 0;
    double tempVal;
    
    while (ii < noRows && jj < noCols - 1)
    {
        // Find the vector consisting of the entries in the jj-th column starting from the ii-th row.
        double jjthColumn[noRows - ii];
        for (int i = 0; i < noRows - ii; i++)
        {
            jjthColumn[i] = augmentedMatrix[ii + i][jj];
        }
        
        // Set the pivot in the jj-th column as the non-zero entry with the maximum absolute value (for numerical stability).
        int i_max = argmaxAbs(noRows - ii, jjthColumn) + ii;
        if (augmentedMatrix[i_max][jj] == 0)
        {
            // No pivots in this column.
            jj++;
        }
        else
        {
            // Swap the current row with the row with largest absolute pivot value.
            if (i_max != ii)
            {
                for (int j = jj; j < noCols; j++)
                {
                    tempVal = augmentedMatrix[i_max][j];
                    augmentedMatrix[i_max][j] = augmentedMatrix[ii][j];
                    augmentedMatrix[ii][j] = tempVal;
                }
            }
            // Do for all the rows below the pivot row:
            for (int i = ii + 1; i < noRows; i++)
            {
                // No need to handle rows with 0 entries in the pivot column.
                if (augmentedMatrix[i][jj] == 0)
                {
                    continue;
                }
                // Find the scale factor for the row.
                double factor = augmentedMatrix[i][jj] / augmentedMatrix[ii][jj];
                // Fill the lower part of the pivot with zeros.
                augmentedMatrix[i][jj] = 0;
                // Do for all remaining entries of the current row:
                for (int j = jj + 1; j < noCols; j++)
                {
                    augmentedMatrix[i][j] = augmentedMatrix[i][j] - factor * augmentedMatrix[ii][j];
                }
            }
            // Increase the pivot row and column.
            ii++;
            jj++;
        }
    }
}



// Find the solution to the linear system of equations given in row echelon form via backward substitution in the case of a
// unique solution.
static void
backwardSubstitute (int noVars, double rowEchelonMatrix[][noVars + 1], double solution[noVars])
{
    for (int i = noVars - 1; i >= 0; i--)
    {
        double value = rowEchelonMatrix[i][noVars];
        for (int j = noVars - 1; j > i; j--)
        {
            value = value - rowEchelonMatrix[i][j] * solution[j];
        }
        solution[i] = value / rowEchelonMatrix[i][i];
    }
}



// Solve the given linear system of equations via the Gaussian elimination method and return an integer encoding the case
// for the solution, i.e., -2 for trivial equation, -1 for no solutions, 0 for a unique solution, and 1 for inifinetely many
// solutions.
static int
solveWithGaussianElimination (int noEqs, int noVars, double coefficientMatrix[][noVars], double * valueVector, double * solution)
{
    // Construct the augmented matrix by concatenating the coefficient matrix and the value vector.
    // Number of rows and columns in the augmented matrix.
    int noRows = noEqs, noCols = noVars + 1;
    // The augmented matrix
    double augmentedMatrix[noRows][noCols];
    for (int i = 0; i < noRows; i++)
    {
        for (int j = 0; j < noCols; j++)
        {
            if (j == noCols - 1)
            {
                augmentedMatrix[i][j] = valueVector[i];
            }
            else
            {
                augmentedMatrix[i][j] = coefficientMatrix[i][j];
            }
        }
    }
    
    // Reduce the augmented matrix to row Echelon form.
    reduceToRowEchelon(noRows, noCols, augmentedMatrix);
    
    // Identify the number of independent equations by identifying the all-zero rows of the coefficients part of the reduced
    // augmented matrix by sweeping through the rows from the end. After reducing to the row echelon form, all-zero rows of
    // the coeffcient part of the matrix (if any) are the most bottom ones.
    int noIndependentEquations;
    int stopSweep = 0;
    for (int i = noRows - 1; i >= 0; i--)
    {
        for (int j = 0; j < noCols - 1; j++)
        {
            if (augmentedMatrix[i][j] != 0)
            {
                stopSweep = 1;
                //
                noIndependentEquations = i + 1;
                break;
            }
        }
        if (stopSweep == 1)
        {
            // We have found a non-all-zero row of the coefficient part of the reduced augmented matrix.
            // The rows above are guaranteed to be non-all-zero as well, so we exit the sweep.
            break;
        }
        else
        {
            // The coefficients part of the current row is all zero.
            if (augmentedMatrix[i][noCols - 1] != 0)
            {
                // If the value part of the row is not zero, then there are no solutions to the system.
                return -1;
            }
            else
            {
                // If i == 0, then all entries of the augmented matrix are zero, i.e., the system is trivial.
                if (i == 0)
                {
                    return -2;
                }
            }
        }
    }
    
    // If the number of independent equations is less than the number of variables, then there are infinitely many solutions.
    // I will provide a characterization for these solutions in the second pass!!!
    if (noIndependentEquations < noVars)
    {
        return 1;
    }
    // Else, the number of independent equations must equal the number of variables, which is the case of a uniqe solution.
    else
    {
        // Find the solution from the reduced augmented matrix via backward substitution.
        backwardSubstitute(noVars, augmentedMatrix, solution);
        return 0;
    }
}



// Find the solution to the given linear system of equations described by the augmented matrix via Gaussian elimination.
int
main ()
{
    // Initialize the linear system of equations Ax = y:
    // The expected values of the entries in the coefficient function A and the value vector y.
    double coefficientMatrixExpValues[3][3] = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
    double valueVectorExpValues[3] = {8, -11, -3};
    
    // Number of equations and variables in the system.
    int noEqs = sizeof(coefficientMatrixExpValues) / sizeof(coefficientMatrixExpValues[0]);
    int noVars = sizeof(coefficientMatrixExpValues[0]) / sizeof(double);
    
    // Size of the value vector should be equal to the number of equations.
    int sizeValueVector = sizeof(valueVectorExpValues) / sizeof(double);
    if (sizeValueVector != noEqs)
    {
        printf("\nInvalid system of equations!\n");
        printf("The number of rows of the coefficient matrix must equal the size of the value vector.\n\n");
        return 0;
    }
    
    // Load the distributions.
    double coefficientMatrix[noEqs][noVars];
    double valueVector[noEqs];
    loadDistributions (noEqs, noVars, coefficientMatrixExpValues, valueVectorExpValues, coefficientMatrix, valueVector);
    
    // Print the inputs:
    // Print the coefficient matrix A.
    printf("\nInputs to the system of linear equations Ax = y:\n");
    printf("\nCoefficient Matrix (A)\n");
    printf("----------------------\n");
    for (int i = 0; i < noEqs; i++)
    {
        for (int j = 0; j < noVars; j++)
        {
            printf("%.2f\t", coefficientMatrix[i][j]);
        }
        printf("\n");
    }
    // Print the value vector y.
    printf("\nValue Vector (y)\n");
    printf("----------------\n");
    for (int i = 0; i < noEqs; i++)
    {
        printf("y_%d = %.2f\n", i, valueVector[i]);
    }
    
    // Solve the system via Gaussian elimination and get the case for the solution.
    double solution[noVars];
    int solutionCase;
    solutionCase = solveWithGaussianElimination(noEqs, noVars, coefficientMatrix, valueVector, solution);
    
    // Print the case for the solution, and in the case of a unique solution, print the solution as well as
    // the predicted value vector obtained when the solution is plugged into the equation Ax = y.
    // The case of the trivial equation, that is, A = y = 0:
    if (solutionCase == -2)
    {
        printf("\nThe given system is trivial!\n");
    }
    // The case of no solutions:
    else if (solutionCase == -1)
    {
        printf("\nThere are no solutions!\n");
    }
    // The case of a unique solution:
    else if (solutionCase == 0)
    {
        printf("\nThere is a unique solution:\n");
        
        // Print the solution.
        printf("\nSolution\n");
        printf("--------\n");
        for (int i = 0; i < noVars; i++)
        {
            printf("x_%d = %.2f\n", i, solution[i]);
        }
        
        // Find the predicted value vector after plugging in the solution to the equation Ax = y.
        double predictedValueVector[noEqs];
        matrixMultiply(noEqs, noVars, coefficientMatrix, solution, predictedValueVector);
        
        // Print the predicted and the actual value vectors.
        printf("\nValue vector obtained by plugging the solution into the equation Ax = y:\n");
        printf("\nPredicted Value Vector (y')\n");
        printf("---------------------------\n");
        for (int i = 0; i < noEqs; i++)
        {
            printf("y'_%d = %.2f\n", i, predictedValueVector[i]);
        }
    }
    // The case of infinitely many solutions:
    else
    {
        printf("\nThere are infinitely many solutions!\n");
    }
    printf("\n");
    
    return 0;
}
