#include <math.h>
#include <float.h>  // For DBL_MAX

/**
 * Performs Cholesky decomposition of a matrix C = FDF' where:
 *  - F is a lower triangular matrix with 1's on the diagonal
 *  - D is a diagonal matrix
 * 
 * Arguments:
 *     matrix - n x n matrix to be factored (input/output)
 *     n      - size of the matrix
 *     toler  - tolerance threshold to detect singularities
 * 
 * Return value:
 *     The rank of the matrix if it is positive semidefinite (SPD), or -rank if it is not SPD or NND
 *     If a column is redundant, its diagonal element is set to zero.
 *     NaN or infinite diagonal elements are treated as zero.
 *
 *    Terry Therneau
 */
int cholesky2(double **matrix, int n, double toler)
{
    double pivot, eps;
    int rank = 0;
    int nonneg = 1;  // Tracks whether the matrix is non-negative definite

    // Find the largest diagonal element to set a base tolerance level
    eps = 0.0;
    for (int i = 0; i < n; i++) {
        if (matrix[i][i] > eps) {
            eps = matrix[i][i];
        }
        // Ensure symmetry by copying the upper triangle to the lower triangle
        for (int j = i + 1; j < n; j++) {
            matrix[j][i] = matrix[i][j];
        }
    }

    // Set tolerance based on the largest diagonal element
    if (eps == 0.0) {
        eps = toler;  // Use provided tolerance if no positive diagonal element
    } else {
        eps *= toler;  // Adjust tolerance based on largest diagonal element
    }

    // Perform the Cholesky decomposition
    for (int i = 0; i < n; i++) {
        pivot = matrix[i][i];

        // If the pivot is too small or not finite, handle it
        if (!isfinite(pivot) || pivot < eps) {
            matrix[i][i] = 0.0;  // Set diagonal to 0 if pivot is small
            if (pivot < -8 * eps) {
                nonneg = -1;  // Mark as not positive definite if pivot is very negative
            }
        } else {
            rank++;  // Increase rank if pivot is valid
            for (int j = i + 1; j < n; j++) {
                double temp = matrix[j][i] / pivot;
                matrix[j][i] = temp;  // Store the L matrix (lower triangular part)

                // Update the diagonal and other elements
                matrix[j][j] -= temp * temp * pivot;
                for (int k = j + 1; k < n; k++) {
                    matrix[k][j] -= temp * matrix[k][i];
                }
            }
        }
    }

    // Return rank multiplied by the non-negative definiteness status
    return rank * nonneg;
}
