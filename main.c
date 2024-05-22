/*** Abdurrahim Yılmaz ***/
/*** ***/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/*** MATRIX STRUCTURE ***/
/* A matrix structure that is suitable for dynamic matrix operations */
typedef struct {
    float **M;  // Pointer for a 2-dimensional (2D) float type array
    int n, m;   // Numbers of rows and columns for the 2D array
} Matrix;



/*** FUNCTION PROTOTYPES ***/
Matrix construct(Matrix x, int n, int m);   // Returns a Matrix variable by allocating memory for n-by-m sized 2D array
Matrix destruct(Matrix x);                  // Returns a Matrix variable by deallocating the memory allocated for the 2D array
Matrix configure(Matrix x);                 // Returns a Matrix variable by initializing its members with the user-specified parameters
void print(Matrix x);                       // Prints the 2D array elements held in a Matrix variable
float determinant(Matrix x);                // Returns the calculated determinant value for a Matrix variable
Matrix cofactor(Matrix x);                  // Returns a Matrix variable that holds the cofactors of the 2D array of a Matrix variable
Matrix transpose(Matrix x);                 // Returns a Matrix variable that holds the transpose of the 2D array of a Matrix variable
Matrix inverse(Matrix x);                   // Returns a Matrix variable that holds the inverse of the 2D array of a Matrix variable
Matrix multiply(Matrix x, Matrix y);        // Returns a Matrix variable that holds the multiplication of 2D array of 2 Matrix variables
Matrix solve(Matrix A, Matrix b);           // Returns the solution in Matrix type for a linear equations system specified by 2 Matrix variables



/*** MAIN FUNCTION ***/
// Define and solve a linear equations system (Ax=b)
int main() {
    Matrix A = configure((Matrix){0});                                        // Define and configure the coefficient matrix (n-by-n)
    Matrix b = configure((Matrix){0});                                        // Define and configure the constant (right-hand side) vector (n-by-1)
    Matrix x = solve(A, b);                                                   // Solve the system and store the solution in the variable vector (n-by-1)
    getchar();                                                                // Press enter to see the solution
    for(int i = 0; i < x.n; i++) printf("x%d = %g\t", i+1, x.M[i][0]);        // Print the solution
    putchar('\n');                                                            // Give a new line after printing the solution
    getchar();                                                                // Press enter to exit the program
    destruct(A);
    destruct(b);
    destruct(x);                                                               
    return 0;
}



/*** FUNCTION IMPLEMENTATIONS ***/

/* Returns a Matrix variable by allocating memory for n-by-m sized 2D array */
Matrix construct(Matrix x, int n, int m) {
    x.n = n; // Initialize the number of rows
    x.m = m; // Initialize the number of columns
    x.M = (float **)malloc(n * sizeof(float *)); // Allocate memory for rows
    for (int i = 0; i < n; i++) {
        x.M[i] = (float *)malloc(m * sizeof(float)); // Allocate memory for columns
    }
    return x; // Return the constructed Matrix
}

/* Returns a Matrix variable by deallocating the memory allocated for the 2D array */
Matrix destruct(Matrix x) {
    for (int i = 0; i < x.n; i++) {
        free(x.M[i]); // Deallocate the memory for columns
    }
    free(x.M); // Deallocate the memory for rows
    return (Matrix){0}; // Return an empty Matrix
}


/* Returns a Matrix variable by initializing its members with the user-specified parameters */
Matrix configure(Matrix x) {
    printf("Enter the matrix size as numbers of rows and columns: "); // Echo instructions
    scanf("%d %d", &x.n, &x.m); // Read numbers of rows and columns
    x = construct(x, x.n, x.m); // Call construct() with read values
    printf("Enter the matrix elements row by row:\n"); // Echo instructions
    for (int i = 0; i < x.n; i++) {
        for (int j = 0; j < x.m; j++) {
            scanf("%f", &x.M[i][j]); // Read matrix elements
        }
    }
    printf("\n"); // Put a new line character
    return x; // Return the configured Matrix
}

/* Prints the 2D array elements held in a Matrix variable */
void print(Matrix x) {
    printf("\n"); // Put a new line character
    for (int i = 0; i < x.n; i++) {
        for (int j = 0; j < x.m; j++) {
            printf("%g\t", x.M[i][j]); // Print each element separated by tab characters
        }
        printf("\n"); // Put a new line character after each row
    }
}

/* Returns the calculated determinant value for a Matrix variable */
float determinant(Matrix x) {
    float det = 0; // Define a float variable for the determinant value
    int n = x.n; // Define matrix size with a single variable by assuming a square matrix

    // If the matrix is not square then echo an error and return NAN
    if (x.n != x.m) {
        printf("\nThe matrix is non-square and has no determinant!\n");
        return NAN;
    }

    // If the matrix size is not positive then echo an error and return NAN
    if (n <= 0) {
        printf("\nThe matrix does not exist!\n");
        return NAN;
    }

    // If the matrix has a single element then the determinant value will be equal to its value
    if (n == 1) {
        return x.M[0][0];
    }

    // If the matrix size is 2-by-2 then calculate the determinant by immediate calculation
    if (n == 2) {
        return x.M[0][0] * x.M[1][1] - x.M[0][1] * x.M[1][0];
    }

    // Initialize the determinant to 0 for accumulation
    det = 0;

    // Define a minor matrix by reducing the original matrix size by 1
    Matrix minor = construct((Matrix){0}, n - 1, n - 1);

    // Iterate through the first row elements to create minors
    for (int k = 0; k < n; k++) {
        // Iterate through rows of the original matrix, skipping the first row
        for (int i = 1; i < n; i++) {
            int col = 0; // Define a manual index to iterate the columns of the minor matrix
            // Iterate through columns of the original matrix
            for (int j = 0; j < n; j++) {
                if (j == k) continue; // Skip the column that corresponds to the element the minor is looked for
                minor.M[i - 1][col] = x.M[i][j]; // Assign the original matrix elements to the minor matrix
                col++; // Increment the manual column iterator of the minor matrix
            }
        }
        // Accumulate the minor determinants using recursive call of determinant()
        det += pow(-1, k) * x.M[0][k] * determinant(minor);
    }

    // Clean up the minor matrix
    minor = destruct(minor);

    // Return the calculated determinant value
    return det;
}

/* Returns a Matrix variable that holds the cofactors of the 2D array of a Matrix variable */
Matrix cofactor(Matrix x) {
    int n = x.n; // Define the matrix size with a single variable by assuming a square matrix
    Matrix cfc = construct((Matrix){0}, n, n); // Define and construct a square Matrix variable for the cofactor matrix

    // If the matrix is not square then echo an error and return an empty matrix
    if (x.n != x.m) {
        printf("\nThe matrix is non-square and has no cofactor!\n");
        return (Matrix){0};
    }

    // If the matrix size is not positive then echo an error and return an empty matrix
    if (n <= 0) {
        printf("\nThe matrix does not exist!\n");
        return (Matrix){0};
    }

    // If the matrix has a single element then its cofactor becomes 1
    if (n == 1) {
        cfc.M[0][0] = 1;
        return cfc;
    }

    // Define a minor matrix by reducing the original matrix size by 1
    Matrix minor = construct((Matrix){0}, n - 1, n - 1);

    // Calculate minor matrices for each element of the original matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int minorRow = 0; // Define a manual index for minor matrix rows
            // Iterate through rows of the original matrix
            for (int row = 0; row < n; row++) {
                if (row == i) continue; // Skip the row of the element
                int minorCol = 0; // Define a manual index for minor matrix columns
                // Iterate through columns of the original matrix
                for (int col = 0; col < n; col++) {
                    if (col == j) continue; // Skip the column of the element
                    minor.M[minorRow][minorCol] = x.M[row][col]; // Assign the original matrix elements to the minor matrix
                    minorCol++; // Increment the minor column index
                }
                minorRow++; // Increment the minor row index
            }
            // Calculate the cofactor element value
            cfc.M[i][j] = pow(-1, i + j) * determinant(minor);
        }
    }

    // Clean up the minor matrix
    minor = destruct(minor);

    // Return the cofactor matrix
    return cfc;
}

/* Returns a Matrix variable that holds the transpose of the 2D array of a Matrix variable */
Matrix transpose(Matrix x) {
    Matrix tr = construct((Matrix){0}, x.m, x.n); // Define and construct a reversed-sized (columns-by-rows) Matrix variable for transpose
    // Assign rows of the original matrix to columns of the transpose matrix
    for (int i = 0; i < x.n; i++) {
        for (int j = 0; j < x.m; j++) {
            tr.M[j][i] = x.M[i][j];
        }
    }
    return tr; // Return the transposed Matrix variable
}

/* Returns a Matrix variable that holds the inverse of the 2D array of a Matrix variable */
Matrix inverse(Matrix x) {
    int n = x.n; // Define the matrix size with a single variable by assuming a square matrix
    float det = determinant(x); // Calculate the determinant value of the matrix
    Matrix inv = construct((Matrix){0}, n, n); // Define and construct a square Matrix variable for the inverse matrix

    // If the matrix determinant is 0, then echo an error
    if (det == 0) {
        printf("\nThe matrix is non-square or singular and has no inverse!\n");
        return (Matrix){0};
    }

    // Define and obtain the adjoint matrix as the transpose of the cofactor matrix
    Matrix adj = transpose(cofactor(x));

    // Calculate the inverse matrix elements
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv.M[i][j] = adj.M[i][j] / det;
        }
    }

    // Adjoint matrix may be destructed after all the inverse matrix elements calculated
    adj = destruct(adj);

    // Return the inverted Matrix variable
    return inv;
}

/* Returns a Matrix variable that holds the multiplication of 2D arrays of 2 Matrix variables */
Matrix multiply(Matrix x, Matrix y) {
    Matrix mlt = construct((Matrix){0}, x.n, y.m); // Define and construct a proper-sized Matrix variable for the matrix multiplication

    // If sizes of the matrices match for the matrix multiplication
    if (x.m == y.n) {
        // Iterate through the rows and columns of the multiplication matrix
        for (int i = 0; i < x.n; i++) {
            for (int j = 0; j < y.m; j++) {
                mlt.M[i][j] = 0; // Initialize each element of the multiplication matrix as 0
                for (int k = 0; k < x.m; k++) {
                    mlt.M[i][j] += x.M[i][k] * y.M[k][j]; // Sum of products
                }
            }
        }
    } else {
        // If sizes of matrices do not match then echo an error
        printf("\nSizes of matrices do not match for the matrix multiplication!\n");
    }
    return mlt; // Return the Matrix variable that holds the multiplication of matrices
}

/* Returns the solution in Matrix type for a linear equations system specified by 2 Matrix variables */
Matrix solve(Matrix A, Matrix b) {
    // Return the calculated solution for the linear equations system in Matrix type via multiply() and inverse()
    return multiply(inverse(A), b);
}
