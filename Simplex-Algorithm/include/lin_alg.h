#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error.h"

typedef struct {
    unsigned int size_r;    // Number of rows of matrix
    unsigned int size_c;    // Number of columns of matrix
    double **entries;       // Actual matrix entries
    unsigned int rank;      // Rank of the matrix
} Matrix;

typedef struct {
    unsigned int size;  // Number of entries
    double* entries;    // Actual vector entries
} Vector;

// Compute a QR-decompostion of A and store it in Q, R. 
// Uses classical Gram-Schmidt
void QR_decomp(const Matrix *A, Matrix *Q, Matrix *R);

// Return a pointer to a vector of given size with zero entries
Vector *zero_vector(int size);
// Return a pointer to a matrix of given size with zero entries
Matrix *zero_matrix(int size_r, int size_c);

// Return a pointer to a matrix equal to AB, uses a naieve algorithm
Matrix *mult_matrix(const Matrix *A, const Matrix *B);
// Return a pointer to the transpose of the given matrix
Matrix *trans_matrix(const Matrix *matrix);

// Return a pointer to matrix * vector
Vector *mult_vector(const Matrix *matrix, const Vector *vector);

// Return a pointer to a vector equal to a+b
Vector *add_vector(const Vector *a, const Vector *b);
// Return a pointer to a vector equal to a-b
Vector *sub_vector(const Vector *a, const Vector *b);
// Set a equal to a+b
void add_to_vector(Vector *a, const Vector *b);
// Set a equal to a-b
void sub_to_vector(Vector *a, const Vector *b);
// Set a equal to l*a
void scalar_to_vector(Vector *a, const double l);


// Replace the i-th column of given matrix by given vector
void replace_col(Matrix *matrix, const Vector *vector, unsigned int col);
// Place the entries of a column in a matrix in a vector
void copy_col(const Matrix *matrix, Vector *vector, unsigned int col);

// Set all entries in given vector to zero
void reset_vector(const Vector *vector);

// Print a matrix to stdout
void print_matrix(const Matrix *matrix);
// Print a vector to stdout
void print_vector(const Vector *vector);

// Return the inner product of two vectors
double inner_product(const Vector *a, const Vector *b);
// Return the norm of a vector
double norm(const Vector *vector);

// Attempt to solve a linear system and store solution in ret.
int solve_system(const Matrix *matrix, const Vector *vect, Vector *ret);

// Free structures
void free_vector(Vector *vector);
void free_matrix(Matrix *matrix);
