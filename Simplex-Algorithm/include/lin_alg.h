#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error.h"

#define ZERO_TOL 0.001

// Structure for a matrix
typedef struct {
    unsigned int size_r;    // Number of rows of matrix
    unsigned int size_c;    // Number of columns of matrix
    double **entries;       // Actual matrix entries
    unsigned int rank;      // Rank of the matrix
} Matrix;

// Structure for a vector
typedef struct {
    unsigned int size;  // Number of entries
    double *entries;    // Actual vector entries
} Vector;

// Return whether l is very close to zero (see ZERO_TOL)
int is_zero(double l);

// Compute a QR-decompostion of A and store it in Q, R. 
// Uses classical Gram-Schmidt and only gives an orthogonal Q
// if A is of full rank.
void QR_decomp(const Matrix *A, Matrix *Q, Matrix *R);

// Return a pointer to a vector of given size with zero entries
Vector *zero_vector(int size);
// Return a pointer to a matrix of given size with zero entries
Matrix *zero_matrix(int size_r, int size_c);
// Return a pointer to a copy of a vector
Vector *copy_vector(const Vector *vector);
// Copy entries of a vector into another vector
void copy_to_vector(Vector *vector, const Vector *ret);

// Return a pointer to vector whose entries are zero iff
// the corresponding entries in the given vector are
Vector *swap_zero_nonzero(const Vector *vector);

// Return a pointer to a submatrix of matrix consisting of the columns
// corresponding to the entries in bitmask unequal to 0
Matrix *subind_matrix(const Matrix *matrix, const Vector *bitmask);
// Return a pointer to a submatrix of matrix consisting of the columns
// corresponding to the entries in bitmask unequal to 0
Vector *subind_vector(const Vector *vector, const Vector *bitmask);

// Return a pointer to a matrix equal to AB, uses a naieve algorithm
Matrix *mult_matrix(const Matrix *A, const Matrix *B);
// Return a pointer to the transpose of the given matrix
Matrix *trans_matrix(const Matrix *matrix);

//Return whether given vector <= 0
int is_smaller_zero_vect(const Vector *vect);

// Set A equal to l*A
void scalar_to_matrix(Matrix *A, const double l);

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

// Return a pointer to the inverse of A
Matrix *inverse_matrix(const Matrix *A);
// Return solution to Ax = b
Vector *solve_system(const Matrix *A, const Vector *b);
// Return solution to QRx = b, where Q unitary, R upper triangular. 
Vector *solve_system_QR(const Matrix *Q, const Matrix *R, const Vector *b);
// Return the rank of a square matrix
int rank(const Matrix *A);
// Return the rank of an upper triangular matrix
int rank_R(const Matrix *R);

// Free structures
void free_vector(Vector *vector);
void free_matrix(Matrix *matrix);
