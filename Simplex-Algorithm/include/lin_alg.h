#include <stdio.h>
#include <stdlib.h>

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

void print_matrix(Matrix matrix);
void print_vector(Vector vector);

int solve_system(const Matrix matrix, const Vector vect, Vector *ret);