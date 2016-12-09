#include "lin_alg.h"

int solve_system(const Matrix *A, const Vector *b, Vector *ret) {       
    if (A->size_r != A->size_c)
        runtime_error("solve_system: matrix should be square");
    if (A->size_r != b->size || A->size_c != ret->size)
        runtime_error("solve_system: incompatible sizes");
    if (rank(A) < A->size_c)
        runtime_error("solve_system: matrix should be of full rank");



    // Find a QR-decomposition of A and multiply Q^t b
    Matrix *Q = zero_matrix(A->size_r, A->size_c);
    Matrix *R = zero_matrix(A->size_c, A->size_c);
    QR_decomp(A, Q, R);
    Matrix *Qt = trans_matrix(Q);
    Vector *Qtb = mult_vector(Qt, b);
    free_matrix(Qt);

    // Find a solution by back-substitution
    double val;
    int i; // current column
    int j; // current row
    for (i = ret->size - 1; i >= 0; i--) {
        val = Qtb->entries[i];
        for (j = ret->size - 1; j > i; j--) {
            val -= ret->entries[j] * R->entries[i][j];    
        }
        if (val == 0) {
            ret->entries[i] = 0;
        }
        else {
            // No solution
            if (is_zero(R->entries[i][i])) {
                free_matrix(Q);
                free_matrix(R);
                free_vector(Qtb);
                return 0;
            }
            ret->entries[i] = val / R->entries[i][i];
        }
    }
    free_matrix(Q);
    free_matrix(R);
    free_vector(Qtb);
    return 1;
}

void QR_decomp(const Matrix *A, Matrix *Q, Matrix *R) {
    if (A->size_r != Q->size_r || A->size_c != Q->size_c)
        runtime_error("QR_decomp: A, Q are of unequal size");
    if (A->size_c != R->size_c || A->size_c != R->size_c)
        runtime_error("QR_decomp: R of wrong size (should be nxn)");
    
    Vector *q = zero_vector(A->size_r);
    Vector *dummy = zero_vector(A->size_r);
    unsigned int i; // current column
    unsigned int j; // current row
    for (i = 0; i < A->size_c; i++) {
        copy_col(A, q, i);
        for (j = 0; j < i; j++) {
            copy_col(Q, dummy, j);
            R->entries[j][i] = inner_product(dummy, q); 
            scalar_to_vector(dummy, R->entries[j][i]);
            sub_to_vector(q, dummy);
        }
        R->entries[i][i] = norm(q);
        if (!is_zero(R->entries[i][i])) {
            scalar_to_vector(q, 1. / R->entries[i][i]);
            replace_col(Q, q, i);
        }
        else 
            R->entries[i][i] = 0;
    }
    free_vector(q);
    free_vector(dummy);
}

int is_zero(double l) {
    return (abs(l) < ZERO_TOL);
}

int rank(Matrix *A) {
    if (A->size_r != A->size_c)
        runtime_error("rank: matrix should be square");
    
    // Find a QR-decomposition of A
    Matrix *Q = zero_matrix(A->size_r, A->size_c);
    Matrix *R = zero_matrix(A->size_c, A->size_c);
    QR_decomp(A, Q, R);
    
    unsigned int i, r;
    r = 0;
    // Count number of non-zero diagonal elements of R
    for (i = 0; i < R->size_c; i++)
        r += (!is_zero(R->entries[i][i]));
    free_matrix(Q);
    free_matrix(R);
    return r;
}

void reset_vector(const Vector *vector) {
    unsigned int i;
    for (i = 0; i < vector->size; i++)
        vector->entries[i] = 0;
}

double inner_product(const Vector *a, const Vector *b) {
    if (a->size != b->size)
        runtime_error("inner_product: vertices are of unequal size");

    unsigned int i;
    double res = 0;
    for (i = 0; i < a->size; i++) 
        res += a->entries[i] * b->entries[i];
    return res;
}

double norm(const Vector *vector) {
    return sqrt(inner_product(vector, vector));
}

void copy_col(const Matrix *matrix, Vector *vector, unsigned int col) {
    if (vector->size != matrix->size_r)
        runtime_error("copy_col: vertex and matrix size incompatible");
    if (matrix->size_c < col)
        runtime_error("copy_col: col index too high");

    unsigned int i;
    for (i = 0; i < vector->size; i++)
        vector->entries[i] = matrix->entries[i][col];    
}

void replace_col(Matrix *matrix, const Vector *vector, unsigned int col) {
    if (vector->size != matrix->size_r)
        runtime_error("replace_col: vertex and matrix size incompatible");
    if (matrix->size_c < col)
        runtime_error("replace_col: col index too high");

    unsigned int i;
    for (i = 0; i < vector->size; i++)
        matrix->entries[i][col] = vector->entries[i];
}

Matrix *trans_matrix(const Matrix *matrix) {
    Matrix *res = zero_matrix(matrix->size_c, matrix->size_r);
    unsigned int i, j;
    for (i = 0; i < matrix->size_r; i++) {
        for (j = 0; j < matrix->size_c; j++) {
            res->entries[j][i] = matrix->entries[i][j];
        }    
    }
    return res;
}


Matrix *mult_matrix(const Matrix *A, const Matrix *B) {
    if (A->size_c != B->size_r)
        runtime_error("mult_matrix: matrices of incompatible sizes");

    Matrix *res = zero_matrix(A->size_r, B->size_c);
    unsigned int i, j, k;
    for (i = 0; i < A->size_r; i++) {
        for (j = 0; j < B->size_c; j++) {
            for (k = 0; k < A->size_c; k++) {
                res->entries[i][j] += A->entries[i][k] * B->entries[k][j];    
            }
        }    
    }
    return res;
}

void add_to_vector(Vector *a, const Vector *b) {
    if (a->size != b->size)
        runtime_error("sub_vector: vertices are of unequal size");

    unsigned int i;
    for (i = 0; i < a->size; i++)
        a->entries[i] += b->entries[i];
}

void sub_to_vector(Vector *a, const Vector *b) {
    if (a->size != b->size)
        runtime_error("sub_vector: vertices are of unequal size");

    unsigned int i;
    for (i = 0; i < a->size; i++)
        a->entries[i] -= b->entries[i];
}

void scalar_to_vector(Vector *a, const double l) {
    unsigned int i;
    for (i = 0; i < a->size; i++)
        a->entries[i] *= l;    
}

Vector *mult_vector(const Matrix *matrix, const Vector *vector) {
    if (matrix->size_c != vector->size)
        runtime_error("mult_vector: vertex and matrix size incompatible");

    Vector *res = zero_vector(matrix->size_r);
    
    unsigned int i,j;
    for (i = 0; i < matrix->size_r; i++) {
        for (j = 0; j < vector->size; j ++) {
            res->entries[i] += matrix->entries[i][j] * vector->entries[j];
        }
    }
    return res;
}

Vector *add_vector(const Vector *a, const Vector *b) {
    if (a->size != b->size)
        runtime_error("sub_vector: vertices are of unequal size");

    Vector *res = zero_vector(a->size);

    unsigned int i;
    for (i = 0; i < a->size; i++)
        res->entries[i] = a->entries[i] + b->entries[i];
    return res;
}

Vector *sub_vector(const Vector *a, const Vector *b) {
    if (a->size != b->size)
        runtime_error("sub_vector: vertices are of unequal size");

    Vector *res = zero_vector(a->size);

    unsigned int i;
    for (i = 0; i < a->size; i++)
        res->entries[i] = a->entries[i] - b->entries[i];
    return res;
}

Vector *zero_vector(int size) {
    Vector *res = calloc(1, sizeof(Vector));
    res->size = size;
    res->entries = calloc(size, sizeof(double));
    return res;
}

Matrix *zero_matrix(int size_r, int size_c) {
    Matrix *res = calloc(1, sizeof(Matrix));
    res->size_r = size_r;
    res->size_c = size_c;
    res->entries = calloc(size_r, sizeof(double*));
    unsigned int i;
    for (i = 0; i < size_r; i++) 
        res->entries[i] = calloc(size_c, sizeof(double));
    return res;
}

void print_matrix(const Matrix *matrix) {
    unsigned int i,j;
    for (i = 0; i < matrix->size_r; i++) {
        printf("[");
        for (j = 0; j < matrix->size_c; j++) {
            printf("%5.1lf%s", matrix->entries[i][j], (j == matrix->size_c-1 ? "]\n" : ", "));
        }
    }
}

void print_vector(const Vector *vector) {
    unsigned int j;
    printf("[");
    for (j = 0; j < vector->size; j++) {
         printf("%.1lf%s", vector->entries[j], (j == vector->size-1 ? "]\n" : ", "));
    }
}   

void free_vector(Vector *vector) {
    free(vector->entries);
    free(vector);
}

void free_matrix(Matrix *matrix) {
    unsigned int i;
    for (i = 0; i < matrix->size_r; i++)
        free(matrix->entries[i]);
    free(matrix->entries);
    free(matrix);
}