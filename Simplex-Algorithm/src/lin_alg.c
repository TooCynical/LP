#include "lin_alg.h"

Matrix *inverse_matrix(const Matrix *A) {
    if (A->size_r != A->size_c)
        runtime_error("inverse_matrix: matrix should be square");

    // Find a QR-decomposition of A
    Matrix *Q = zero_matrix(A->size_r, A->size_c);
    Matrix *R = zero_matrix(A->size_c, A->size_c);
    QR_decomp(A, Q, R);

    if (rank_R(R) < A->size_c) {
        runtime_error("inverse_matrix: matrix singular");
    }
    Matrix *ret = zero_matrix(A->size_r, A->size_c);

    // Solve Ax = e_i for all basis vectors
    unsigned int i;
    for (i = 0; i < A->size_c; i ++) {
        Vector *e = zero_vector(A->size_c);
        e->entries[i] = 1;
        Vector *sol = solve_system_QR(Q, R, e);
        replace_col(ret, sol, i);
        // This freeing is suboptimal
        free_vector(e);
        free_vector(sol);
    }
    free_matrix(Q);
    free_matrix(R);
    return ret;
}

Vector *solve_system(const Matrix *A, const Vector *b) {       
    if (A->size_r != A->size_c)
        runtime_error("solve_system: matrix should be square");
    if (A->size_r != b->size)
        runtime_error("solve_system: incompatible sizes");

    // Find a QR-decomposition of A and multiply Q^t b
    Matrix *Q = zero_matrix(A->size_r, A->size_c);
    Matrix *R = zero_matrix(A->size_c, A->size_c);
    QR_decomp(A, Q, R);

    if (rank_R(R) < A->size_c)
        runtime_error("solve_system: matrix should be of full rank");

    Vector *ret = solve_system_QR(Q, R, b);
    free_matrix(Q);
    free_matrix(R);
    return ret;
}

Vector *solve_system_QR(const Matrix *Q, const Matrix *R, const Vector *b) {
    
    if (Q->size_c != R->size_r)
        runtime_error("solve_system_QR: Q, R dimension mismatch");

    Matrix *Qt = trans_matrix(Q);
    Vector *Qtb = mult_vector(Qt, b);
    free_matrix(Qt);

    // Find a solution by back-substitution
    Vector *ret = zero_vector(R->size_c);
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
                runtime_error("solve_system_QR: system has no solution");
            }
            ret->entries[i] = val / R->entries[i][i];
        }
    }
    free_vector(Qtb);
    return ret;
}

void QR_decomp(const Matrix *A, Matrix *Q, Matrix *R) {
    if (A->size_r != Q->size_r || A->size_c != Q->size_c)
        runtime_error("QR_decomp: A, Q are of unequal size");
    if (A->size_c != R->size_c || A->size_c != R->size_c)
        runtime_error("QR_decomp: R of wrong size (should be mxm)");
    
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
    if (l >= 0)
        return (l < ZERO_TOL);
    else 
        return (-l < ZERO_TOL);
}

int rank(const Matrix *A) {
    if (A->size_r != A->size_c)
        runtime_error("rank: matrix should be square");
    
    // Find a QR-decomposition of A
    Matrix *Q = zero_matrix(A->size_r, A->size_c);
    Matrix *R = zero_matrix(A->size_c, A->size_c);
    QR_decomp(A, Q, R);
    
    int r = rank_R(R);

    free_matrix(Q);
    free_matrix(R);
    return r;
}

int rank_R(const Matrix *R) {
    unsigned int i, r;
    r = 0;
    // Count number of non-zero diagonal elements of R
    for (i = 0; i < R->size_c; i++)
        r += (!is_zero(R->entries[i][i]));
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

int is_smaller_zero_vect(const Vector *vector) {
    unsigned int i;
    for (i = 0; i < vector->size; i++) {
        if (vector->entries[i] > ZERO_TOL)
            return 0;
    }
    return 1;
}

Vector *subind_vector(const Vector *vector, const Vector *bitmask) {
    if (vector->size != bitmask->size)
        runtime_error("subind_vector: bitmask should be of same size as vector");
    unsigned int i, count;
    count = 0;
    // Set which indices are not zero in bitmask
    Vector *indices = zero_vector(vector->size);
    for (i = 0; i < bitmask->size; i++) {
        if (bitmask->entries[i] != 0) {
            indices->entries[count] = i;
            count ++;
        }
    }

    // Now create the subvector
    Vector *res = zero_vector(count);
    unsigned int j, d;
    for (j = 0; j < count; j++) {
        d = (int)indices->entries[j];
        res->entries[j] = vector->entries[d];
    }
    free_vector(indices);
    return res;
}

Matrix *subind_matrix(const Matrix *matrix, const Vector *bitmask) {
    if (matrix->size_c != bitmask->size)
        runtime_error("subind_matrix: bitmask should have an entry for each column of matrix");
    unsigned int i, count;
    count = 0;
    // Set which indices are not zero in bitmask
    Vector *indices = zero_vector(matrix->size_c);
    for (i = 0; i < bitmask->size; i++) {
        if (bitmask->entries[i] != 0) {
            indices->entries[count] = i;
            count ++;
        }
    }

    // Now create the submatrix
    Matrix *res = zero_matrix(matrix->size_r, count);
    unsigned int j, d;
    for (i = 0; i < matrix->size_r; i++) {
        for (j = 0; j < count; j++) {
            d = (int)indices->entries[j];
            res->entries[i][j] = matrix->entries[i][d];
        }
    }
    free_vector(indices);
    return res;
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

void scalar_to_matrix(Matrix *A, const double l) {
    unsigned int i, j;
    for (i = 0; i < A->size_r; i++) {
        for (j = 0; j < A->size_c; j++) {
            A->entries[i][j] *= l;    
        }
    }
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

Vector *swap_zero_nonzero(const Vector *vector) {
    Vector *ret = zero_vector(vector->size);
    unsigned int i;
    for (i = 0; i < vector->size; i++)
        ret->entries[i] = (vector->entries[i] == 0);
    return ret;
}

void copy_to_vector(Vector *vector, const Vector *ret) {
    if (vector->size != ret->size)
        runtime_error("copy_to_vector: vectors of unequal size");
    unsigned int i;
    for (i = 0; i < vector->size; i++)
        ret->entries[i] = vector->entries[i];
}

Vector *copy_vector(const Vector *vector) {
    Vector *ret = zero_vector(vector->size);
    unsigned int i;
    for (i = 0; i < vector->size; i++)
        ret->entries[i] = vector->entries[i];
    return ret;
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
            printf("%5.2lf%s", matrix->entries[i][j], (j == matrix->size_c-1 ? "]\n" : ", "));
        }
    }
}

void print_vector(const Vector *vector) {
    unsigned int j;
    printf("[");
    for (j = 0; j < vector->size; j++) {
         printf("%.2lf%s", vector->entries[j], (j == vector->size-1 ? "]\n" : ", "));
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