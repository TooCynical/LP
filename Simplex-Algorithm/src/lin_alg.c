#include "lin_alg.h"

int solve_system(const Matrix matrix, const Vector vect, Vector *ret) {    	
    
	// If Ax = y has solution, put it in ret and
	return 0;
	// If Ax = y has no solution
	return 1;
}

void print_matrix(const Matrix matrix) {
    unsigned int i,j;
    for (i = 0; i < matrix.size_r; i++) {
        printf("A[%d] = [", i);
        for (j = 0; j < matrix.size_c; j++) {
            printf("%5.1lf%s", matrix.entries[i][j], (j == matrix.size_c-1 ? "]\n" : ", "));
        }
    }
}

void print_vector(const Vector vector) {
    unsigned int j;
    printf("[");
    for (j = 0; j < vector.size; j++) {
         printf("%.1lf%s", vector.entries[j], (j == vector.size-1 ? "]\n" : ", "));
    }
}   