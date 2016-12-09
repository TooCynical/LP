#include <stdio.h>
#include <stdlib.h>

#include "LP.h"

int main(int argc, const char *argv[]){
    
    LP *P;

    // Check if enough arguments were given
    if (argc < 2) {
        //fprintf(stderr, "Usage:  %s  <lp file> [verbose]\n", argv[0]);
        //return EXIT_FAILURE;
        P = get_LP("./instances/test_1");
    }
    else {
        /* Read LP from the file and print */
        P = get_LP(argv[1]);
    }
    print_LP(P);
    // print_vector(add_vector(&P->b, &P->b));
    // print_vector(sub_vector(&P->b, &P->b));
    // printf("norm(b) : %f\n", norm(&P->b));
    // replace_col(&P->A, &P->b, 0);
    // print_matrix(&P->A);
    // print_matrix(trans_matrix(&P->A));
    // Matrix *Q = zero_matrix(P->A.size_r, P->A.size_c);
    // Matrix *R = zero_matrix(P->A.size_c, P->A.size_c);
    // QR_decomp(&P->A, Q, R);

    // print_matrix(Q);
    // print_matrix(R);
    // printf("A:\n");
    // print_matrix(mult_matrix(Q, R));

    Vector *sol = zero_vector(P->A->size_c);
    Vector *bb = zero_vector(P->A->size_r);
    bb->entries[0] = 2;
    
    printf("Rank A: %d\n", rank(P->A));
    print_vector(bb);

    if (solve_system(P->A, bb, sol)) {
        print_vector(sol);
    }
    else
        printf("No solution\n");

    free_vector(sol);
    free_vector(bb);
    free_LP(P);


    return EXIT_SUCCESS;
}