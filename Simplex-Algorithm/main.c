#include <stdio.h>
#include <stdlib.h>

#include "simplex_algorithm.h"

#define EQUALITY_FORM   1
#define INEQUALITY_FORM 0

int main(int argc, const char *argv[]){
    
    LP *P;

    int form = INEQUALITY_FORM;
    int pivot_rule = LCR;
    if (argc < 2 || argc > 5) {
        fprintf(stderr, "Usage: %s <lp file> <form> <pivot_rule> \n", argv[0]);
        return EXIT_SUCCESS;
    }
    if (argc == 2) {
        P = get_LP(argv[1]);
    }
    if (argc == 3) {
        P = get_LP(argv[1]);
        form = atoi(argv[2]);
    }
    if (argc == 4) {
        P = get_LP(argv[1]);
        form = atoi(argv[2]);   
        pivot_rule = atoi(argv[3]);
    }

    // If LP was given in inequality form, add slack variables
    int orig_size;
    if (form == INEQUALITY_FORM) {
        orig_size = P->A->size_c;
        transform_LP_to_equality(P);
    }

    // Solve the LP using the simplex algorithm
    Vector *sol = zero_vector(P->A->size_c);
    int result = simplex_solve_LP(P, pivot_rule, sol);

    // Don't display the slack variables added:
    if (form == INEQUALITY_FORM) {
        Vector *temp = sol;
        sol = zero_vector(orig_size);
        unsigned int i;
        for (i = 0; i < orig_size; i++)
            sol->entries[i] = temp->entries[i];
        free_vector(temp);
    }

    // Display result
    switch (result) {
        case 0:
            print_vector(sol);
            break;
        case 1:
            printf("LP is infeasible.\n");
            break;
        case 2:
            printf("LP is unbounded.\n");
            break;
        case 3:
            printf("LP is of wrong form (m > n).\nMaybe use inequalities?\n");
            break;
    }

    free_LP(P);
    free_vector(sol);
    return EXIT_SUCCESS;
}