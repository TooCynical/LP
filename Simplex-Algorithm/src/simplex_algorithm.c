#include "simplex_algorithm.h"

int choose_alpha(Tableaux *T, int pivot_rule) {
    switch (pivot_rule) {
        case BLAND:
            return choose_alpha_BLAND(T);
            break;
        case LCR:
            return choose_alpha_LCR(T);
            break;
        default:
            return choose_alpha_BLAND(T);
    }
}

int choose_alpha_BLAND(Tableaux *T) {
    int alpha;
    for (alpha = 0; alpha < T->r->size; alpha++) {
        if (T->r->entries[alpha] > 0)
            return alpha;
    }
    runtime_error("choose_alpha_SC: no valid alpha could be found");
    return -1;
}

int choose_alpha_LCR(Tableaux *T) {
    int alpha;
    int best_alpha;
    double largest_coef = -1;
    for (alpha = 0; alpha < T->r->size; alpha++) {
        if (T->r->entries[alpha] > largest_coef) {
            best_alpha = alpha;
            largest_coef = T->r->entries[alpha];
        }
    }
    if (largest_coef > 0) 
        return best_alpha;
    runtime_error("choose_alpha_SC: no valid alpha could be found");
    return -1;
}

int choose_beta(Tableaux *T, int alpha) {
    int i, best_beta, best_set;
    double best_value, cur_value;

    best_set = 0;
    for (i = 0; i < T->Q->size_r; i++) {
        if (T->Q->entries[i][alpha] < -ZERO_TOL) {
            cur_value = T->p->entries[i] / T->Q->entries[i][alpha];
            if (cur_value > best_value || !best_set) {
                best_beta = i;
                best_value = cur_value;
                best_set = 1;
            }
        }
    }
    if (best_set == 1)
        return best_beta;
    else
        return -1;
}

void swap_basis(Vector* B, const Tableaux *T, int alpha, int beta) {
    B->entries[(int)T->indices_N->entries[alpha]] = 1;
    B->entries[(int)T->indices_B->entries[beta]] = 0;
}

int simplex_solve_LP_basis(LP *P, Vector *basis, int pivot_rule, Vector *ret) {

    // Compute an initial simplex-tableaux
    Tableaux *T = build_tableaux(P, basis);

    // Keep swapping basis elements according to the pivot rule until r <= 0
    while (!is_smaller_zero_vect(T->r)) {
        int alpha = choose_alpha(T, pivot_rule);
        int beta = choose_beta(T, alpha);
        
        // LP unbounded
        if (beta == -1) {
            free_tableaux(T);
            return 1;
        }
        swap_basis(basis, T, alpha, beta);
        set_tableaux(T, basis);
    }

    copy_to_vector(T->x, ret);
    free_tableaux(T);
    return 0;
}

// This should maybe be modularized
LP *initial_basis_LP(LP *P) {
    // Set a new A equal to (A | Im)
    Matrix *new_A = zero_matrix(P->A->size_r, P->A->size_c + P->A->size_r);
    unsigned int i, j;
    // Copy entries
    for (i = 0; i < P->A->size_r; i++) {
        for (j = 0; j < P->A->size_c; j++) {
            new_A->entries[i][j] = P->A->entries[i][j];
        }
    }
    // Add identity entries
    for (i = 0; i < P->A->size_r; i++)
        new_A->entries[i][P->A->size_c + i] = 1;

    // Set a new c equal to (0, 0, ... , 0, -1, ... , -1)
    Vector *new_c = zero_vector(P->c->size + P->A->size_r);
    for (i = P->A->size_c; i < new_c->size; i++)
        new_c->entries[i] = -1;

    LP *new_LP = calloc(1, sizeof(LP));
    new_LP->A = new_A;
    new_LP->c = new_c;
    new_LP->b = copy_vector(P->b);

    return new_LP;
}

int find_initial_basis(LP *P, int pivot_rule, Vector* ret) {
    
    // Find initial basis LP I and set initial basis (for I!)
    LP *I = initial_basis_LP(P);
    Vector *basis = zero_vector(I->A->size_c);
    unsigned int i;
    for (i = P->A->size_c; i < I->A->size_c; i++)
        basis->entries[i] = 1;

    // Find an optimal solution for I
    Vector *sol = zero_vector(I->A->size_c);
    simplex_solve_LP_basis(I, basis, pivot_rule, sol);
    free_vector(basis); 

    // Check if it has nonnegative value
    if (inner_product(I->c, sol) < 0) {
        free_LP(I);
        free_vector(sol);
        return 1;
    }

    free_LP(I);

    // Copy basis
    for (i = 0; i < ret->size; i++)
        ret->entries[i] = sol->entries[i];
    free_vector(sol);
    return 0;
}

int simplex_solve_LP(LP *P, int pivot_rule, Vector *ret) {

    // Check m <= n
    if (P->A->size_r > P->A->size_c)
        return WRONG_FORM;
    // Find an initial basis
    Vector *init_basis = zero_vector(ret->size);

    int init_basis_result = find_initial_basis(P, pivot_rule, init_basis);

    // System infeasible
    if (init_basis_result == 1) {
        free_vector(init_basis);
        return INFEASIBLE;
    }


    // Solve the LP and return whether it is bounded
    if (simplex_solve_LP_basis(P, init_basis, pivot_rule, ret) == 0) {
        free_vector(init_basis);    
        return SOLVABLE;
    }
    else {
        free_vector(init_basis);
        return UNBOUNDED;
    }
}