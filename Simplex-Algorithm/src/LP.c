#include "LP.h"

Tableaux *build_tableaux(LP *P, Vector *basis) {
    Tableaux *ret = calloc(1, sizeof(Tableaux));
    ret->P = P;
    set_tableaux(ret, basis);
    ret->initialized = 1;
    return ret;
}

LP *empty_LP() {
    LP *P = calloc(1, sizeof(LP));
    P->A = calloc(1, sizeof(Matrix));
    P->b = calloc(1, sizeof(Vector));
    P->c = calloc(1, sizeof(Vector));
    return P;
}

// Negate a row in an LP - not exported
void negate_row_LP(LP *P, int row) {
    if (row > P->A->size_r)
        runtime_error("negate_row_LP: row number too large");
    unsigned int i;
    for (i = 0; i < P->A->size_c; i++)
        P->A->entries[row][i] *= -1;
    P->b->entries[row] *= -1;
}

void transform_LP_to_equality(LP *P) {
    // Set a new A equal to (A | Im)
    Matrix *new_A = zero_matrix(P->A->size_r, P->A->size_c + P->A->size_r);
    // Copy entries
    unsigned int i, j;
    for (i = 0; i < P->A->size_r; i++) {
        for (j = 0; j < P->A->size_c; j++) {
            new_A->entries[i][j] = P->A->entries[i][j];
        }
    }
    // Add identity entries
    for (i = 0; i < P->A->size_r; i++)
        new_A->entries[i][P->A->size_c + i] = 1;

    // Set a new c equal to (c 0)
    Vector *new_c = zero_vector(P->c->size + P->A->size_r);
    for (i = 0; i < P->c->size; i++)
        new_c->entries[i] = P->c->entries[i];

    free_matrix(P->A);
    free_vector(P->c);
    P->A = new_A;
    P->c = new_c;

    // Negate rows where b < 0
    for (i = 0; i < P->A->size_r; i++) {
        if (P->b->entries[i] < 0)
            negate_row_LP(P, i);
    }
}

LP *get_LP(const char *filename) {
	LP *P = empty_LP();
    int n,m;
	read_LP(filename, &m, &n,
            &P->A->entries, 
            &P->b->entries, 
            &P->c->entries);
    P->A->size_r = m;
    P->A->size_c = n;
    P->b->size = m;
    P->c->size = n;

	return P;
}

// Helper function for set_tableaux - not exported
void set_basic_solution(Tableaux *T) {
    int i, real_index;
    T->x = zero_vector(T->size);
    for (i = 0; i < T->size_B; i++) {
        real_index = (int) T->indices_B->entries[i];
        T->x->entries[real_index] = T->p->entries[i];
        if (T->x->entries[real_index] < -ZERO_TOL)
            runtime_error("set_basic_solution: negative entry in basic solution");
    }
}

// Helper function for set_tableaux - not exported
void set_conversion_tables(Tableaux *T) {
    T->indices_B = zero_vector(T->size_B);
    T->indices_N = zero_vector(T->size_N);
    unsigned int i, count;
    count = 0;
    for (i = 0; i < T->size; i++) {
        if (T->B->entries[i] != 0) {
            T->indices_B->entries[count] = i;
            count ++;
        }
    }
    count = 0;
    for (i = 0; i < T->size; i++) {
        if (T->N->entries[i] != 0) {
            T->indices_N->entries[count] = i;
            count ++;
        }
    }
}

// Helper function for set_tableaux - not exported
void set_r_and_sizes(Tableaux *T) { 
    // Compute helper vectors c_B and c_N
    Vector *cB = subind_vector(T->P->c, T->B);
    Vector *cN = subind_vector(T->P->c, T->N);

    // Compute r = c_N - (c_B^t * Q)^t = Q^t * c_B
    T->r = copy_vector(cN);
    Matrix *Qt = trans_matrix(T->Q);
    Vector *temp = mult_vector(Qt, cB);
    add_to_vector(T->r, temp);
    free_matrix(Qt);
    free_vector(temp);

    // Compute z0
    T->z0 = inner_product(cB, T->p);

    // Set sizes for later referencing
    T->size_B = cB->size;
    T->size_N = cN->size;
    T->size = T->size_B + T->size_N;
    
    // Free structures
    free_vector(cB);
    free_vector(cN);
}

void set_p_and_Q(Tableaux *T) {
    // Compute helper matrices A_B, A_N and (A_B)^-1
    Matrix *A = T->P->A;
    Matrix* AB = subind_matrix(A, T->B);
    Matrix* AN = subind_matrix(A, T->N);
    Matrix* ABinv = inverse_matrix(AB);

    // Compute p and Q
    T->p = mult_vector(ABinv, T->P->b);
    T->Q = mult_matrix(ABinv, AN);
    scalar_to_matrix(T->Q, -1);

    // Free structures 
    free_matrix(AB);
    free_matrix(AN);
    free_matrix(ABinv);
}

void set_tableaux(Tableaux *T, Vector *basis) {
    
    // Release previous structures
    if (T->initialized == 1)
        partial_free_tableaux(T);
    
    // Set basic / non-basic variables
    T->B = copy_vector(basis);
    T->N = swap_zero_nonzero(T->B);

    // Set p and Q
    set_p_and_Q(T);

    // Set r and some sizes for later reference
    set_r_and_sizes(T); 

    // Set index conversion tables
    set_conversion_tables(T);

    // Compute the basic solution
    set_basic_solution(T);
}

void print_tableaux(Tableaux *T) {
    printf("Underlying LP:\n");
    print_LP(T->P);
    
    printf("trans(B):\n");
    print_vector(T->B);
    
    printf("trans(N):\n");
    print_vector(T->N);
    
    printf("trans(p):\n");
    print_vector(T->p);
    
    printf("trans(r):\n");
    print_vector(T->r);
    
    printf("Q:\n");
    print_matrix(T->Q);
    
    printf("z0: \n%f\n", T->z0);

    printf("x:\n");
    print_vector(T->x); 
}

void print_LP(LP *P) {
    printf("A = \n");
    print_matrix(P->A);
    printf("trans(b) = \n");
    print_vector(P->b);
    printf("trans(c) = \n");
    print_vector(P->c);
}

void partial_free_tableaux(Tableaux *T) {
    free_vector(T->B);
    free_vector(T->N);
    free_vector(T->p);
    free_vector(T->r);
    free_matrix(T->Q);
    
    free_vector(T->indices_B);
    free_vector(T->indices_N);

    free_vector(T->x);
}

void free_tableaux(Tableaux *T) {
    free_vector(T->B);
    free_vector(T->N);
    free_vector(T->p);
    free_vector(T->r);
    free_matrix(T->Q);
    
    free_vector(T->indices_B);
    free_vector(T->indices_N);
    
    free_vector(T->x);
    
    free(T);
}

void free_LP(LP *P) {
    free_matrix(P->A);
    free_vector(P->b);
    free_vector(P->c);
    free(P);       
}
