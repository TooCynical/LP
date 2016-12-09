#include "LP.h"

LP *empty_LP() {
    LP *P = calloc(1, sizeof(LP));
    P->A = calloc(1, sizeof(Matrix));
    P->b = calloc(1, sizeof(Vector));
    P->c = calloc(1, sizeof(Vector));
    return P;
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

void print_LP(LP *P) {
    printf("A = \n");
    print_matrix(P->A);
    printf("trans(b) = \n");
    print_vector(P->b);
    printf("trans(c) = \n");
    print_vector(P->c);
}

void free_LP(LP *P) {
    free_matrix(P->A);
    free_vector(P->b);
    free_vector(P->c);
    free(P);       
}