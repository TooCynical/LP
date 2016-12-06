#include "LP.h"

LP *empty_LP() {
    LP *P = calloc(1, sizeof(LP));
    return P;
}

LP *get_LP(const char *filename) {
	LP *P = empty_LP();
    int n,m;
	read_LP(filename, &m, &n,
            &P->A.entries, 
            &P->b.entries, 
            &P->c.entries);
    P->A.size_r = m;
    P->A.size_c = n;
    P->b.size = m;
    P->c.size = n;

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