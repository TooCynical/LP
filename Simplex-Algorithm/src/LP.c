#include "LP.h"

void free_LP(LP* P) {
	int i;
	for (i = 0; i< (*P).m; i++)
		free((*P).A[i]);
	free((*P).A);	
	free((*P).c);
	free((*P).b);
	free(P);
}

void print_LP(LP *P) {
	int i, j;

	printf("c = [");
      	for (j = 0; j < (*P).n; j++) {
        	printf("%.1lf%s", (*P).c[j], (j == (*P).n-1 ? "]\n" : ", "));
    	}

    printf("A has %d row%s and %d column%s.\n",
            (*P).m, ((*P).m == 1 ? "" : "s"),
            (*P).n, ((*P).n == 1 ? "" : "s"));
    if ((*P).m) {
	    printf("transpose(b) = [");
	    for (i = 0; i < (*P).m; i++) {
	        printf("%.1lf%s", (*P).b[i], (i == (*P).m-1 ? "]\n" : ", "));
	    }
	}
	else
		printf("b is empty.\n");

    for (i = 0; i < (*P).m; i++) {
    	if ((*P).n > 0)
   			printf("A[%d] = [", i);
     	for (j = 0; j < (*P).n; j++) {
        	printf("%5.1lf%s", (*P).A[i][j], (j == (*P).n-1 ? "]\n" : ", "));
     	}
    }
}

LP *get_LP(const char *filename) {
	LP *P = calloc(1, sizeof(LP));
	read_LP(filename, &((*P).m),
			&(*P).n, &(*P).A,
			&(*P).b, &(*P).c);
	return P;
}