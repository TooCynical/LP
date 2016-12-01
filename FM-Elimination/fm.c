/* 	Lucas Slot - lfh.slot@gmail.com
 *  Antonio de la Torre - delatvan@gmail.com
 *	University of Bonn
 *
 *	fm.c
 * 	
 * 	Contains an implementation of Fourier-Motzkin-
 *  elimination for solving Linear Programs (LPs).
 *  Provided a correctly formatted LP in a file,
 *  this program will determine whether the LP
 *  is feasible, and print either a solution or
 *  a certficitate of infeasibility, accordingly.
 *
 *  WARNING: input should exactly as specified in 
 *  the problem: otherwise GIGO applies.
 */

#include "fm.h"

/* Free memory allocated to an LP,
 * assuming c is not used. */
void free_LP(LP* P) {
	int i, j, k;
	for (i = 0; i< (*P).m; i++) {
		free((*P).A[i]);
		free((*P).C[i]);
	}
	free((*P).A);	
	free((*P).C);
	free((*P).b);
	free(P);
}

/* Free a sequence of LPs. */
void free_sequence(LP **R) {
	int i, k;
	k = (*R[0]).n + 1;
	for (i = 0; i < k; i++)
		free_LP(R[i]);
	free(R);
}

/* Print the word 'empty' followed by a certificate
 * showing the fully reduced LP P is infeasible, 
 * given an i such that b[i] < 0, using the i-th
 * row of C. */
void print_certificate(LP* P, int c) {
	int i;
	printf("empty ");
    printf("[");
	for (i = 0; i < (*P).orig_m; i++)
        printf("%.3lf%s", (*P).C[c][i], (i == ((*P).orig_m-1) ? "]\n" : ", "));
}

/* Print a solution vector provided by find_solution. */
void print_solution(double *sol, int n) {
    int i;
    printf("[");
    for (i = 0; i < n; i++) {
        printf("%.3lf%s", sol[i], (i == n - 1? "]\n" : ", "));
    }
    /* Free sol. */
    free(sol);
}

/* Find a solution given a sequence of LPs by back-substitution. */
double *find_solution(LP **R) {
	int i, j, k, l, n, m, max_set, min_set;
	double min, max;
	double *b;
	LP *P = R[0];
	n = (*P).n;
	double *sol = calloc(n, sizeof(double));

	/* Starting at the second to last LP, backtrack
	 * through the entire sequence R of LPs. */
	for (i = 1; i <= n; i++){
		min = 0;
		max = 0;
		min_set = 0;
		max_set = 0;

		/* Update b in the n-i-th LP by back-substitution. */
		b = calloc((*R[n-i]).m, sizeof(double));
		for (k = 0; k < (*R[n-i]).m; k++)
			b[k] += (*R[n-i]).b[k];

		/* Update the values using previously found solutions 
		 * for the variables x_{n-i+1}, ..., x_n. */
		for (j = 1; j < i; j++) {
			for (k = 0; k < (*R[n-i]).m; k++)
				b[k] -= (*R[n-i]).A[k][j] * sol[n - (i-j)];
		}

		/* Find a value for x_{n-i} using the updated b. */
		for (l = 0; l < (*R[n-i]).m; l++) {
			/* Update min if i-th to last variable negative. */
			if ((*R[n-i]).A[l][0] < 0) {
				if (-b[l] > min || min_set == 0) {
					min = -b[l];
					min_set = 1;
				}
			}
			/* Update max if i-th to last variable positive. */
			if ((*R[n-i]).A[l][0] > 0) {
				if (b[l] < max || max_set == 0) {
					max = b[l];
					max_set = 1;
				}
			}
		}
		/* We have to take some care in case the max/min was never set */
		if (!max_set)
			max = min+1;
		if (!min_set)
			min = max-1;
		sol[n-i] = (max + min) / 2.0;
		free(b);
	}
	return sol;
}

/* Create a new LP, equivalent to the one given,
 * that has one fewer variables using FM-elimination */
LP *eliminate_variable(LP *P) {
	int i, j, k, count;
	count = 0;

	/* Allocate memory */
	LP *new_P = calloc(1, sizeof(LP));
	(*new_P).n = (*P).n - 1;
	(*new_P).m = (*P).n_pos * (*P).n_neg + (*P).n_zero;
	if((*new_P).m < 0 || (*new_P).m > 1000000) {
		printf("Too many inequalities: exiting...\n");
		exit(0);
	}
	(*new_P).orig_m = (*P).orig_m;
	(*new_P).A = calloc((*new_P).m, sizeof(double*));
	(*new_P).C = calloc((*new_P).m, sizeof(double*));
	
	for (i = 0; i < (*new_P).m; i++) {
		(*new_P).A[i] = calloc((*new_P).n, sizeof(double));
		(*new_P).C[i] = calloc((*new_P).orig_m, sizeof(double));
	}
	
	(*new_P).b = calloc((*new_P).m, sizeof(double));

	/* Set new equations */
	for (i = 0; i < (*P).m; i++) {
		for (j = 0; j < (*P).m; j++) {
			/* Positive/negative pair, add equations. */
			if ((*P).A[i][0] > 0 && (*P).A[j][0] < 0) {
				/* Update A and b */
				for (k = 0; k < (*new_P).n; k++) 
					(*new_P).A[count][k] = (*P).A[i][k+1] + (*P).A[j][k+1];
				(*new_P).b[count] = (*P).b[i] + (*P).b[j];
				/* Update C */
				for (k = 0; k < (*new_P).orig_m; k++)
					(*new_P).C[count][k] = (*P).C[i][k] + (*P).C[j][k];
				count ++;
			}
			/* Zero coefficient, just copy equation. */
			if ((*P).A[i][0] == 0 && i == j) {
				/* Update A and b */
				for (k = 0; k < (*new_P).n; k++) 
					(*new_P).A[count][k] = (*P).A[i][k+1];
				(*new_P).b[count] = (*P).b[i];
				/* Update C */
				for (k = 0; k < (*new_P).orig_m; k++)
					(*new_P).C[count][k] = (*P).C[i][k];
				count ++;
			}
		}
	}
	normalize_LP(new_P);
	count_types_LP(new_P);
	return new_P;
}

/* Prints the matrix A, C and vector b of an LP, as well 
 * as the number of positive, zero and negative coefficients
 * corresponding to the first variable */
void print_LP(LP *P) {
	int i, j;

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
    for (i = 0; i < (*P).m; i++) {
    	if ((*P).orig_m > 0)
   			printf("C[%d] = [", i);
     	for (j = 0; j < (*P).orig_m; j++) {
        	printf("%5.1lf%s", (*P).C[i][j], (j == (*P).orig_m-1 ? "]\n" : ", "));
     	}
  	}

    printf("Types: %d - %d - %d\n", (*P).n_pos, (*P).n_zero, (*P).n_neg);
}

/* Check whether an LP is feasible by reducing it to
 * an equivalent system in zero variables and then 
 * checking b >= 0. Depending on feasibility, a
 * solution or a certificate will be printed. */
void check_feasibility(LP *P, int verbose) {
	int i;
	double *sol, cert;

	/* red will contain a sequence of n+1 LPs equivalent
	 * to P, such that red[i] has one fewer variable 
	 * than red[i-1]. */
	LP **red = calloc((*P).n + 1, sizeof(LP*));
	red[0] = P;

	if (verbose) 
		print_LP(red[0]);

	/* Reduce P, save each intermediate LP */
	for (i = 1; i < (*P).n + 1; i++) {
		red[i] = eliminate_variable(red[i-1]);
		if (verbose) 
			print_LP(red[i]);
	}

	/* Check b >= 0 */
	for (i = 0; i < (*red[(*P).n]).m; i++) {
		/* If b[i] < 0, the system is infeasible,
		 * and a certificate is provided by the 
		 * i-th row of C */
		if ((*red[(*P).n]).b[i] < 0) {
			print_certificate(red[(*P).n], i);
			free_sequence(red);
			return;
		}
	}
	/* If b>=0, a solution can be found by back-substitution */
	print_solution(find_solution(red), (*P).n);
	free_sequence(red);
}

/* Count the negative, positive, and zero coefficients
 * in A corresponding to the first variable. */
void count_types_LP(LP *P) {
	int i,j;

    /* Only count types if there are actually any equations */
    if (!(*P).n)
    	return;

	for (i = 0; i < (*P).m; i++) {
		if ((*P).A[i][0] > 0)
			(*P).n_pos ++;
		if ((*P).A[i][0] == 0)
			(*P).n_zero ++;
		if ((*P).A[i][0] < 0)
			(*P).n_neg ++;
	}
}

/* Normalize b and the rows of A with respect
 * to the coefficient of the first variable. */
void normalize_LP(LP *P) {
	int i,j;
   	double t;

 	/* only normalize if there are actually any equations */
    if (!(*P).n)
    	return;

	for (i = 0; i < (*P).m; i++) {
		/* If first coefficient is non-zero, we can normalize. */
		if (t = abs((*P).A[i][0])) {
			(*P).b[i] /= t;
			for (j = 0; j < (*P).n; j++)
				(*P).A[i][j] /= t;
			for (j = 0; j < (*P).orig_m; j++)
				(*P).C[i][j] /= t;
		}
	}
}

/* Read an LP from a file and set up 
 * data structures */
LP *get_LP(const char *filename) {
	int i;
	LP *P = calloc(1, sizeof(LP));
	
	/* Read file */
	read_LP(filename, &((*P).m),
			&(*P).n, &(*P).A,
			&(*P).b, &(*P).c);
	/* We don't need c, free it right away. */
	free((*P).c);

	/* Set C to mxm identity. */
	(*P).orig_m = (*P).m;
	(*P).C = calloc((*P).m, sizeof(double*));
	for (i = 0; i < (*P).m; i++) {
		(*P).C[i] = calloc((*P).m, sizeof(double*));
		(*P).C[i][i] = 1;
	}

	normalize_LP(P);
	count_types_LP(P);
	return P;
}

int main(int argc, const char *argv[]){
	
	/* Check if enough arguments were given */
	if (argc < 2) {
    	fprintf(stderr, "Usage:  %s  <lp file> [verbose]\n", argv[0]);
      	return EXIT_FAILURE;
   	}

   	/* If a third argument is given, 
   	 * print all reduction steps */
   	int verbose;
   	if (argc > 2)
   		verbose = 1;
   	else
   		verbose = 0;

   	/* Read LP from the file and check if feasible */
    LP *P = get_LP(argv[1]);
    check_feasibility(P, verbose);
    return EXIT_SUCCESS;
}