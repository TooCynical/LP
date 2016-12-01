/* 	Lucas Slot - lfh.slot@gmail.com
 *  Antonio de la Torre - delatvan@gmail.com
 *	University of Bonn
 *
 *	fm.h
 * 	
 * 	Header for fm.c
 */

#include <stdio.h>
#include <stdlib.h>

#include "LP_reader.h"

/* Basic struct for LPs. */
typedef struct {
	
	/* Normal LP data */
	int m;
	int n;
    double **A;
    double *b;
    double *c;

    /* Matrix used to keep track of the linear
     * combinations used to find new inequalities.
     * Needed to find certificates */
    double **C;
	int orig_m;

	/* Amount of positive, zero, and negative 
	 * coefficients in A corresponding to the
	 * first variable */
    int n_pos;
    int n_zero;
    int n_neg;
} LP;

LP *eliminate_variable(LP*);
LP *get_LP(const char*);

void build_LP(LP*);
void print_LP(LP*);
void count_types_LP(LP*);
void normalize_LP(LP*);
void print_solution(double*, int);
void print_certificate(LP*, int);
void check_feasibility(LP*, int);

double *find_solution(LP**);
