#include <stdio.h>
#include <stdlib.h>

#include "LP_reader.h"
#include "lin_alg.h"
#include "error.h"

typedef struct {
    Vector *b; 	// inequality vector
    Vector *c; 	// costs vector
    Matrix *A; 	// inequality matrix

} LP;

typedef struct {
	Vector B;
	Vector N;
	Vector p;	
	Vector r;
	Matrix Q;
	double z0;
} Tableaux;

LP *get_LP(const char*);
void print_LP(LP *P);
void free_LP(LP *P);