#include <stdio.h>
#include <stdlib.h>

#include "LP_reader.h"
#include "lin_alg.h"
#include "error.h"

// Structure for an LP
typedef struct {
    Vector *b; 	// inequality vector
    Vector *c; 	// costs vector
    Matrix *A; 	// inequality matrix

} LP;

// Structure for a simplex-tableaux
typedef struct {
	Vector *B;
	Vector *N;
	Vector *p;
	Vector *r;
	Matrix *Q;
	double z0;
} Tableaux;

// Return pointer to LP read from file
LP *get_LP(const char *filename);

// Print relevant info of LP to stdout
void print_LP(LP *P);

// Free structures
void free_LP(LP *P);