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
    LP *P;      // Underlying LP 
	
    Vector *B; // bitmask indicating which variables are in basis (1 = in basis)
	Vector *N; // bitmask indicating which variables are not in basis (1 = not in basis)
    Vector *p; // As in lecture notes
	Vector *r; // As in lecture notes
	Matrix *Q; // As in lecture notes
    double z0; // As in lecture notes
    
    int size;   // Amount of variables
    int size_B; // Amount of basic variables
    int size_N; // Amount of non-basic variables

    Vector *indices_B; // Conversion table to real indices
    Vector *indices_N; // Conversion table to real indices

    Vector *x; // Basic solution

    int initialized; // Is this tableaux already initialized?
} Tableaux;

// Return pointer to a tableaux with underlying LP P computed for a given basis
Tableaux *build_tableaux(LP *P, Vector *basis);

// Set all values in a tableaux for a given basis, which is assumed to be feasible
void set_tableaux(Tableaux *T, Vector *basis);

// Print relevant info of tableaux to stdout
void print_tableaux(Tableaux *T);

// Transform LP to equality form by adding slack variables and
// multiplying any rows where b < 0 by -1
void transform_LP_to_equality(LP *P);

// Return pointer to LP read from file
LP *get_LP(const char *filename);

// Print relevant info of LP to stdout
void print_LP(LP *P);

// Free only those structures redefined in set_partition
void partial_free_tableaux(Tableaux *T);

// Free structures
void free_LP(LP *P);
void free_tableaux(Tableaux *T);
