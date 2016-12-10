#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error.h"
#include "LP.h"

#define SOLVABLE    0
#define INFEASIBLE  1
#define UNBOUNDED   2
#define WRONG_FORM  3

#define BLAND 0
#define LCR 1

// Solve an LP using the simplex algorithm, provided an initial feasible basis.
// Return 0 if system feasible and bounded and store optimal solution in ret
// Return 1 if system unbounded
int simplex_solve_LP_basis(LP *P, Vector *basis, int pivot_rule, Vector *ret);

// Different methods for choosing an alpha
int choose_alpha(Tableaux *T, int pivot_rule);
// Bland's Rule
int choose_alpha_BLAND(Tableaux *T);
// Largest Coeffient Rule
int choose_alpha_LCR(Tableaux *T);

// Choose beta as in lecture notes. If no valid beta exists,
// return -1, indicating the LP is unbounded
int choose_beta(Tableaux *T, int alpha);

// Swap alpha for beta in B using the index conversion tables
void swap_basis(Vector* B, const Tableaux *T, int alpha, int beta);

// Return a pointer to the LP used to find an initial basis for the LP
LP *initial_basis_LP(LP *P);

// Find an initial basis for an LP using the method from the lecture notes
// Return 0 if system feasible and store the basis in ret
// Return 1 if no basis can be found
int find_initial_basis(LP *P, int pivot_rule, Vector* ret);

// Solve an LP using the simplex algorithm.
// Return 0 if system feasible and bounded and store optimal solution in ret
// Return 1 if system infeasible
// Return 2 if system unbounded
// Return 3 if system cannot be solved (m > n)
int simplex_solve_LP(LP *P, int pivot_rule, Vector *ret);
