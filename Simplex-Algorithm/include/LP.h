#include <stdio.h>
#include <stdlib.h>

#include "LP_reader.h"

typedef struct {
	int m;
	int n;
    double **A;
    double *b;
    double *c;
} LP;

LP *get_LP(const char*);
void print_LP(LP*);
void free_LP(LP*);