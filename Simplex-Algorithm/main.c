#include <stdio.h>
#include <stdlib.h>

#include "LP.h"

int main(int argc, const char *argv[]){
    /* Check if enough arguments were given */
    if (argc < 2) {
        fprintf(stderr, "Usage:  %s  <lp file> [verbose]\n", argv[0]);
        return EXIT_FAILURE;
    }
    /* Read LP from the file and print */
    LP *P = get_LP(argv[1]);
    print_LP(P);
    free_LP(P);
    return EXIT_SUCCESS;
}