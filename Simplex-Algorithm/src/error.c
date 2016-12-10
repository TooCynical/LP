#include "error.h"

void runtime_error(char errorstr[]) {
    printf("ERROR: %s\n", errorstr);
    exit(1);
}