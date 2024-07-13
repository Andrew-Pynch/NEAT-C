#include "common.h"
#include "game.h"
#include <stdlib.h>

void safe_free(void **pp) {
    if (pp != NULL && *pp != NULL) {
        // printf("Freeing memory at address: %p\n", *pp);
        free(*pp);
        *pp = NULL;
    }
}

void log_string(char *message) {
    if (DEBUG) {
        printf("%s\n", message);
    }
}

void log_int(int message) {
    if (DEBUG) {
        printf("%d\n", message);
    }
}

void log_float(float message) {
    if (DEBUG) {
        printf("%f\n", message);
    }
}
