#include <stdio.h>
#include <stdlib.h>

#include "antenna.h"


int main(int argc, char * argv[])
{
        if (argc < 2) {
                fprintf(stderr, "error: missing input file\n");
                exit(EXIT_FAILURE);
        }

        struct grand_antenna * antenna = NULL;
        int i;
        for (i = 0; i < 10; i++) {
                if (antenna != NULL) free(antenna);
                grand_antenna_load(&antenna, argv[1]);
        }

        double complex z;
        if (antenna->impedance(antenna, GRAND_ANTENNA_ARM_SN, 100, &z) !=
            EXIT_SUCCESS) {
                fprintf(stderr, "error: %s\n", grand_error_get());
                exit(EXIT_FAILURE);
        };

        printf("Z = (%.3E, %.3E) Ohm\n", creal(z), cimag(z));

        if (antenna != NULL) antenna->destroy(&antenna);
        exit(EXIT_SUCCESS);
}
