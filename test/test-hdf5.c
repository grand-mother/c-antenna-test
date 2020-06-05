#include <stdio.h>
#include <stdlib.h>

#include "antenna.h"


int main(int argc, char * argv[])
{
        if (argc < 2) {
                fprintf(stderr, "missing input file\n");
                exit(EXIT_FAILURE);
        }

        struct grand_antenna * antenna;
        grand_antenna_load(&antenna, argv[1]);

        double complex z = antenna->impedance(antenna, GRAND_ANTENNA_ARM_SN,
                                              100.);
        printf("Z = (%.3E, %.3E) Ohm\n", creal(z), cimag(z));

        if (antenna != NULL) antenna->destroy(&antenna);
        exit(EXIT_SUCCESS);
}
