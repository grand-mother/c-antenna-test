#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

        double z_mag,z_phase;
        if (antenna->impedance(antenna, GRAND_ANTENNA_ARM_SN, 100, &z_mag,&z_phase) !=
            EXIT_SUCCESS) {
                fprintf(stderr, "error: %s\n", grand_error_get());
                exit(EXIT_FAILURE);
        };

        printf("Z = (%.3E, %.3E) Ohm\n", z_mag*cos(z_phase), z_mag*sin(z_phase));
  double Effl_mag[3],Effl_phase[3];
  for(i=0;i<3;i++) {
    Effl_mag[i] = 0;
    Effl_phase[i] = 0;
  }
  if (antenna->effective_length(antenna, GRAND_ANTENNA_ARM_EW, 100.3, 0,60,Effl_mag,Effl_phase) !=
       EXIT_SUCCESS) {
           fprintf(stderr, "error: %s\n", grand_error_get());
           exit(EXIT_FAILURE);
   };
  printf("Effective length = (");
  for(i=0;i<3;i++) printf("%.3E+%.3EI, ", Effl_mag[i]*cos(Effl_phase[i]),Effl_mag[i]*sin(Effl_phase[i]));
  printf(")\n");

        if (antenna != NULL) antenna->destroy(&antenna);
        exit(EXIT_SUCCESS);
}
