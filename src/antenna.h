#pragma unique

#include <complex.h>

enum grand_antenna_arm {
        GRAND_ANTENNA_ARM_SN = 0,
        GRAND_ANTENNA_ARM_EW,
        GRAND_ANTENNA_ARM_Z
};

struct grand_antenna {
        double complex (*impedance) (struct grand_antenna * antenna,
                                     enum grand_antenna_arm arm,
                                     double frequency);

        double complex (*effective_length) (struct grand_antenna * antenna,
                                            enum grand_antenna_arm arm,
                                            double frequency,
                                            double phi,
                                            double theta);

        void (*destroy) (struct grand_antenna ** antenna);
};


int grand_antenna_load(struct grand_antenna ** antenna, const char * path);
