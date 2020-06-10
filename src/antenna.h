#pragma unique


void grand_error_clear(void);
const char * grand_error_get(void);

enum grand_antenna_arm {
        GRAND_ANTENNA_ARM_SN = 0,
        GRAND_ANTENNA_ARM_EW,
        GRAND_ANTENNA_ARM_Z
};

struct grand_antenna {
        int (*impedance) (struct grand_antenna * antenna,
                          enum grand_antenna_arm arm,
                          double frequency,
                          double  * res_mag, double * res_phase
                         );

        int (*effective_length) (struct grand_antenna * antenna,
                                 enum grand_antenna_arm arm,
                                 double frequency,
                                 double phi,
                                 double theta,
                                 double res_mag[3], double res_phase[3]
                                );

        void (*destroy) (struct grand_antenna ** antenna);
};


int grand_antenna_load(struct grand_antenna ** antenna, const char * path);
