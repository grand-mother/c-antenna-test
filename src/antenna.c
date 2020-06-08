#include <stdarg.h>
#include <stdlib.h>

#include "hdf5.h"
#include "antenna.h"


struct tabulated_antenna
{
        struct grand_antenna api;

        int n_arm;
        int n_freq;
        int n_phi;
        int n_theta;

        float * freq;
        float * phi;
        float * theta;
        float * leff_phi;
        float * leff_theta;
        float * phase_phi;
        float * phase_theta;
        float * resistance;
        float * reactance;

        float data[];
};


#define ERROR_MSG_SIZE 1024
static char last_error_msg[ERROR_MSG_SIZE] = {0x0};

static int error(char * format, ...)
{
        va_list args;
        va_start(args, format);
        vsnprintf(last_error_msg, ERROR_MSG_SIZE, format, args);
        va_end(args);

        return EXIT_FAILURE;
}
#undef ERROR_MSG_SIZE


void grand_error_clear(void)
{
        last_error_msg[0] = 0x0;
}


const char * grand_error_get(void)
{
        return (last_error_msg[0] == 0x0) ? NULL : last_error_msg;
}


static int tabulated_antenna_impedance(
    struct grand_antenna * antenna, enum grand_antenna_arm arm,
    double frequency, double complex * result)
{
        if (antenna == NULL) {
                return error("bad antenna (NULL)");
        }
        struct tabulated_antenna * t = (void *)antenna;

        if ((arm < 0) || (arm >= t->n_arm)) {
                return error("bad arm (expected a number in [0, %d], got %d)",
                             t->n_arm -1);
        }

        const int n = t->n_freq;
        const float * const f = t->freq;
        const float * const r = t->resistance + n * arm;
        const float * const x = t->reactance + n * arm;

        if ((frequency < f[0]) || (frequency > f[n - 1])) {
                return error("bad frequency (expected a value in [%g, %g], "
                             "got %g)", f[0], f[n - 1], frequency);
        }
        else if (frequency == f[n -1]) {
                *result = r[n - 1] + I * x[n - 1];
        } else {
                double h1 = (frequency - f[0]) / (f[n - 1] - f[0]) * (n - 1);
                const int i0 = (int)h1;
                h1 -= i0;
                const double h0 = 1 - h1;
                const int i1 = i0 + 1;
                *result = r[i0] * h0 + r[i1] * h1 +
                          I * (x[i0] * h0 + x[i1] * h1);
        }

        return EXIT_SUCCESS;
}


static int tabulated_antenna_effective_length(
    struct grand_antenna * antenna, enum grand_antenna_arm arm,
    double frequency, double phi, double theta, double complex result[3])
{
        /* XXX interpolate */

        return EXIT_FAILURE;
}


static void tabulated_antenna_destroy (struct grand_antenna ** antenna)
{
        if (antenna == NULL) return;
        free(*antenna);
        *antenna = NULL;
}


static int tabulated_antenna_load(struct tabulated_antenna ** antenna,
                                  const char * fname)
{
        hid_t file;
        hid_t group[3];
        hid_t dataset[3], dataspace[3], datatype[3]; //freq,phi,theta
        int rank, idata, igroup, ifreq;
        herr_t status;
        hsize_t dims[3], max_dims;
        const int n_arm = 3;
        struct tabulated_antenna * ant;

        int rc = EXIT_FAILURE;
        if (antenna == NULL) return rc;
        *antenna = NULL;

        if((file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
                goto clean_and_exit;
        if ((group[GRAND_ANTENNA_ARM_SN] =
            H5Gopen(file, "SN", H5P_DEFAULT)) < 0) goto clean_and_exit;
        if ((group[GRAND_ANTENNA_ARM_EW] =
            H5Gopen(file, "EW", H5P_DEFAULT)) < 0) goto clean_and_exit;
        if ((group[GRAND_ANTENNA_ARM_Z] =
            H5Gopen(file, "Z", H5P_DEFAULT)) < 0) goto clean_and_exit;

        //1. read frequency info (same for all, read for NS arm
        if ((dataset[0] =
            H5Dopen(group[GRAND_ANTENNA_ARM_SN], "frequency", H5P_DEFAULT)) < 0)
                goto clean_and_exit;

        //2. read phi info
        if ((dataset[1] =
            H5Dopen(group[GRAND_ANTENNA_ARM_SN], "phi", H5P_DEFAULT )) < 0)
                goto clean_and_exit;

        //3. read theta info
        if ((dataset[2] =
            H5Dopen(group[GRAND_ANTENNA_ARM_SN], "theta", H5P_DEFAULT)) < 0)
                goto clean_and_exit;

        for (idata = 0; idata < 3; idata++) {
                datatype[idata] = H5Dget_type(dataset[idata]); /* datatype identifier */
                dataspace[idata] = H5Dget_space(dataset[idata]); /* dataspace identifier */
                rank = H5Sget_simple_extent_dims(dataspace[idata], dims + idata,
                                                 &max_dims);
                if (rank != 1) goto clean_and_exit;
        }

        size_t size = sizeof(*ant) + sizeof(float) * (
            dims[0] + //list frequencies
            dims[1] + //list phi
            dims[2] + //list theta
            n_arm * (
                dims[0] * (
                    2 + //resistance+reactance
                    dims[1] * dims[2] * 4 //complex phi,theta leff
                )
            )
        );
        if ((ant = malloc(size)) == NULL)
                goto clean_and_exit; //cannot allocate free memory

        //set pointers to all starting points of the data
        ant->n_arm = n_arm;
        ant->n_freq = dims[0];
        ant->n_phi = dims[1];
        ant->n_theta = dims[2];
        ant->freq = ant->data;
        ant->phi = &(ant->freq[ant->n_freq]);
        ant->theta = &(ant->phi[ant->n_phi]);
        ant->leff_phi = &(ant->theta[ant->n_theta]);
        idata = n_arm * ant->n_freq * ant->n_phi * ant->n_theta; //4 dimensional arrays!
        ant->leff_theta = &(ant->leff_phi[idata]);
        ant->phase_phi = &(ant->leff_theta[idata]);
        ant->phase_theta = &(ant->phase_phi[idata]);
        ant->resistance = &(ant->phase_theta[idata]);
        ant->reactance =&(ant->resistance[n_arm * ant->n_freq]);
        status = H5Dread(dataset[0], datatype[0], H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         ant->freq);
        status = H5Dread(dataset[1], datatype[1], H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         ant->phi);
        status = H5Dread(dataset[2], datatype[2], H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         ant->theta);
        for (idata = 0; idata < 3; idata++)
                H5Dclose(dataset[idata]);
        idata = ant->n_freq * ant->n_phi * ant->n_theta; //amount for each arm (group in HDF5 terminology)

        /* start with something really ugly; read the resistance and reactance
         * in the leff dataspace, simply because it is not stored properly in
         * the hdf5. Afterwards, move in the right part of dataspace and
         * overwrite the leff
         */
        /* XXX the HDF5 data format needs to be changed for the new layout */
        for (igroup = 0; igroup < 3; igroup++) {
                dataset[0] = H5Dopen(group[igroup], "resistance", H5P_DEFAULT);
                if (dataset[0] < 0) goto clean_and_exit;
                datatype[0] = H5Dget_type(dataset[0]); /* datatype identifier */
                dataspace[0] = H5Dget_space(dataset[0]); /* dataspace identifier */
                status = H5Dread(dataset[0], datatype[0], H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                 &(ant->leff_phi[idata * igroup]));
                for (ifreq=0; ifreq < ant->n_freq; ifreq++) {
                        ant->resistance[igroup * ant->n_freq + ifreq] =
                            ant->leff_phi[idata * igroup + ifreq *
                                          (ant->n_phi * ant->n_theta)];
                }
                H5Dclose(dataset[0]);

                if ((dataset[0] = H5Dopen(group[igroup], "reactance",
                                          H5P_DEFAULT) ) < 0)
                        goto clean_and_exit;

                datatype[0] = H5Dget_type (dataset[0]); /* datatype identifier */
                dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
                status = H5Dread(dataset[0], datatype[0], H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT, &(ant->leff_phi[idata*igroup]));
                for(ifreq = 0; ifreq < ant->n_freq; ifreq++) {
                        ant->reactance[igroup * ant->n_freq + ifreq] =
                            ant->leff_phi[idata * igroup + ifreq *
                                          (ant->n_phi*ant->n_theta)];
                }
                H5Dclose(dataset[0]);
        }
        // that was the really ugly part; next decode the remainder of the file

        for (igroup = 0; igroup < 3; igroup++) {
                dataset[0] = H5Dopen(group[igroup], "leff_phi", H5P_DEFAULT);
                if (dataset[0] < 0) goto clean_and_exit;
                datatype[0] = H5Dget_type(dataset[0]); /* datatype identifier */
                dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
                status = H5Dread(dataset[0], datatype[0], H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT, &(ant->leff_phi[idata * igroup]));
                H5Dclose(dataset[0]);
        }

        for (igroup = 0; igroup < 3; igroup++) {
                dataset[0] = H5Dopen(group[igroup], "leff_theta", H5P_DEFAULT);
                if (dataset[0] < 0) goto clean_and_exit;
                datatype[0] = H5Dget_type(dataset[0]); /* datatype identifier */
                dataspace[0] = H5Dget_space(dataset[0]); /* dataspace identifier */
                status = H5Dread(dataset[0], datatype[0], H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT,
                                 &(ant->leff_theta[idata * igroup]));
                H5Dclose(dataset[0]);
        }

        for (igroup = 0; igroup < 3; igroup++) {
                dataset[0] = H5Dopen(group[igroup], "phase_phi", H5P_DEFAULT);
                if (dataset[0] < 0) goto clean_and_exit;
                datatype[0] = H5Dget_type(dataset[0]); /* datatype identifier */
                dataspace[0] = H5Dget_space(dataset[0]); /* dataspace identifier */
                status = H5Dread(dataset[0], datatype[0], H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT,
                                 &(ant->phase_phi[idata * igroup]));
                H5Dclose(dataset[0]);
        }

        for (igroup = 0; igroup < 3; igroup++) {
                dataset[0] = H5Dopen(group[igroup], "phase_theta", H5P_DEFAULT);
                if (dataset[0] < 0) goto clean_and_exit;
                datatype[0] = H5Dget_type(dataset[0]); /* datatype identifier */
                dataspace[0] = H5Dget_space(dataset[0]); /* dataspace identifier */
                status = H5Dread(dataset[0], datatype[0], H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT,
                                 &(ant->phase_theta[idata * igroup]));
                H5Dclose(dataset[0]);
        }


        /* Map the API and acknowledge success */
        ant->api.impedance = &tabulated_antenna_impedance;
        ant->api.effective_length = &tabulated_antenna_effective_length;
        ant->api.destroy = &tabulated_antenna_destroy;
        *antenna = ant;
        rc = EXIT_SUCCESS;

clean_and_exit:
        H5Fclose(file);
        if (*antenna == NULL) free(ant);

        return rc;
}


int grand_antenna_load(struct grand_antenna ** antenna, const char * path)
{
        return tabulated_antenna_load((struct tabulated_antenna **)antenna,
                                      path);
}
