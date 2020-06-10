#include <stdarg.h>
#include <stdlib.h>
#include<math.h>

#include "hdf5.h"
#include "antenna.h"
#define RADDEG 57.295779513082321

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
    double frequency, double * res_mag, double * res_phase)
{
        if (antenna == NULL) {
                return error("bad antenna (NULL)");
        }
        struct tabulated_antenna * t = (struct tabulated_antenna *)antenna;

        if ((arm < 0) || (arm >= t->n_arm)) {
                return error("bad arm (expected a value in [0, %d], got %d)",
                             t->n_arm -1);
        }

        const int n = t->n_freq;
        const float * const f = t->freq;
        const float * const r = t->resistance + n * arm;
        const float * const x = t->reactance + n * arm;

        double re, im;
        if ((frequency < f[0]) || (frequency > f[n - 1])) {
                return error("bad frequency (expected a value in [%g, %g], "
                             "got %g)", f[0], f[n - 1], frequency);
        } else if (frequency == f[n -1]) {
                re = r[n - 1];
                im = x[n - 1];
        } else {
                double h1 = (frequency - f[0]) / (f[n - 1] - f[0]) * (n - 1);
                const int i0 = (int)h1;
                h1 -= i0;
                const double h0 = 1 - h1;
                const int i1 = i0 + 1;
                re = r[i0] * h0 + r[i1] * h1;
                im = x[i0] * h0 + x[i1] * h1;
        }

        if (res_mag != NULL) *res_mag = sqrt(re * re + im * im);
        if (res_phase != NULL) {
                if ((re != 0) || (im != 0))
                        *res_phase = atan2(im, re);
                else
                        *res_phase = 0;
        }

        return EXIT_SUCCESS;
}


static int tabulated_antenna_effective_length(
    struct grand_antenna * antenna, enum grand_antenna_arm arm,
    double frequency, double phi, double theta,
    double res_mag[3],double res_phase[3])
{
        double weight[8]; // (0-3: low freq, even: low theta, 0,1,4,5: low phi)
        int index[8]; // to make it easier...
        double Lsph_real[3], Lsph_imag[3];
        double Lcart_real[3], Lcart_imag[3];
        // phi and theta are in degrees!
        double phi_rad = phi/RADDEG;
        double theta_rad = theta/RADDEG;

        if (antenna == NULL)
                return error("bad antenna (NULL)");
        struct tabulated_antenna * t = (struct tabulated_antenna *)antenna;

        if ((arm < 0) || (arm >= t->n_arm)) {
                return error("bad arm (expected a number in [0, %d], got %d)",
                             t->n_arm -1);
        }

        // freq. dependent weights; easy interpolation
        const int nf = t->n_freq;
        const float * const f = t->freq;
        if ((frequency < f[0]) || (frequency > f[nf - 1]))
          return error("bad frequency (expected a value in [%g, %g], got %g)",
                       f[0], f[nf - 1], frequency);
        double h1f = (frequency - f[0]) / (f[nf - 1] - f[0]) * (nf - 1);
        int i0f = (int)h1f;
        int i1f = i0f+1;
        if(i1f == nf) i1f = i0f;

        // Theta dependent weights; easy interpolation
        const int nt = t->n_theta;
        const float * const th = t->theta;
        if ((theta < th[0]) || (theta > th[nt - 1])) {
                return error("bad theta (expected a value in [%g, %g], got %g)",
                             th[0], th[nt - 1], theta);
        }
        double h1t = (theta- th[0]) / (th[nt - 1] - th[0]) * (nt - 1);
        int i0t = (int)h1t;
        int i1t = i0t+1;
        if(i1t == nt) i1t=i0t;

        // Phi dependent weights; easy interpolation Note that phi 0==360!
        const int np = t->n_phi;
        const float * const p = t->phi;
        double h1p = (fmod(phi,360) - p[0]) / (p[np - 1] - p[0]) * (np - 1);
        int i0p = (int)h1p;
        int i1p = i0p+1;
        if(i1p == np) i1p = i0p;

        weight[0] = (i1f - h1f) * (i1t - h1t) * (i1p - h1p);
        index[0] = (i0f * np + i0p) * nt + i0t;
        weight[1] = (i1f - h1f) * (h1t - i0t) * (i1p - h1p);
        index[1] = (i0f * np + i0p) * nt + i1t;
        weight[2] = (i1f - h1f) * (i1t - h1t) * (h1p - i0p);
        index[2] = (i0f * np + i1p) * nt + i0t;
        weight[3] = (i1f - h1f) * (h1t - i0t) * (h1p - i0p);
        index[3] = (i0f * np + i1p) * nt + i1t;
        weight[4] = (h1f - i0f) * (i1t - h1t) * (i1p - h1p);
        index[4] = (i1f * np + i0p) * nt + i0t;
        weight[5] = (h1f - i0f) * (h1t - i0t) * (i1p - h1p);
        index[5] = (i1f * np + i0p) * nt + i1t;
        weight[6] = (h1f - i0f) * (i1t - h1t) * (h1p - i0p);
        index[6] = (i1f * np + i1p) * nt + i0t;
        weight[7] = (h1f - i0f) * (h1t - i0t) * (h1p - i0p);
        index[7] = (i1f * np + i1p) * nt + i1t;
        Lsph_real[0] = 0; //not used anyway (L_R)
        Lsph_imag[0] = 0; //not used anyway (L_R)
        const float * const lefp = t->leff_phi + nf * nt * np * arm;
        const float * const lpp = t->phase_phi + nf * nt * np * arm;
        Lsph_real[1] = 0;
        Lsph_imag[1] = 0;
        for(int i = 0; i < 8; i++) {
                Lsph_real[1] += weight[i] * lefp[index[i]] * cos(lpp[index[i]]);
                Lsph_imag[1] += weight[i] * lefp[index[i]] * sin(lpp[index[i]]);
        }
        Lsph_real[2] = 0;
        Lsph_imag[2] = 0;
        const float * const left = t->leff_theta + nf * nt * np * arm;
        const float * const lpt = t->phase_theta + nf * nt * np * arm;
        for(int i = 0; i < 8; i++) {
                Lsph_real[2] += weight[i] * left[index[i]] * cos(lpt[index[i]]);
                Lsph_imag[2] += weight[i] * left[index[i]] * sin(lpt[index[i]]);
        }
        Lcart_real[0] = -Lsph_real[1] * sin(phi_rad) +
                         Lsph_real[2] * cos(phi_rad) * cos(theta_rad);
        Lcart_real[1] = Lsph_real[1] * cos(phi_rad) +
                        Lsph_real[2] * sin(phi_rad) * cos(theta_rad);
        Lcart_real[2] = -Lsph_real[2] * sin(theta_rad);
        Lcart_imag[0] = -Lsph_imag[1] * sin(phi_rad) +
                         Lsph_imag[2] * cos(phi_rad) * cos(theta_rad);
        Lcart_imag[1] = Lsph_imag[1] * cos(phi_rad) +
                        Lsph_imag[2] * sin(phi_rad) * cos(theta_rad);
        Lcart_imag[2] = -Lsph_imag[2] * sin(theta_rad);

        if (res_mag != NULL) {
                for(int i = 0; i < 3; i++)
                        res_mag[i] = sqrt(Lcart_real[i] * Lcart_real[i] +
                                          Lcart_imag[i] * Lcart_imag[i]);
        }

        if (res_phase != NULL) {
                for(int i = 0; i < 3; i++) {
                        if ((Lcart_real[i] != 0) || (Lcart_imag[i] != 0)) {
                                res_phase[i] = atan2(Lcart_imag[i],
                                                     Lcart_real[i]);
                        } else {
                                res_phase[i] = 0;
                        }
                }
        }

        return EXIT_SUCCESS;
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
        size_t size;

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

        size = sizeof(*ant) + sizeof(float) * (
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
        if ((ant = (struct tabulated_antenna *)malloc(size)) == NULL)
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
