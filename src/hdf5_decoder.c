#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"hdf5.h"
#include "GRAND_antenna.h"



Antenna HorizonAntenna;

int read_antenna(char *fname, Antenna *ant)
{
  hid_t file;
  hid_t group[3];
  hid_t dataset[3],dataspace[3],datatype[3]; //freq,phi,theta
  int rank,idata,igroup,ifreq;
  herr_t status;
  hsize_t dims,max_dims;
  
  if(ant->freq != NULL){
    free(ant->freq);
    memset(ant,0,sizeof(Antenna));
  }
  if((file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT))<0) return(-1);
  if((group[ANT_NS] = H5Gopen(file, "SN", H5P_DEFAULT))<0) return(-1);
  if((group[ANT_EW] = H5Gopen(file, "EW", H5P_DEFAULT))<0) return(-1);
  if((group[ANT_Z] = H5Gopen(file, "Z", H5P_DEFAULT))<0) return(-1);
  ant->dim[ANT_N_ARM] = 3; // 3 arms
  //1. read frequency info (same for all, read for NS arm
  if((dataset[0] = H5Dopen(group[ANT_NS], "frequency", H5P_DEFAULT))<0) {
    memset(ant,0,sizeof(Antenna));
    return(-1);
  }
  //2. read phi info
  if((dataset[1] = H5Dopen(group[ANT_NS], "phi", H5P_DEFAULT))<0) {
    memset(ant,0,sizeof(Antenna));
    return(-1);
  }
  //3. read theta info
  if((dataset[2] = H5Dopen(group[ANT_NS], "theta", H5P_DEFAULT))<0) {
    memset(ant,0,sizeof(Antenna));
    return(-1);
  }
  for(idata=0;idata<3;idata++){
    datatype[idata] = H5Dget_type (dataset[idata]); /* datatype identifier */
    dataspace[idata] = H5Dget_space (dataset[idata]); /* dataspace identifier */
    rank = H5Sget_simple_extent_dims (dataspace[idata], &(ant->dim[idata+1]), &max_dims);
    if(rank!= 1 ) {
      memset(ant,0,sizeof(Antenna));
      return(-2); //file does not match description
    }
  }
  ant->freq = malloc((ant->dim[ANT_N_FREQ]+ //list frequencies
                      ant->dim[ANT_N_PHI]+   //list phi
                      ant->dim[ANT_N_THETA]+ //list theta
                      ant->dim[ANT_N_ARM]*(
                      ant->dim[ANT_N_FREQ]*(2+ //resistance+reactance
                           ant->dim[ANT_N_PHI]*ant->dim[ANT_N_THETA]*4 //complex phi,theta leff
                                            ))
                      )*sizeof(float));
  if(ant->freq == NULL){
    memset(ant,0,sizeof(Antenna));
    return(-3); //cannot allocate free memory
  }
  //set pointers to all starting points of the data
  ant->phi = &(ant->freq[ant->dim[ANT_N_FREQ]]);
  ant->theta = &(ant->phi[ant->dim[ANT_N_PHI]]);
  ant->leff_phi = &(ant->theta[ant->dim[ANT_N_THETA]]);
  idata =ant->dim[ANT_N_ARM]*ant->dim[ANT_N_FREQ]*
    ant->dim[ANT_N_PHI]*ant->dim[ANT_N_THETA]; //4 dimensional arrays!
  ant->leff_theta = &(ant->leff_phi[idata]);
  ant->phase_phi = &(ant->leff_theta[idata]);
  ant->phase_theta = &(ant->phase_phi[idata]);
  ant->resistance = &(ant->phase_theta[idata]);
  ant->reactance =&(ant->resistance[ant->dim[ANT_N_ARM]*ant->dim[ANT_N_FREQ]]);
  status = H5Dread(dataset[0],datatype[0],H5S_ALL,H5S_ALL,H5P_DEFAULT,ant->freq);
  status = H5Dread(dataset[1],datatype[1],H5S_ALL,H5S_ALL,H5P_DEFAULT,ant->phi);
  status = H5Dread(dataset[2],datatype[2],H5S_ALL,H5S_ALL,H5P_DEFAULT,ant->theta);
  for(idata=0;idata<3;idata++) H5Dclose(dataset[idata]);
  idata=ant->dim[ANT_N_FREQ]*ant->dim[ANT_N_PHI]*ant->dim[ANT_N_THETA]; //amount for each arm (group in HDF5 terminology)
  //start with something really ugly; read the resistance and reactance in the leff dataspace, simply because it is not stored properly in the hdf5. Afterwards, move in the right part of dataspace and overwrite the leff
  for(igroup=0;igroup<3;igroup++){
    if((dataset[0] = H5Dopen(group[igroup], "resistance", H5P_DEFAULT))<0) {
      free(ant->freq);
      memset(ant,0,sizeof(Antenna));
      return(-1);
    }
    datatype[0] = H5Dget_type (dataset[0]); /* datatype identifier */
    dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
    status = H5Dread(dataset[0],datatype[0],H5S_ALL,H5S_ALL,H5P_DEFAULT,
                     &(ant->leff_phi[idata*igroup]));
    for(ifreq=0;ifreq<ant->dim[ANT_N_FREQ];ifreq++){
      ant->resistance[igroup*ant->dim[ANT_N_FREQ]+ifreq] =ant->leff_phi[idata*igroup+ifreq*(ant->dim[ANT_N_PHI]*ant->dim[ANT_N_THETA])];
    }
    H5Dclose(dataset[0]);
    if((dataset[0] = H5Dopen(group[igroup], "reactance", H5P_DEFAULT))<0) {
      free(ant->freq);
      memset(ant,0,sizeof(Antenna));
      return(-1);
    }
    datatype[0] = H5Dget_type (dataset[0]); /* datatype identifier */
    dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
    status = H5Dread(dataset[0],datatype[0],H5S_ALL,H5S_ALL,H5P_DEFAULT,
                     &(ant->leff_phi[idata*igroup]));
    for(ifreq=0;ifreq<ant->dim[ANT_N_FREQ];ifreq++){
      ant->reactance[igroup*ant->dim[ANT_N_FREQ]+ifreq] =ant->leff_phi[idata*igroup+ifreq*(ant->dim[ANT_N_PHI]*ant->dim[ANT_N_THETA])];
    }
    H5Dclose(dataset[0]);
  }
// that was the really ugly part; next decode the remainder of the file
  for(igroup=0;igroup<3;igroup++){
    if((dataset[0] = H5Dopen(group[igroup], "leff_phi", H5P_DEFAULT))<0) {
      free(ant->freq);
      memset(ant,0,sizeof(Antenna));
      return(-1);
    }
    datatype[0] = H5Dget_type (dataset[0]); /* datatype identifier */
    dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
    status = H5Dread(dataset[0],datatype[0],H5S_ALL,H5S_ALL,H5P_DEFAULT,
                     &(ant->leff_phi[idata*igroup]));
    H5Dclose(dataset[0]);
  }
  for(igroup=0;igroup<3;igroup++){
    if((dataset[0] = H5Dopen(group[igroup], "leff_theta", H5P_DEFAULT))<0){
      free(ant->freq);
      memset(ant,0,sizeof(Antenna));
      return(-1);
    } 
    datatype[0] = H5Dget_type (dataset[0]); /* datatype identifier */
    dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
    status = H5Dread(dataset[0],datatype[0],H5S_ALL,H5S_ALL,H5P_DEFAULT,
                     &(ant->leff_theta[idata*igroup]));
    H5Dclose(dataset[0]);
  }
  for(igroup=0;igroup<3;igroup++){
    if((dataset[0] = H5Dopen(group[igroup], "phase_phi", H5P_DEFAULT))<0){
      free(ant->freq);
      memset(ant,0,sizeof(Antenna));
      return(-1);
    }
    datatype[0] = H5Dget_type (dataset[0]); /* datatype identifier */
    dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
    status = H5Dread(dataset[0],datatype[0],H5S_ALL,H5S_ALL,H5P_DEFAULT,
                     &(ant->phase_phi[idata*igroup]));
    H5Dclose(dataset[0]);
  }
  for(igroup=0;igroup<3;igroup++){
    if((dataset[0] = H5Dopen(group[igroup], "phase_theta", H5P_DEFAULT))<0){
      free(ant->freq);
      memset(ant,0,sizeof(Antenna));
      return(-1);
    }
    datatype[0] = H5Dget_type (dataset[0]); /* datatype identifier */
    dataspace[0] = H5Dget_space (dataset[0]); /* dataspace identifier */
    status = H5Dread(dataset[0],datatype[0],H5S_ALL,H5S_ALL,H5P_DEFAULT,
                     &(ant->phase_theta[idata*igroup]));
    H5Dclose(dataset[0]);
  }
  for(igroup=0;igroup<3;igroup++)H5Gclose(group[igroup]);
  H5Fclose(file);
  return(0);
}

int main()
{
  int ifreq,iphi,ith;
  //initialization
  memset(&HorizonAntenna,0,sizeof(Antenna));
  read_antenna("HorizonAntenna_leff_loaded.hdf5",&HorizonAntenna);
  printf("Dimensions: %llu %llu %llu %llu\n",
         HorizonAntenna.dim[0],HorizonAntenna.dim[1],
         HorizonAntenna.dim[2],HorizonAntenna.dim[3]);
  for(ifreq=0;ifreq<HorizonAntenna.dim[ANT_N_FREQ];ifreq++) printf("\t%g",HorizonAntenna.freq[ifreq]);
  printf("\n");
  for(ifreq=0;ifreq<HorizonAntenna.dim[ANT_N_FREQ];ifreq++) printf("\t%g",HorizonAntenna.resistance[HorizonAntenna.dim[ANT_N_FREQ]+ifreq]);
  printf("\n");
  for(iphi=0;iphi<HorizonAntenna.dim[ANT_N_PHI];iphi++) printf("\t%g",HorizonAntenna.phi[iphi]);
  printf("\n");
  for(ith=0;ith<HorizonAntenna.dim[ANT_N_THETA];ith++) printf("\t%g",HorizonAntenna.theta[ith]);
  printf("\n");
  /*
  ith = 45;
  printf("%g MHz and %g Zenith\n",HorizonAntenna.arm[ANT_EW].freqdata[80],HorizonAntenna.arm[ANT_EW].thetadata[ith]);
  for(iphi=0;iphi<ANT_N_PHI;iphi++){
    if((iphi%5)==0)printf("\n%g",HorizonAntenna.arm[ANT_EW].phidata[iphi]);
    printf("\t%g",HorizonAntenna.arm[ANT_EW].leff_phi[80][iphi][ith]);
  }
  printf("\n");
  for(iphi=0;iphi<ANT_N_PHI;iphi++){
    if((iphi%5)==0)printf("\n%g",HorizonAntenna.arm[ANT_NS].phidata[iphi]);
    printf("\t%g",HorizonAntenna.arm[ANT_NS].leff_phi[80][iphi][ith]);
  }
  printf("\n");
  for(iphi=0;iphi<ANT_N_PHI;iphi++){
    if((iphi%5)==0)printf("\n%g",HorizonAntenna.arm[ANT_V].phidata[iphi]);
    printf("\t%g",HorizonAntenna.arm[ANT_V].leff_phi[80][iphi][ith]);
  }
  printf("\n");
*/
}
