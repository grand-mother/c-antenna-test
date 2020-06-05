// number of bins for all parameters of the antenna model
// Note that phi and theta are wrt GRAND conventions!
// phi wrt North, theta = Zenith angle (all in degrees)
// the phases are in radians though!
#define ANT_N_ARM 0
#define ANT_N_FREQ 1
#define ANT_N_PHI 2
#define ANT_N_THETA 3

#define ANT_NS 0 //North-South arm
#define ANT_EW 1 //East-West Arm
#define ANT_V 2 //Vertical Arm
#define ANT_X 0 //X arm
#define ANT_Y 1 //Y Arm
#define ANT_Z 2 //Z Arm

typedef struct{
  hsize_t dim[4]; //Orientation,frequency,phi,theta
  float *freq,*phi,*theta;
  float *leff_phi;//arm,freq,phi,theta
  float *leff_theta; //arm,freq,phi,theta
  float *phase_phi;//arm,freq,phi,theta
  float *phase_theta;//arm,freq,phi,theta
  float *resistance; //arm,freq
  float *reactance; //arm,freq
}Antenna;
