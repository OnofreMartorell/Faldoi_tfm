#ifndef PARAMETERS_H
#define PARAMETERS_H

//METHODS
#define M_TVL1       0
#define M_TVL1_W     1
#define M_NLTVL1     2
#define M_NLTVL1_W   3
#define M_TVCSAD     4
#define M_TVCSAD_W   5
#define M_NLTVCSAD   6
#define M_NLTVCSAD_W 7
#define M_TVL1_OCC   8

//IMAGE_PARAMETERS
#define PRESMOOTHING_SIGMA  0.90


//lambda = 0.25;
//theta = 0.3;
//const float beta = 1;
//const float alpha = 0.01;
//const float tau_u = 0.125;
//const float tau_eta = 0.125;
//const float tau_chi = 0.125;


//OPTICAL FLOW PARAMETERS
#define PAR_DEFAULT_LAMBDA  1//40
#define PAR_DEFAULT_THETA   0.3
#define PAR_DEFAULT_TAU     0.125 //0.25
#define PAR_DEFAULT_BETA    1
#define PAR_DEFAULT_ALPHA   0.01
#define PAR_DEFAULT_TAU_U   0.125
#define PAR_DEFAULT_TAU_ETA 0.125
#define PAR_DEFAULT_TAU_CHI 0.1
#define PAR_DEFAULT_TOL_D   0.01
#define PAR_DEFAULT_VERBOSE 1  //0

#define PAR_DEFAULT_GAMMA 0.05  //0

#define MAX_ITERATIONS_OF 4
#define MAX_ITERATIONS 400

#define GRAD_IS_ZERO 1E-8
#define GRAD_IS_ZERO_GLOBAL 1E-10


#define PAR_DEFAULT_NPROC   0    //0

#define PAR_DEFAULT_NWARPS_LOCAL  1  //5
#define PAR_DEFAULT_NWARPS_GLOBAL  5  //5

//Specific stuff for NLTV

#define NL_SPATIAL 2
#define NL_INTENSITY 2
#define NL_BETA  2//Neighboor
#define NL_DUAL_VAR (2*NL_BETA + 1)*(2*NL_BETA + 1) - 1 // 5x5

//Specific Stuff for the CSAD
#define DT_R  3//Neighboor 7x7
#define DT_NEI (2*DT_R + 1)*(2*DT_R + 1) - 1 // 7x7




#endif // PARAMETERS_H
