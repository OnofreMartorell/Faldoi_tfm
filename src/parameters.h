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


//OPTICAL FLOW PARAMETERS
#define PAR_DEFAULT_LAMBDA  40 //41.8887721851 //40
#define PAR_DEFAULT_THETA   0.3 //0.683018503834 //0.3
#define PAR_DEFAULT_TAU     0.125 //0.25

#define PAR_DEFAULT_BETA    0.3213159893 //0.025
#define PAR_DEFAULT_ALPHA   0.0262151490259
#define PAR_DEFAULT_TAU_U   0.032086987536
#define PAR_DEFAULT_TAU_ETA 0.0144537526023
#define PAR_DEFAULT_TAU_CHI 0.048638480293
#define PAR_DEFAULT_MU      1.4113023753 //1

#define PAR_DEFAULT_TOL_D   0.01
#define PAR_DEFAULT_VERBOSE 0  //0
#define SAVE_RESULTS        0

//Lambda	Theta	         Beta	        Alpha	        Tau u	         Tau eta	Tau chi  	mu
//41,8887721851	0,683018503834	0,3213159893	0,0262151490259	0,032086987536	0,0144537526023	0,048638480293	1,4113023753



#define PAR_DEFAULT_GAMMA 0.05

#define MAX_ITERATIONS_LOCAL 4 //4
#define MAX_ITERATIONS_GLOBAL 400 //400

#define GRAD_IS_ZERO 1E-8
#define GRAD_IS_ZERO_GLOBAL 1E-10


#define PAR_DEFAULT_NPROC   0    //0

#define PAR_DEFAULT_NWARPS_LOCAL  1  //1
#define PAR_DEFAULT_NWARPS_GLOBAL  5  //5

#define ITER_XI 25 //12
#define ITER_CHI 25 //12
#define THRESHOLD_DELTA 0.6

#define GLOBAL_STEP 1
#define LOCAL_STEP 0

//Parameters for faldoi and prunning
#define LOCAL_ITER 3 //3
#define TU_TOL 0.01
#define FB_TOL 2

//Parameters for bilateral filter
#define PATCH_BILATERAL_FILTER 2
#define SIGMA_BILATERAL_DIST   4.0
#define SIGMA_BILATERAL_COLOR  0.08
#define ITER_BILATERAL_FILTER  10

//Specific stuff for NLTV

#define NL_SPATIAL 2
#define NL_INTENSITY 2
#define NL_BETA  2 //Neighboor
#define NL_DUAL_VAR (2*NL_BETA + 1)*(2*NL_BETA + 1) - 1 // 5x5

//Specific Stuff for the CSAD
#define DT_R  3 //Neighboor 7x7
#define DT_NEI (2*DT_R + 1)*(2*DT_R + 1) - 1 // 7x7

#define MAX_PATCH 50


#endif // PARAMETERS_H
