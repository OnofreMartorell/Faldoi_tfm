#ifndef TVL2_MODEL_H
#define TVL2_MODEL_H

#include "energy_structures.h"
#include "aux_energy_model.h"


//OPTICAL FLOW PARAMETERS
#define TVL2_LAMBDA  40//40
#define TVL2_THETA   0.3
#define TVL2_TAU     0.125 //0.25
#define TVL2_NWARPS  1  //5
#define TVL2_TOL_D   0.01
#define TVL2_VERBOSE 0  //0


////INITIALIZATION OF EACH METHOD
void  intialize_stuff_tvl2coupled(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore
          );
void  free_stuff_tvl2coupled(SpecificOFStuff *ofStuff);

void eval_tvl2coupled(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    Tvl2CoupledOFStuff *tvl2,
    float *ener_N,
    const int ii, // initial column
    const int ij, // initial row
    const int ei, // end column
    const int ej, // end row
    const float lambda,  // weight of the data term
    const float theta
    );

// Variational Optical flow method based on initial fixed values
// It minimizes the energy of \int_{B(x)} ||J(u)|| + |I_{1}(x+u)-I_{0}(x)|
// s.t u = u_0 for i.seeds
// J(u) = (u_x, u_y; v_x, v_y)
void guided_tvl2coupled(
    const float *I0,           // source image
    const float *I1,           // target image
    OpticalFlowData *ofD,
    Tvl2CoupledOFStuff *tvl2,
    float *ener_N,
    const int ii, // initial column
    const int ij, // initial row
    const int ei, // end column
    const int ej, // end row
    const float lambda,  // weight of the data term
    const float theta,   // weight of the data term
    const float tau,     // time step
    const float tol_OF,  // tol max allowed
    const int   warps,   // number of warpings per scale
    const bool  verbose  // enable/disable the verbose mode
    );

#endif //TVL2-L1 functional
