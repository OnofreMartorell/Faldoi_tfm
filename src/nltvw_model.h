#ifndef NLTVL1_MODEL_W_H
#define NLTVL1_MODEL_W_H



void  intialize_stuff_nltvl1_w(
          SpecificOFStuff *ofStuff,
          OpticalFlowData *ofCore);



void  free_stuff_nltvl1_w(SpecificOFStuff *ofStuff);


void eval_nltvl1_w(
    float *I0,           // source image
    float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTVL1Stuff_W *nltvl1w,
    float *ener_N,
    const int ii, // initial column
    const int ij, // initial row
    const int ei, // end column
    const int ej, // end row
    const float lambda,  // weight of the data term
    const float theta
    );


void guided_nltvl1_w(
    float *I0,           // source image
    float *I1,           // target image
    OpticalFlowData *ofD,
    NonLocalTVL1Stuff_W *nltvl1w,
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

#endif //NLTVL1_MODEL_W