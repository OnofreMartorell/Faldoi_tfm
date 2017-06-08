#ifndef MATCH_VARIATIONAL_CORE_DATA_H
#define MATCH_VARIATIONAL_CORE_DATA_H

#include <vector>

#define MAX(x,y) ((x)>(y)?(x):(y))

#define M_TVL1       0
#define M_TVL1_W     1
#define M_NLTVL1     2 
#define M_NLTVL1_W   3 
#define M_TVCSAD     4
#define M_TVCSAD_W   5
#define M_NLTVCSAD   6
#define M_NLTVCSAD_W 7
#define M_TVL1_OCC   8

//Specific stuff for NLTV

#define NL_SPATIAL 2
#define NL_INTENSITY 2
#define NL_BETA  2//Neighboor
#define NL_DUAL_VAR (2*NL_BETA + 1)*(2*NL_BETA + 1) -1 // 5x5

//Specific Stuff for the CSAD
#define DT_R  3//Neighboor 7x7
#define DT_NEI (2*DT_R + 1)*(2*DT_R + 1) -1 // 7x7


struct BilateralWeight{
    float wp[NL_DUAL_VAR]; // weight of non local
    int   apj[NL_DUAL_VAR]; //absolute positon of p(y,x) (row)
    int   api[NL_DUAL_VAR]; //absolute position of p(y,x) (colum)
    float wt = 0.0;
};

struct OpticalFlowData{
    /* data */
    //TODO: This should be outside of this structure
    float *u1;
    float *u2;
    float *u1_ba;
    float *u2_ba;
    float *chi;
    int   *fixed_points;
    int   *trust_points;
    float *saliency; //It stores the saliency value for each pixel.
    BilateralWeight *weight;
    float *u1_ini; //Optical flow values to initialize the local patch
    float *u2_ini; //Optical flow values to initialize the local patch
    int wr;
    int w;
    int h;
    int method; //
};


struct DualVariables{
    float sc[NL_DUAL_VAR]; // value of p(x,y)
    float wp[NL_DUAL_VAR]; // weight of non local
    int   apj[NL_DUAL_VAR]; //absolute positon of p(y,x) (row)
    int   api[NL_DUAL_VAR]; //absolute position of p(y,x) (colum)
    int   rp[NL_DUAL_VAR]; //relative position of p(y,x) in the structure
    float wt = 0.0;
};


//Struct
struct PosNei{
    int   api[DT_NEI]; //absolute positon of Intensity (row)
    int   apj[DT_NEI]; //absolute position of intensity (colum)
    float b[DT_NEI];
    std::vector<float>  ba;
    int n;
};

////Specific struct for the different functionals
struct  Tvl2CoupledOFStuff{
    //Dual variables
    float *xi11 = nullptr;
    float *xi12 = nullptr;
    float *xi21 = nullptr;
    float *xi22 = nullptr;

    float *u1x = nullptr;
    float *u1y = nullptr;
    float *u2x = nullptr;
    float *u2y = nullptr;

    float *v1 = nullptr;
    float *v2 = nullptr;

    float *rho_c = nullptr;
    float *grad = nullptr;

    float *u1_ = nullptr;
    float *u2_ = nullptr;

    float *u1Aux = nullptr;
    float *u2Aux = nullptr;

    float *I1x = nullptr;
    float *I1y = nullptr;
    float *I1w = nullptr;
    float *I1wx = nullptr;
    float *I1wy = nullptr;

    float *div_xi1 = nullptr;
    float *div_xi2 = nullptr;
    float *u_N = nullptr;
};

// ////Specific struct for the different functionals
// struct  Tvl1OFStuff
// {
//   float *xi11;
//   float *xi12;
//   float *xi21;
//   float *xi22;
//   float *u1x;
//   float *u1y;
//   float *u2x;
//   float *u2y;
//   float *v1;
//   float *v2;
//   float *rho_c;
//   float *grad;
//   float *u1_;
//   float *u2_;
//   float *u1Aux;
//   float *u2Aux;
//   float *I1x;
//   float *I1y;
//   float *I1w;
//   float *I1wx;
//   float *I1wy;
//   float *div_xi1;
//   float *div_xi2;
//   float *u_N;
// };

// struct RotInvariantStuff
// {
//   float *xi11;
//   float *xi12;
//   float *xi22;
//   float *u1x;
//   float *u1y;
//   float *u2x;
//   float *u2y;
//   float *v1;
//   float *v2;
//   float *rho_c;
//   float *grad;
//   float *u1_;
//   float *u2_;
//   float *u1Aux;
//   float *u2Aux;
//   float *I1x;
//   float *I1y;
//   float *I1w;
//   float *I1wx;
//   float *I1wy;
//   float *div_xi1;
//   float *div_xi2;
//   float *u_N;
// };



struct NonLocalTVL1Stuff
{
    DualVariables *p;
    DualVariables *q;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};


struct  TvCsadStuff
{
    PosNei *pnei;
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;
    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_xi1;
    float *div_xi2;
};

struct NonLocalTvCsadStuff
{
    DualVariables *p;
    DualVariables *q;
    PosNei *pnei;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};

///////////////////PESOS///////////////////////////////////////////////////////
struct  Tvl2CoupledOFStuff_W
{
    int iiw;
    int ijw;
    float *weight;
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;
    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1Aux;
    float *u2Aux;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_xi1;
    float *div_xi2;
    float *u_N;
};

struct NonLocalTvCsadStuff_W
{
    int iiw;
    int ijw;
    float *weight;
    DualVariables *p;
    DualVariables *q;
    PosNei *pnei;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};


struct NonLocalTVL1Stuff_W
{
    int iiw;
    int ijw;
    float *weight;
    DualVariables *p;
    DualVariables *q;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_p;
    float *div_q;
    float *u_N;
};

struct  TvCsadStuff_W
{
    int iiw;
    int ijw;
    float *weight;
    PosNei *pnei;
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;
    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;
    float *v1;
    float *v2;
    float *rho_c;
    float *grad;
    float *u1_;
    float *u2_;
    float *u1_tmp;
    float *u2_tmp;
    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;
    float *div_xi1;
    float *div_xi2;
};
////////////////////////
/////////OCCLUSIONS/////
///////////////////////
struct  Tvl2CoupledOFStuff_occ
{
    //Occlusion variable
    float *chi;
    float *chix;
    float *chiy;

    //Weigth
    float *g;

    float *diff_u_N;

    //Dual variables
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;

    float *u1x;
    float *u1y;
    float *u2x;
    float *u2y;
    //Dual variables for chi
    float *eta1;
    float *eta2;


    float *v1;
    float *v2;

    float *rho_c1;
    float *rho_c_1;
    float *grad_1;
    float *grad__1;


    float *I1x;
    float *I1y;
    float *I1w;
    float *I1wx;
    float *I1wy;

    float *I_1x;
    float *I_1y;
    float *I_1w;
    float *I_1wx;
    float *I_1wy;


    float *vi_div1;
    float *grad_x1;
    float *grad_y1;
    float *vi_div2;
    float *grad_x2;
    float *grad_y2;
    float *g_xi11;
    float *g_xi12;
    float *g_xi21;
    float *g_xi22;
    float *div_g_xi1;
    float *div_g_xi2;


    float *F;
    float *G;

    float *div_u;
    float *g_eta1;
    float *g_eta2;
    float *div_g_eta;
};



/////////////////////////////////////


//General Struct for the auxiliar stuff.
//Each pointer contains the auxiliar necessary information to estimate the of. 
struct  SpecificOFStuff
{ //TODO: Should think a best option. To link the things
    //Creo que el problema viene por como declaramos los punteros y todo eso.

    Tvl2CoupledOFStuff  tvl2;
    NonLocalTVL1Stuff   nltvl1;
    TvCsadStuff         tvcsad;
    NonLocalTvCsadStuff nltvcsad;

    Tvl2CoupledOFStuff_W  tvl2w;
    NonLocalTvCsadStuff_W nltvcsadw;
    NonLocalTVL1Stuff_W   nltvl1w;
    TvCsadStuff_W         tvcsadw;

    Tvl2CoupledOFStuff_occ  tvl2_occ;

};


#endif// ENERGY_STRUCTURES_H
