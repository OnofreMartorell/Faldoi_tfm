#define MAX_ITERATIONS_OF 4
#define GRAD_IS_ZERO 1E-8
/////////////////////////////


//Compute the divergence (backward differences) over a patch
void divergence_patch(
    const float *v1, // x component of the vector field
    const float *v2, // y component of the vector field
    float *div,      // output divergence
    const int ii,     // initial column
    const int ij,     // initial row
    const int ei,     // end column
    const int ej,     // end row
    const int nx    // image width
    );


//TODO: Poner explicacion
//Compute the forward gradient (forward difference) over a patch
void forward_gradient_mixed_bound(
    const float *f, //input image
    float *fx,      //computed x derivative
    float *fy,      //computed y derivative
    const int ii,     // initial column
    const int ij,     // initial row
    const int ei,     // end column
    const int ej,     // end row
    const int nx,   //image width
    const int ny   //image height
    );


//Compute the forward gradient (forward difference) over a patch
void forward_gradient_patch(
    const float *f, //input image
    float *fx,      //computed x derivative
    float *fy,      //computed y derivative
    const int ii,     // initial column
    const int ij,     // initial row
    const int ei,     // end column
    const int ej,     // end row
    const int nx   //image width
    );

//Create 
void bicubic_interpolation_warp_patch(
  const float *input,     // image to be warped
  const float *u,         // x component of the vector field
  const float *v,         // y component of the vector field
  float       *output,    // image warped with bicubic interpolation
  const int    ii,     // initial column
  const int    ij,     // initial row
  const int    ei,     // end column
  const int    ej,     // end row
  const int    nx,        // image width
  const int    ny,        // image height
  bool         border_out // if true, put zeros outside the region
  );
////////////////////////////////////////////////////////////////////////////////
////////////////////////AUXILIAR FUNCTIONS NLTVL1///////////////////////////////
 bool positive(int val);
 float getsample_0(float *x, int w, int h, int pd, int i, int j, int l);

 int validate_ap(int w, int h, int i, int j, int di, int dj);

 int validate_ap_patch(
                const int ii,     // initial column
                const int ij,     // initial row
                const int ei,     // end column
                const int ej,     // end row
                const int w,
                int i, 
                int j
  );

 float get_wspatial(int l, int k);
 float get_wcolor(float *a, int w, int h, int i, int j, int l, int k,int pd);

 float get_weight(float *a, int w, int h, int i, int j, int l, int k,int pd);

 void nltv_ini_dual_variables(
                float *a,
                const int pd,
                const int w,
                const int h,
                const int n_d,
                const int radius,
                DualVariables *p,
                DualVariables *q
  );

 void non_local_divergence(
            DualVariables *p,
            const int ii, // initial column
            const int ij, // initial row
            const int ei, // end column
            const int ej, // end row
            const int w,
            int n_d,
            float *div_p
    );

 //////////////////////////CSAD///////////////////////////////////////////////
 void csad_ini_pos_nei(
                const int w,
                const int h,
                const int ndt,
                const int rdt,
                PosNei *pnei
  );
