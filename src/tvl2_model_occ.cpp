#ifndef TVL2_MODEL_OCC
#define TVL2_MODEL_OCC

#include <cmath>
#include <cstdio>
#include <cassert>
#include "energy_structures.h"
#include "aux_energy_model.h"


////INITIALIZATION OF EACH METHOD
void  intialize_stuff_tvl2coupled_occ(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore) {

    const int w = ofCore->w;
    const int h = ofCore->h;

    //Occlusion variable
    ofStuff->tvl2_occ.chi =  new float[w*h];
    ofStuff->tvl2_occ.chix =  new float[w*h];
    ofStuff->tvl2_occ.chiy =  new float[w*h];

    //Weigth
    ofStuff->tvl2_occ.g =  new float[w*h];

    //Dual variables
    ofStuff->tvl2_occ.xi11 =  new float[w*h];
    ofStuff->tvl2_occ.xi12 =  new float[w*h];
    ofStuff->tvl2_occ.xi21 =  new float[w*h];
    ofStuff->tvl2_occ.xi22 =  new float[w*h];

    //Derivatives of u1,u2
    ofStuff->tvl2_occ.u1x =  new float[w*h];
    ofStuff->tvl2_occ.u1y =  new float[w*h];
    ofStuff->tvl2_occ.u2x =  new float[w*h];
    ofStuff->tvl2_occ.u2y =  new float[w*h];

    ofStuff->tvl2_occ.v1 =  new float[w*h];
    ofStuff->tvl2_occ.v2 =  new float[w*h];

    ofStuff->tvl2_occ.rho_c1 =  new float[w*h];
    ofStuff->tvl2_occ.rho_c2 =  new float[w*h];
    ofStuff->tvl2_occ.grad =  new float[w*h];
    ofStuff->tvl2_occ.grad_ba =  new float[w*h];

    ofStuff->tvl2_occ.u1Aux =  new float[w*h];
    ofStuff->tvl2_occ.u2Aux =  new float[w*h];

    ofStuff->tvl2_occ.I2x =  new float[w*h];
    ofStuff->tvl2_occ.I2y =  new float[w*h];
    ofStuff->tvl2_occ.I2w =  new float[w*h];
    ofStuff->tvl2_occ.I2wx =  new float[w*h];
    ofStuff->tvl2_occ.I2wy =  new float[w*h];

    ofStuff->tvl2_occ.I0x =  new float[w*h];
    ofStuff->tvl2_occ.I0y =  new float[w*h];
    ofStuff->tvl2_occ.I0w =  new float[w*h];
    ofStuff->tvl2_occ.I0wx =  new float[w*h];
    ofStuff->tvl2_occ.I0wy =  new float[w*h];

    ofStuff->tvl2_occ.div_xi1 =  new float[w*h];
    ofStuff->tvl2_occ.div_xi2 =  new float[w*h];
    ofStuff->tvl2_occ.u_N =  new float[w*h];
}

void  free_stuff_tvl2coupled_occ(SpecificOFStuff *ofStuff){

    delete [] ofStuff->tvl2_occ.chi;
    delete [] ofStuff->tvl2_occ.chix;
    delete [] ofStuff->tvl2_occ.chiy;
    delete [] ofStuff->tvl2_occ.g;
    delete [] ofStuff->tvl2_occ.xi11;
    delete [] ofStuff->tvl2_occ.xi12;
    delete [] ofStuff->tvl2_occ.xi21;
    delete [] ofStuff->tvl2_occ.xi22;
    delete [] ofStuff->tvl2_occ.u1x;
    delete [] ofStuff->tvl2_occ.u1y;
    delete [] ofStuff->tvl2_occ.u2x;
    delete [] ofStuff->tvl2_occ.u2y;
    delete [] ofStuff->tvl2_occ.v1;
    delete [] ofStuff->tvl2_occ.v2;
    delete [] ofStuff->tvl2_occ.rho_c1;
    delete [] ofStuff->tvl2_occ.rho_c2;
    delete [] ofStuff->tvl2_occ.grad;
    delete [] ofStuff->tvl2_occ.u1Aux;
    delete [] ofStuff->tvl2_occ.u2Aux;
    delete [] ofStuff->tvl2_occ.I2x;
    delete [] ofStuff->tvl2_occ.I2y;
    delete [] ofStuff->tvl2_occ.I2w;
    delete [] ofStuff->tvl2_occ.I2wx;
    delete [] ofStuff->tvl2_occ.I2wy;
    delete [] ofStuff->tvl2_occ.I0x;
    delete [] ofStuff->tvl2_occ.I0y;
    delete [] ofStuff->tvl2_occ.I0w;
    delete [] ofStuff->tvl2_occ.I0wx;
    delete [] ofStuff->tvl2_occ.I0wy;
    delete [] ofStuff->tvl2_occ.div_xi1;
    delete [] ofStuff->tvl2_occ.div_xi2;
    delete [] ofStuff->tvl2_occ.u_N;
}

//////////////////////////////////////////////////////////////
////TV-l2 COUPLED OPTICAL FLOW PROBLEM WITH OCCLUSIONS////////
/////////////////////////////////////////////////////////////
//Dual variable
static void tvl2coupled_getD(
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        float *u1x,
        float *u1y,
        float *u2x,
        float *u2y,
        float tau,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const int nx
        ){
    //Compute the value of xi
#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            const int i = l*nx + k;
            const float g11 = xi11[i]*xi11[i];
            const float g12 = xi12[i]*xi12[i];
            const float g21 = xi21[i]*xi21[i];
            const float g22 = xi22[i]*xi22[i];

            float xi_N = sqrt(g11 + g12 + g21 + g22);

            xi_N = MAX(1, xi_N);

            xi11[i] = (xi11[i] + tau*u1x[i])/xi_N;
            xi12[i] = (xi12[i] + tau*u1y[i])/xi_N;
            xi21[i] = (xi21[i] + tau*u2x[i])/xi_N;
            xi22[i] = (xi22[i] + tau*u2y[i])/xi_N;
        }
    }
}


//Primal variable
static void tvl2coupled_getP(
        float *u1,
        float *u2,
        float *v1,
        float *v2,
        float *div_xi1,
        float *div_xi2,
        float *u_N,
        int *mask,
        float theta,
        float tau,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const int nx,
        float *err
        ){
    float err_D = 0.0;

#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++){
            const int i = l*nx + k;
            const float u1k = u1[i];
            const float u2k = u2[i];

            //Only modify the inpainting domain
            // if (mask[i]==0){
            u1[i] = u1k - tau*(-div_xi1[i] + (u1k - v1[i])/theta);
            u2[i] = u2k - tau*(-div_xi2[i] + (u2k - v2[i])/theta);

            u_N[i]= (u1[i] - u1k) * (u1[i] - u1k) +
                    (u2[i] - u2k) * (u2[i] - u2k);
            // }else {
            //   u_N[i] = 0.0;
            // }
        }


    //Get the max val
    err_D = 0;
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            if (err_D < u_N[l*nx + k]){
                err_D = u_N[l*nx + k];
            }
        }
    }
    (*err) = err_D;

}

void eval_tvl2coupled_occ(
        float *I0,           // source image
        float *I1,           // target image
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_occ *tvl2_occ,
        float *ener_N,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        const float lambda,  // weight of the data term
        const float theta
        ){

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;


    //Columns and Rows
    const int nx = ofD->w;
    const int ny = ofD->h;

    //Optical flow derivatives
    float *v1   = tvl2_occ->v1;
    float *v2   = tvl2_occ->v2;
    float *u1x  = tvl2_occ->u1x;
    float *u1y  = tvl2_occ->u1y;
    float *u2x  = tvl2_occ->u2x;
    float *u2y  = tvl2_occ->u2y;

    float *I2w = tvl2_occ->I2w;

    float ener = 0.0;

    //forward_gradient_mixed_bound(u1,u1x,u1y,ii,ij,ei,ej,nx,ny);
    //forward_gradient_mixed_bound(u2,u2x,u2y,ii,ij,ei,ej,nx,ny);
    forward_gradient_patch(u1, u1x, u1y, ii, ij, ei, ej, nx);
    forward_gradient_patch(u2, u2x, u2y, ii, ij, ei, ej, nx);
    bicubic_interpolation_warp_patch(I1,  u1, u2, I2w,
                                     ii, ij, ei, ej, nx, ny, false);

    //Energy for all the patch. Maybe it useful only the 8 pixel around the seed.
    int m  = 0;
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            const int i = l*nx + k;
            float dt = lambda*fabs(I2w[i] - I0[i]);
            float dc = (1/(2*theta))*
                    ((u1[i] - v1[i])*(u1[i]- v1[i]) + (u2[i] - v2[i])*(u2[i] - v2[i]));
            float g1  = u1x[i]*u1x[i];
            float g12 = u1y[i]*u1y[i];
            float g21 = u2x[i]*u2x[i];
            float g2  = u2y[i]*u2y[i];
            float g  = sqrt(g1 + g12 + g21 + g2);
            if (!std::isfinite(dt)){
                std::printf("Datos corruptos\n");
            }
            if (!std::isfinite(g)){
                std::printf("Regularizacion corrupta\n");
            }
            ener += dc + dt + g;
            m++;
        }
    }
    ener /=(m*1.0);
    (*ener_N) = ener;
    assert(ener >= 0.0);
}

// Variational Optical flow method based on initial fixed values
// It minimizes the energy of \int_{B(x)} ||J(u)|| + |I_{1}(x+u)-I_{0}(x)|
// s.t u = u_0 for i.seeds
// J(u) = (u_x, u_y; v_x, v_y)
void guided_tvl2coupled_occ(
        float *I1,           // source image
        float *I2,           // target image
        float *I0,
        OpticalFlowData *ofD,
        Tvl2CoupledOFStuff_occ *tvl2_occ,
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
        ) {

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    float *u1_ba = ofD->u1_ba;
    float *u2_ba = ofD->u2_ba;
    int *mask = ofD->fixed_points;

    //Columns and Rows
    const int nx = ofD->w;
    const int ny = ofD->h;

    //Occlusion variable
    float *chi = tvl2_occ->chi;
    float *chix = tvl2_occ->chix;
    float *chiy = tvl2_occ->chiy;

    //Weigth
    float *g = tvl2_occ->g;

    //Dual variables
    float *xi11 = tvl2_occ->xi11;
    float *xi12 = tvl2_occ->xi12;
    float *xi21 = tvl2_occ->xi21;
    float *xi22 = tvl2_occ->xi22;

    //Derivatives of u1,u2
    float *u1x = tvl2_occ->u1x;
    float *u1y = tvl2_occ->u1y;
    float *u2x = tvl2_occ->u2x;
    float *u2y = tvl2_occ->u2y;

    float *v1 = tvl2_occ->v1;
    float *v2 = tvl2_occ->v2;

    float *rho_c1 = tvl2_occ->rho_c1;
    float *rho_c2 = tvl2_occ->rho_c2;
    float *grad = tvl2_occ->grad;

    float *u1Aux = tvl2_occ->u1Aux;
    float *u2Aux = tvl2_occ->u2Aux;

    float *I2x = tvl2_occ->I2x;
    float *I2y = tvl2_occ->I2y;
    float *I2w = tvl2_occ->I2w;
    float *I2wx = tvl2_occ->I2wx;
    float *I2wy = tvl2_occ->I2wy;

    float *I0x = tvl2_occ->I0x;
    float *I0y = tvl2_occ->I0y;
    float *I0w = tvl2_occ->I0w;
    float *I0wx = tvl2_occ->I0wx;
    float *I0wy = tvl2_occ->I0wy;

    float *div_xi1 = tvl2_occ->div_xi1;
    float *div_xi2 = tvl2_occ->div_xi2;
    float *u_N = tvl2_occ->u_N;

    const float l_t = lambda * theta;

    //Initialization of dual variables and backward flow
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            const int  i = l*nx + k;            
            xi11[i] = xi12[i] = xi21[i] = xi22[i] = 0.0;

            u1_ba[i] = -u1[i];
            u2_ba[i] = -u2[i];
        }
    }

    for (int warpings = 0; warpings < warps; warpings++) {
        // Compute the warping of I2 and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp_patch(I2,  u1, u2, I2w, ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I2x, u1, u2, I2wx, ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I2y, u1, u2, I2wy, ii, ij, ei, ej, nx, ny, false);

        // Compute the warping of I0 and its derivatives Ir(x - u1o), Irx (x - u1o) and Iry (x - u2o)
        bicubic_interpolation_warp_patch(I0,  u1_ba, u2_ba, I0w, ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I0x, u1_ba, u2_ba, I0wx, ii, ij, ei, ej, nx, ny, false);
        bicubic_interpolation_warp_patch(I0y, u1_ba, u2_ba, I0wy, ii, ij, ei, ej, nx, ny, false);
//////Continue here/////
        //Compute values that will not change during the whole wraping
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for (int l = ij; l < ej; l++)
            for (int k = ii; k < ei; k++){

                const int i = l*nx + k;

                const float Ix2 = I1wx[i] * I1wx[i];
                const float Iy2 = I1wy[i] * I1wy[i];

                // store the |Grad(I2)|^2
                grad[i] = (Ix2 + Iy2);

                // Compute the constant part of the rho function
                rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                            - I1wy[i] * u2[i] - I0[i]);
            }

#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for (int l = ij; l < ej; l++){
            for (int k = ii; k < ei; k++){
                const int i = l*nx + k;
                u1_[i] = u1[i];
                u2_[i] = u2[i];
            }
        }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS_OF){

            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
            for (int l = ij; l < ej; l++){
                for (int k = ii; k < ei; k++){
                    const int i = l*nx + k;
                    const float rho = rho_c[i]
                            + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
                    float d1, d2;

                    if (rho < - l_t * grad[i]){
                        d1 = l_t * I1wx[i];
                        d2 = l_t * I1wy[i];
                    }
                    else{
                        if (rho > l_t * grad[i]){
                            d1 = -l_t * I1wx[i];
                            d2 = -l_t * I1wy[i];
                        }
                        else{
                            // if gradient is too small, we treat it as zero
                            if (grad[i] < GRAD_IS_ZERO)
                                d1 = d2 = 0;
                            else{
                                float fi = -rho/grad[i];
                                d1 = fi * I1wx[i];
                                d2 = fi * I1wy[i];
                            }
                        }
                    }
                    v1[i] = u1[i] + d1;
                    v2[i] = u2[i] + d2;
                }
            }
            // estimate the values of the variable (u1, u2)

            //Compute dual variables
            forward_gradient_patch(u1_, u1x, u1y, ii, ij, ei, ej, nx);
            forward_gradient_patch(u2_, u2x, u2y, ii, ij, ei, ej, nx);
            tvl2coupled_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y,
                             tau, ii, ij, ei, ej, nx);

            //Primal variables
            divergence_patch(xi11, xi12, div_xi1, ii, ij, ei, ej, nx);
            divergence_patch(xi21, xi22, div_xi2, ii, ij, ei, ej, nx);

            //Save previous iteration
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
            for (int l = ij; l < ej; l++){
                for (int k = ii; k < ei; k++){
                    const int i = l*nx + k;
                    u1Aux[i] = u1[i];
                    u2Aux[i] = u2[i];
                }
            }

            tvl2coupled_getP(u1, u2, v1, v2, div_xi1, div_xi2, u_N,
                             mask, theta, tau, ii, ij, ei, ej, nx, &err_D);

            //(aceleration = 1);
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
            for (int l = ij; l < ej; l++){
                for (int k = ii; k < ei; k++){
                    const int i = l*nx + k;
                    u1_[i] = 2*u1[i] - u1Aux[i];
                    u2_[i] = 2*u2[i] - u2Aux[i];
                }
            }


        }
        if (verbose)
            std::printf("Warping: %d, Iter: %d "
                        "Error: %f\n", warpings, n, err_D);
    }
    eval_tvl2coupled_occ(I1, I2, ofD, tvl2_occ, ener_N, ii, ij, ei, ej, lambda, theta);
}

#endif //TVL2-L1 functional
