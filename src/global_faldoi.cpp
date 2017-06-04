// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license athis program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2013, Roberto P.Palomares <roberto.palomares@upf.edu>
// All rights reserved.

#ifndef GLOBAL_FALDOI
#define GLOBAL_FALDOI

#include <cmath>
#include <cstdio>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <vector>
#include <algorithm> 

extern "C" {
//#include "mask.h"
#include "bicubic_interpolation.h"
#include "iio.h"
}

#include "tvl2_model_occ.h"
#include "utils.h"


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>

using namespace std;

#define GRAD_IS_ZERO_GLOBAL 1E-10
#define MAX_ITERATIONS 400

#define PRESMOOTHING_SIGMA  0.90
#define PAR_DEFAULT_NPROC   0    //0
#define PAR_DEFAULT_LAMBDA  40//40
#define PAR_DEFAULT_THETA   0.1
#define PAR_DEFAULT_TAU     0.125 //0.25
#define PAR_DEFAULT_NWARPS  5  //5
#define PAR_DEFAULT_TOL_D   0.01
#define PAR_DEFAULT_NPROC   0    //0
#define PAR_DEFAULT_VERBOSE 1  //0

#define M_TVL1       0
#define M_TVL1_W     1
#define M_NLTVL1     2
#define M_NLTVL1_W   3
#define M_TVCSAD     4
#define M_TVCSAD_W   5
#define M_NLTVCSAD   6
#define M_NLTVCSAD_W 7
#define M_TVL1_OCC 8


//#ifndef DISABLE_OMP
//#include <omp.h>
//#endif//DISABLE_OMP



#define MAX(x,y) ((x)>(y)?(x):(y))


//////////
///////////WARNING
/**
 *
 * Function to compute the optical flow in one scale
 *
 **/
void Dual_TVL1_optic_flow(
        float *I0,           // source image
        float *I1,           // target image
        float *u1,           // x component of the optical flow
        float *u2,           // y component of the optical flow
        const int   nx,      // image width
        const int   ny,      // image height
        const float tau,     // time step
        const float lambda,  // weight parameter for the data term
        const float theta,   // weight parameter for (u - v)²
        const int   warps,   // number of warpings per scale
        const float epsilon, // tolerance for numerical convergence
        const bool  verbose  // enable/disable the verbose mode
        ) {
    const int   size = nx * ny;
    const float l_t = lambda * theta;

    float *I1x    = new float[size];
    float *I1y    = new float[size];
    float *I1w    = new float[size];
    float *I1wx   = new float[size];
    float *I1wy   = new float[size];
    float *rho_c  = new float[size];
    float *v1     = new float[size];
    float *v2     = new float[size];
    float *p11    = new float[size];
    float *p12    = new float[size];
    float *p21    = new float[size];
    float *p22    = new float[size];
    float *div    = new float[size];
    float *grad   = new float[size];
    float *div_p1 = new float[size];
    float *div_p2 = new float[size];
    float *u1x    = new float[size];
    float *u1y    = new float[size];
    float *u2x    = new float[size];
    float *u2y    = new float[size];

    centered_gradient(I1, I1x, I1y, nx, ny);

    // initialization of p
    for (int i = 0; i < size; i++) {
        p11[i] = p12[i] = 0.0;
        p21[i] = p22[i] = 0.0;
    }

    for (int warpings = 0; warpings < warps; warpings++) {
        // compute the warping of the target image and its derivatives
        bicubic_interpolation_warp(I1,  u1, u2, I1w,  nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

        //#pragma omp parallel for
        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        int n = 0;
        float error = INFINITY;
        while (error > epsilon * epsilon && n < MAX_ITERATIONS) {
            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++) {
                const float rho = rho_c[i]
                        + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);

                float d1, d2;

                if (rho < - l_t * grad[i]) {
                    d1 = l_t * I1wx[i];
                    d2 = l_t * I1wy[i];
                } else {
                    if (rho > l_t * grad[i]) {
                        d1 = -l_t * I1wx[i];
                        d2 = -l_t * I1wy[i];
                    } else {
                        if (grad[i] < GRAD_IS_ZERO)
                            d1 = d2 = 0;
                        else {
                            float fi = -rho/grad[i];
                            d1 = fi * I1wx[i];
                            d2 = fi * I1wy[i];
                        }
                    }
                }

                v1[i] = u1[i] + d1;
                v2[i] = u2[i] + d2;
            }

            // compute the divergence of the dual variable (p1, p2)
            divergence(p11, p12, div_p1, nx ,ny);
            divergence(p21, p22, div_p2, nx ,ny);

            // estimate the values of the optical flow (u1, u2)
            error = 0.0;
            //#pragma omp parallel for reduction(+:error)
            for (int i = 0; i < size; i++) {
                const float u1k = u1[i];
                const float u2k = u2[i];

                u1[i] = v1[i] + theta * div_p1[i];
                u2[i] = v2[i] + theta * div_p2[i];

                error += (u1[i] - u1k) * (u1[i] - u1k) +
                        (u2[i] - u2k) * (u2[i] - u2k);
            }
            error /= size;

            // compute the gradient of the optical flow (Du1, Du2)
            forward_gradient(u1, u1x, u1y, nx ,ny);
            forward_gradient(u2, u2x, u2y, nx ,ny);

            // estimate the values of the dual variable (p1, p2)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++) {
                const float taut = tau / theta;
                const float g1   = hypotf(u1x[i], u1y[i]);
                const float g2   = hypotf(u2x[i], u2y[i]);
                const float ng1  = 1.0 + taut * g1;
                const float ng2  = 1.0 + taut * g2;

                p11[i] = (p11[i] + taut * u1x[i]) / ng1;
                p12[i] = (p12[i] + taut * u1y[i]) / ng1;
                p21[i] = (p21[i] + taut * u2x[i]) / ng2;
                p22[i] = (p22[i] + taut * u2y[i]) / ng2;
            }
        }

        if (verbose)
            fprintf(stderr, "Warping: %d, "
                            "Iterations: %d, "
                            "Error: %f\n", warpings, n, error);
    }

    // delete allocated memory
    free(I1x);
    free(I1y);
    free(I1w);
    free(I1wx);
    free(I1wy);
    free(rho_c);
    free(v1);
    free(v2);
    free(p11);
    free(p12);
    free(p21);
    free(p22);
    free(div);
    free(grad);
    free(div_p1);
    free(div_p2);
    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);
}




////////////WARNING////////
///////////////////////////////////
/*
 * - Name: getP_Du

 * - Output: float *u - New optical flow estimated
 *
*/
void ofDu_getP(
        float *u1,
        float *u2,
        float *v1,
        float *v2,
        float *div_xi1,
        float *div_xi2,
        float *u_N,
        float theta,
        float tau,
        int size,
        float *err
        ){
    float err_D = 0.0;
    float min,max;

    //#pragma omp parallel for reduction(+:err_D)
    for (int i = 0; i < size; i++){

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k  -tau*(-div_xi1[i]  + (u1k - v1[i])/theta);
        u2[i] = u2k  -tau*(-div_xi2[i]  + (u2k - v2[i])/theta);

        u_N[i]= (u1[i] - u1k) * (u1[i] - u1k) +
                (u2[i] - u2k) * (u2[i] - u2k);
    }

    getminmax(&min,&max,u_N,size);

    err_D =max;
    (*err) = err_D;
}


/*
 * - Name: getD_Du

 *
*/
void ofDu_getD(
        float *xi11,
        float *xi12,
        float *xi22,
        float *u1x,
        float *u1y,
        float *u2x,
        float *u2y,
        float tau,
        int size
        ){
    //#pragma omp parallel for

    for (int i = 0; i < size; i++)
    {

        const float g11 = xi11[i]*xi11[i];
        const float g12 = xi12[i]*xi12[i];
        const float g22 = xi22[i]*xi22[i];

        float xi_N = sqrt(g11 + g22 + 2*g12);

        xi_N = MAX(1,xi_N);

        xi11[i] = (xi11[i] + tau*u1x[i])/xi_N;
        xi12[i] = (xi12[i] + 0.5*tau*(u1y[i] + u2x[i]))/xi_N;
        xi22[i] = (xi22[i] + tau*u2y[i])/xi_N;
    }
}


/*
 * - Name: getP_Du

 * - Output: float *u - New optical flow estimated
 *
*/
void ofTVl2_getP(
        float *u1,
        float *u2,
        float *v1,
        float *v2,
        float *div_xi1,
        float *div_xi2,
        float *u_N,
        float theta,
        float tau,
        int size,
        float *err
        )
{
    float err_D = 0.0;
    float min,max;

#pragma omp parallel for
    for (int i = 0; i < size; i++){

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k  -tau*(-div_xi1[i]  + (u1k - v1[i])/theta);
        u2[i] = u2k  -tau*(-div_xi2[i]  + (u2k - v2[i])/theta);

        u_N[i]= (u1[i] - u1k) * (u1[i] - u1k) +
                (u2[i] - u2k) * (u2[i] - u2k);
    }

    getminmax(&min,&max,u_N,size);

    err_D =max;
    (*err) = err_D;
}

/*
 * - Name: ofTVl2_getD

 *
*/
void ofTVl2_getD(
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        float *u1x,
        float *u1y,
        float *u2x,
        float *u2y,
        float tau,
        int size
        ){

#pragma omp parallel for
    for (int i = 0; i < size; i++)
    {

        const float g11 = xi11[i]*xi11[i];
        const float g12 = xi12[i]*xi12[i];
        const float g21 = xi21[i]*xi21[i];
        const float g22 = xi22[i]*xi22[i];

        float xi_N = sqrt(g11 + g12 + g21 + g22);

        xi_N = MAX(1,xi_N);

        xi11[i] = (xi11[i] + tau*u1x[i])/xi_N;
        xi12[i] = (xi12[i] + tau*u1y[i])/xi_N;
        xi21[i] = (xi21[i] + tau*u2x[i])/xi_N;
        xi22[i] = (xi22[i] + tau*u2y[i])/xi_N;
    }
}


void duOF(
        float *I0,           // source image
        float *I1,           // target image
        float *u1,      // x component of the optical flow
        float *u2,      // y component of the optical flow
        float *xi11,
        float *xi12,
        float *xi22,
        const float lambda,  // weight of the data term
        const float theta,   // weight of the data term
        const float tau,     // time step
        const float tol_OF,  // tol max allowed
        const int   nx,      // image width
        const int   ny,      // image height
        const int   warps,   // number of warpings per scale
        const bool  verbose  // enable/disable the verbose mode
        )
{

    const float l_t = lambda * theta;
    const int   size = nx * ny;

    size_t sd = sizeof(float);


    float *u1x    = new float[size];
    float *u1y    = new float[size];
    float *u2x    = new float[size];
    float *u2y    = new float[size];

    float *v1     = new float[size];
    float *v2     = new float[size];

    float *rho_c  = new float[size];
    float *grad   = new float[size];

    float *u1_     = new float[size];
    float *u2_     = new float[size];

    float *u1Aux   = new float[size];
    float *u2Aux   = new float[size];

    float *I1x    = new float[size];
    float *I1y    = new float[size];

    float *I1w    = new float[size];
    float *I1wx   = new float[size];
    float *I1wy   = new float[size];

    //Divergence
    float *div_xi1 = new float[size];
    float *div_xi2 = new float[size];

    float *u_N = new float[size];

    centered_gradient(I1, I1x, I1y, nx, ny);

    for (int warpings = 0; warpings < warps; warpings++)
    {
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1,  u1, u2, I1w,  nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

        for (int i = 0; i < size; i++)
        {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        for (int i = 0; i < nx*ny; i++){
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS)
        {

            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++)
            {
                const float rho = rho_c[i]
                        + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
                float d1, d2;

                if (rho < - l_t * grad[i])
                {
                    d1 = l_t * I1wx[i];
                    d2 = l_t * I1wy[i];
                }
                else
                {
                    if (rho > l_t * grad[i])
                    {
                        d1 = -l_t * I1wx[i];
                        d2 = -l_t * I1wy[i];
                    }
                    else
                    {
                        if (grad[i] < GRAD_IS_ZERO)
                            d1 = d2 = 0;
                        else
                        {
                            float fi = -rho/grad[i];
                            d1 = fi * I1wx[i];
                            d2 = fi * I1wy[i];
                        }
                    }
                }

                v1[i] = u1[i] + d1;
                v2[i] = u2[i] + d2;
            }

            //Dual variables
            forward_gradient(u1_,u1x,u1y,nx,ny);
            forward_gradient(u2_,u2x,u2y,nx,ny);
            ofDu_getD(xi11,xi12,xi22,u1x,u1y,u2x,u2y,tau,size);

            //Primal variables
            divergence(xi11,xi12,div_xi1,nx,ny);
            divergence(xi12,xi22,div_xi2,nx,ny);

            //Almacenamos la iteracion anterior
            for (int i = 0; i < size; i++){
                u1Aux[i] = u1[i];
                u2Aux[i] = u2[i];
            }

            ofDu_getP(u1,u2,v1,v2,div_xi1,div_xi2,u_N,theta,tau,size,&err_D);

            //(aceleration = 1);
            for (int i = 0; i < size; i++){
                u1_[i] = 2*u1[i] - u1Aux[i];
                u2_[i] = 2*u2[i] - u2Aux[i];
            }



        }
        if (verbose)
            fprintf(stderr, "Warping: %d,Iter: %d "
                            "Error: %f\n", warpings,n, err_D);
    }

    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);

    free(v1);
    free(v2);

    free(rho_c);
    free(grad);

    free(u1_);
    free(u2_);

    free(u1Aux);
    free(u2Aux);

    free(I1x);
    free(I1y);

    free(I1w);
    free(I1wx);
    free(I1wy);

    free(div_xi1);
    free(div_xi2);

    free(u_N);
}

void tvl2OF(
        float *I0,           // source image
        float *I1,           // target image
        float *u1,           // x component of the optical flow
        float *u2,           // y component of the optical flow
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        const float lambda,  // weight of the data term
        const float theta,   // weight of the data term
        const float tau,     // time step
        const float tol_OF,  // tol max allowed
        const int   nx,      // image width
        const int   ny,      // image height
        const int   warps,   // number of warpings per scale
        const bool  verbose  // enable/disable the verbose mode
        ) {

    const float l_t = lambda * theta;
    const int   size = nx * ny;


    float *u1x    = new float[size];
    float *u1y    = new float[size];
    float *u2x    = new float[size];
    float *u2y    = new float[size];

    float *v1     = new float[size];
    float *v2     = new float[size];

    float *rho_c  = new float[size];
    float *grad  = new float[size];

    float *u1_     = new float[size];
    float *u2_     = new float[size];

    float *u1Aux   = new float[size];
    float *u2Aux   = new float[size];

    float *I1x    = new float[size];
    float *I1y    = new float[size];

    float *I1w    = new float[size];
    float *I1wx   = new float[size];
    float *I1wy   = new float[size];

    //Divergence
    float *div_xi1 = new float[size];
    float *div_xi2 = new float[size];

    float *u_N = new float[size];

    centered_gradient(I1, I1x, I1y, nx, ny);

    for (int warpings = 0; warpings < warps; warpings++)
    {
        //printf("warpings:%d\n", warpings);
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1,  u1, u2, I1w,  nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);

#pragma omp parallel for
        for (int i = 0; i < size; i++)
        {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        memcpy(u1_, u1,size*sizeof(float));
        memcpy(u2_, u2,size*sizeof(float));
        // for (int i = 0; i < nx*ny; i++){
        //   u1_[i] = u1[i];
        //   u2_[i] = u2[i];
        // }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS)
        {

            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
#pragma omp parallel for
            for (int i = 0; i < size; i++)
            {
                const float rho = rho_c[i]
                        + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
                float d1, d2;

                if (rho < - l_t * grad[i])
                {
                    d1 = l_t * I1wx[i];
                    d2 = l_t * I1wy[i];
                }
                else{
                    if (rho > l_t * grad[i])
                    {
                        d1 = -l_t * I1wx[i];
                        d2 = -l_t * I1wy[i];
                    }
                    else{
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

            //Dual variables
            forward_gradient(u1_, u1x, u1y, nx, ny);
            forward_gradient(u2_, u2x, u2y, nx, ny);
            ofTVl2_getD(xi11, xi12, xi21, xi22, u1x, u1y, u2x, u2y, tau, size);

            //Primal variables
            divergence(xi11, xi12, div_xi1, nx, ny);
            divergence(xi21, xi22, div_xi2, nx, ny);

            //Almacenamos la iteracion anterior
            memcpy(u1Aux, u1, size*sizeof(float));
            memcpy(u2Aux, u2, size*sizeof(float));
            // for (int i = 0; i < size; i++){
            //   u1Aux[i] = u1[i];
            //   u2Aux[i] = u2[i];
            // }

            ofTVl2_getP(u1, u2, v1, v2, div_xi1, div_xi2, u_N, theta, tau, size, &err_D);

            //(aceleration = 1);
            for (int i = 0; i < size; i++){
                u1_[i] = 2*u1[i] - u1Aux[i];
                u2_[i] = 2*u2[i] - u2Aux[i];

            }

        }
        if (verbose)
            fprintf(stderr, "Warping: %d,Iter: %d "
                            "Error: %f\n", warpings, n, err_D);
    }

    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);

    free(v1);
    free(v2);

    free(rho_c);
    free(grad);

    free(u1_);
    free(u2_);

    free(u1Aux);
    free(u2Aux);

    free(I1x);
    free(I1y);

    free(I1w);
    free(I1wx);
    free(I1wy);

    free(div_xi1);
    free(div_xi2);

    free(u_N);
}

////////////////////////////////////NLTVL1//////////////////////////////////////
#define MAX_SPATIAL 2
#define MAX_INTENSITY 5
#define MAX_BETA  2//Neighboor
#define MAX_DUAL_VAR (2*MAX_BETA + 1)*(2*MAX_BETA + 1) -1 // 5x5

//Struct
struct DualVariables_global{
    float sc[MAX_DUAL_VAR]; // value of p(x,y)
    float wp[MAX_DUAL_VAR]; // weight of non local
    int   ap[MAX_DUAL_VAR]; //absolute positon of p(y,x)
    int   rp[MAX_DUAL_VAR]; //relative position of p(y,x) in the structure
    float wt = 0.0;
};

inline bool positive(int val) {
    return val >= 0;
}

float aux_pow2( float f ) {return f*f;}
//Asume que las imagenes no estan normalizadas
void image_to_lab(float *in, int size, float *out) {
    const float T = 0.008856;
    const float color_attenuation = 1.5f;
    for(int i = 0 ; i < size ; i++){
        const float r = in[i]/255.f;
        const float g = in[i + size]/255.f;
        const float b = in[i + 2*size]/255.f;
        float X = 0.412453 * r + 0.357580 * g + 0.180423 * b;
        float Y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
        float Z = 0.019334 * r + 0.119193 * g + 0.950227 * b;
        X /= 0.950456;
        Z /= 1.088754;
        float Y3 = pow(Y,1./3);
        float fX = X>T ? pow(X,1./3) : 7.787 * X + 16/116.;
        float fY = Y>T ? Y3 : 7.787 * Y + 16/116.;
        float fZ = Z>T ? pow(Z,1./3) : 7.787 * Z + 16/116.;
        float L = Y>T ? 116 * Y3 - 16.0 : 903.3 * Y;
        float A = 500 * (fX - fY);
        float B = 200 * (fY - fZ);
        // correct L*a*b*: dark area or light area have less reliable colors
        float correct_lab = exp(-color_attenuation*aux_pow2(aux_pow2(L/100) - 0.6));
        out[i] = L;
        out[i + size] = A*correct_lab;
        out[i + 2*size] = B*correct_lab;
    }
}




static int validate_ap_2(int w, int h, int i, int j, int di, int dj){
    const int r = j + dj; //Row
    const int c = i + di; //Colum
    if ( c < 0 || c >= w || r < 0 || r >= h)
        return -1;
    return r*w + c;
}

static float get_wspatial_2( int l, int k) {
    float ws = MAX_BETA;
    float w_tmp;
    float difS = 0.0;
    difS = hypot(l,k);
    w_tmp = exp(- difS/ws);
    //std::printf("Dif_S: %f W_S: %f\n",difS, w_tmp);

    return w_tmp;
}

static float get_wcolor_2(
        float *a,
        int w,
        int h,
        int i,
        int j,
        int l,
        int k,
        int pd
        ) {
    float wi = MAX_INTENSITY;
    float w_tmp;
    float difI = 0.0;
    for (int m = 0; m < pd; m++){
        float aux = getsample_0(a, w, h, pd, i, j, m)
                - getsample_0(a, w, h, pd, i + l, j + k, m);
        difI += aux*aux;
        //std::printf("valI_%d:%f ",m, difI);
    }
    //std::printf("\n");
    difI = sqrt(difI);
    w_tmp = exp( -difI/wi);

    return w_tmp;
}

static float get_weight_2(
        float *a,
        int w,
        int h,
        int i,
        int j,
        int l,
        int k,
        int pd
        ) {
    float wc = get_wcolor_2(a, w, h, i, j, l, k, pd);
    float ws = get_wspatial_2(l, k);

    return wc*ws;
}

void initialize_dual_variables(
        float *a,
        const int pd,
        const int w,
        const int h,
        const int n_d,
        const int radius,
        DualVariables_global *p,
        DualVariables_global *q
        ) {
    int size = w*h;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < n_d; j++) {
            p[i].sc[j] = -2.0;
            p[i].wp[j] = -2.0;
            p[i].ap[j] = -1; //Indicates that is out
            p[i].rp[j] = -1; //Indicates that it is out.

            q[i].sc[j] = -2.0;
            q[i].wp[j] = -2.0;
            q[i].ap[j] = -1; //Indicates that is out
            q[i].rp[j] = -1; //Indicates that it is out.
        }


    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int   it = 0;
            float ne = 0.0;
            const int pos = j*w + i;
            for (int k = -radius; k < (radius+1); k++)
                for (int l = -radius; l < (radius+1); l++) {
                    //std::printf("OUT Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d\n",radius, j,i,k,l,it);
                    if (!(k == 0 && k == l)) {

                        int ap = validate_ap_2(w, h, i, j, l, k);
                        if (positive(ap)) {
                            p[pos].sc[it] = q[pos].sc[it] = 0.0;
                            p[pos].ap[it] = q[pos].ap[it] = ap;
                            p[pos].rp[it] = q[pos].rp[it] = n_d - (it + 1);
                            assert(p[pos].rp[it] >= 0);
                            assert(q[pos].rp[it] >= 0);

                            //Compute the weigh
                            float wp = sqrt(get_weight_2(a, w, h, i, j, l, k, pd));
                            p[pos].wp[it] = q[pos].wp[it] = wp;

                            ne += wp;
                        }
                        it++;
                        //std::printf("Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d Ne:%d\n",radius, j,i,k,l,it,ne);
                    }
                }
            //TODO: It is used to normalize
            p[pos].wt = ne;
            q[pos].wt = ne;
        }
    // std::printf(" Acaba\n");
}

void non_local_divergence(
        DualVariables_global *p,
        int size,
        int n_d,
        float *div_p
        ) {

    //#pragma omp parallel for
    for (int i = 0; i < size; i++){
        div_p[i] = 0.0;
        for (int j = 0; j < n_d; j++){
            const int ap = p[i].ap[j];
            const int rp = p[i].rp[j];
            if (positive(ap)) {
                assert (p[i].rp[j]>=0);
                const float pxy = p[i].sc[j];
                const float pyx = p[ap].sc[rp];
                const float w = p[i].wp[j];
                div_p[i] += w*(pxy - pyx);
            }
        }
        div_p[i] /=p[i].wt;
    }
}



//Auxiliar Chambolle Scheme functions

/*
 * - Name: getP

 *
*/
void ofnltv_getP(
        float *v1,
        float *v2,
        float *div_p1,
        float *div_p2,
        float theta,
        float tau,
        int size,
        float *u1,
        float *u2,
        float *err
        ) {
    float err_D = 0.0;

    //#pragma omp parallel for reduction(+:err_D)
    for (int i = 0; i < size; i++){

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k  -tau*(div_p1[i]  + (u1k - v1[i])/theta);
        u2[i] = u2k  -tau*(div_p2[i]  + (u2k - v2[i])/theta);

        err_D += (u1[i] - u1k) * (u1[i] - u1k) +
                (u2[i] - u2k) * (u2[i] - u2k);
    }
    err_D /= size;
    (*err) = err_D;
}

/*
 * - Name: getD

 *
*/
void ofnltv_getD(
        float *u1,
        float *u2,
        int size,
        int n_d,
        float tau,
        DualVariables_global *p1,
        DualVariables_global *p2
        )
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < n_d; j++)
        {
            const int ap1 = p1[i].ap[j];
            // const int rp1 = p1[i].rp[j];
            const float wt1 = p1[i].wt;

            const int ap2 = p2[i].ap[j];
            // const int rp2 = p2[i].rp[j];
            const float wt2 = p2[i].wt;

            if (positive(ap1) && positive(ap2))
            {
                // assert(rp1 >=0);
                const float w1 = p1[i].wp[j];
                const float u1x = u1[i];
                const float u1y = u1[ap1];
                const float nlgr1 =  w1 * (u1x -u1y)/wt1;
                const float nl1 = sqrt(nlgr1*nlgr1);
                const float nl1g = 1 + tau * nl1;

                p1[i].sc[j] =  (p1[i].sc[j] + tau *nlgr1)/nl1g;
            }

            if (positive(ap1) && positive(ap2))
            {
                // assert(rp2 >=0);
                const float w2 = p2[i].wp[j];
                const float u2x = u2[i];
                const float u2y = u2[ap2];
                const float nlgr2 =  w2 * (u2x -u2y)/wt2;
                const float nl2 = sqrt(nlgr2*nlgr2);
                const float nl2g = 1 + tau * nl2;

                p2[i].sc[j] =  (p2[i].sc[j] + tau *nlgr2)/nl2g;

            }
        }
}


void nltvl1_PD(
        float *I0,           // source image
        float *I1,           // target image
        float *a,            // source image (color)
        int pd,               //Channel numer
        const float lambda,  // weight of the data term
        const float theta,   // weight of the data term
        const float tau,     // time step
        const float tol_OF,  // tol max allowed
        const int   w,      // image width
        const int   h,      // image height
        const int   warps,   // number of warpings per scale
        const bool  verbose, // enable/disable the verbose mode
        float *u1,       // x component of the optical flow
        float *u2      // y component of the optical flow
        )
{
    const int   size = w * h;
    const float l_t = lambda * theta;

    DualVariables_global *p = new DualVariables_global[size];
    DualVariables_global *q = new DualVariables_global[size];
    float *v1     = new float[size];
    float *v2     = new float[size];
    float *rho_c  = new float[size];
    float *grad   = new float[size];
    float *u1_    = new float[size];
    float *u2_    = new float[size];
    float *u1_tmp = new float[size];
    float *u2_tmp = new float[size];
    float *I1x    = new float[size];
    float *I1y    = new float[size];
    float *I1w    = new float[size];
    float *I1wx   = new float[size];
    float *I1wy   = new float[size];
    float *div_p  = new float[size];
    float *div_q  = new float[size];

    int radius = MAX_BETA;
    int n_d = MAX_DUAL_VAR;


    std::printf("Antes\n");
    //Initialization of the Dual variables.
    initialize_dual_variables(a, pd, w, h, n_d, radius, p, q);
    centered_gradient(I1, I1x, I1y, w, h);

    std::printf("Inicializadon\n");
    for (int warpings = 0; warpings < warps; warpings++)
    {
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1,  u1, u2, I1w,  w, h, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, w, h, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, w, h, true);
        //#pragma omp parallel for
        for (int i = 0; i < size; i++)
        {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        for (int i = 0; i < size; i++)
        {
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        // while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS)
        while (n < MAX_ITERATIONS)
        {

            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++)
            {
                const float rho = rho_c[i]
                        + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);
                float d1, d2;

                if (rho < - l_t * grad[i])
                {
                    d1 = l_t * I1wx[i];
                    d2 = l_t * I1wy[i];
                }
                else
                {
                    if (rho > l_t * grad[i])
                    {
                        d1 = -l_t * I1wx[i];
                        d2 = -l_t * I1wy[i];
                    }
                    else
                    {
                        if (grad[i] < GRAD_IS_ZERO)
                            d1 = d2 = 0;
                        else
                        {
                            float fi = -rho/grad[i];
                            d1 = fi * I1wx[i];
                            d2 = fi * I1wy[i];
                        }
                    }
                }

                v1[i] = u1[i] + d1;
                v2[i] = u2[i] + d2;
            }
            //Dual variables
            ofnltv_getD(u1_, u2_, size, n_d, tau, p, q);
            //Almacenamos la iteracion anterior
            for (int i = 0; i < size; i++)
            {
                u1_tmp[i] = u1[i];
                u2_tmp[i] = u2[i];
            }

            //Primal variables
            non_local_divergence(p, size, n_d, div_p);
            non_local_divergence(q, size, n_d, div_q);
            ofnltv_getP(v1, v2, div_p, div_q, theta, tau, size, u1, u2, &err_D);

            //(aceleration = 1);
            for (int i = 0; i < size; i++)
            {
                u1_[i] = 2*u1[i] - u1_tmp[i];
                u2_[i] = 2*u2[i] - u2_tmp[i];
            }

        }
        if (verbose)
            std::printf("Warping: %d,Iter: %d Error: %f\n", warpings,n, err_D);
    }

    delete [] v1;
    delete [] v2;

    delete [] rho_c;
    delete [] grad;

    delete [] u1_;
    delete [] u2_;

    delete [] u1_tmp;
    delete [] u2_tmp;

    delete [] I1x;
    delete [] I1y;

    delete [] I1w;
    delete [] I1wx;
    delete [] I1wy;

    delete [] div_p;
    delete [] div_q;
}
//////////////////////////////TV-CSAD///////////////////////////////////////////
#define DT_R  3 //Neighboor 7x7
#define DT_NEI (2*DT_R + 1)*(2*DT_R + 1) -1 // 5x5

////Struct
//struct PosNei {
//    int   api[DT_NEI];
//    int   apj[DT_NEI];
//    float b[DT_NEI];
//    std::vector<float>  ba;
//    int n;
//};

////////////////////////////////////////////////////////////////////////////////
void initialize_pos_nei(
        const int w,
        const int h,
        const int n_d,
        const int radius,
        PosNei *p
        )
{ 
    int size = w*h;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < n_d; j++)
        {
            p[i].api[j] = -1; //Indicates that is out
            p[i].apj[j] = -1;
            p[i].b[j] = 0.0;
        }


    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            int   it = 0;
            float ne = 0.0;
            const int pos = j*w + i;
            for (int k = -radius; k < (radius+1); k++)
                for (int l = -radius; l < (radius+1); l++)
                {
                    //std::printf("OUT Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d\n",radius, j,i,k,l,it);
                    if (!(k ==0 && k==l))
                    {

                        int ap = validate_ap_2(w, h, i, j, l, k);
                        if (positive(ap))
                        {
                            p[pos].api[it] = i + l;
                            p[pos].apj[it] = j + k;
                            assert(p[pos].api[it]>=0 && p[pos].apj[it]>=0);
                            ne++;
                        }
                        it++;
                        //std::printf("Radius: %d j:%d, i:%d , k:%d l:%d Iter:%d Ne:%d\n",radius, j,i,k,l,it,ne);
                    }
                }
            p[pos].n = ne;
            for (int k = 0; k < (2*p[pos].n + 1); k++){
                p[pos].ba.push_back(0.0);
            }
        }
}


//Auxiliar Chambolle Scheme functions

//Chambolle functions
/*
 * - Name: getP_Du

 * - Output: float *u - New optical flow estimated
 *
*/
void tvcsad_getP(float *u1,float *u2,float *v1, float *v2,float *div_xi1, float *div_xi2, float *u_N, float theta, float tau, int size, float *err)
{
    float err_D = 0.0;
    float min, max;

    //#pragma omp parallel for reduction(+:err_D)
    for (int i = 0; i < size; i++){

        const float u1k = u1[i];
        const float u2k = u2[i];

        u1[i] = u1k  -tau*(-div_xi1[i]  + (u1k - v1[i])/theta);
        u2[i] = u2k  -tau*(-div_xi2[i]  + (u2k - v2[i])/theta);

        err_D += (u1[i] - u1k) * (u1[i] - u1k) +
                (u2[i] - u2k) * (u2[i] - u2k);

        // u_N[i]= (u1[i] - u1k) * (u1[i] - u1k) +
        //   (u2[i] - u2k) * (u2[i] - u2k);
    }

    // getminmax(&min,&max,u_N,size);

    // err_D =max;
    err_D /= size;
    (*err) = err_D;
}

/*
 * - Name: getD_Du

 *
*/
void tvcsad_getD(float *xi11, float *xi12, float *xi21, float *xi22, float *u1x, float *u1y,float *u2x, float *u2y,float tau, int size){
    //#pragma omp parallel for

    for (int i = 0; i < size; i++)
    {
        float xi1_N = hypot(xi11[i],xi12[i]);
        float xi2_N = hypot(xi21[i],xi22[i]);

        xi1_N = MAX(1,xi1_N);
        xi2_N = MAX(1,xi2_N);

        xi11[i] = (xi11[i] + tau*u1x[i])/xi1_N;
        xi12[i] = (xi12[i] + tau*u1y[i])/xi1_N;

        xi21[i] = (xi21[i] + tau*u2x[i])/xi2_N;
        xi22[i] = (xi22[i] + tau*u2y[i])/xi2_N;
    }
}


void tvcsad_PD(
        float *I0,           // source image
        float *I1,           // target image
        float *xi11,
        float *xi12,
        float *xi21,
        float *xi22,
        const float lambda,  // weight of the data term
        const float theta,   // weight of the data term
        const float tau,     // time step
        const float tol_OF,  // tol max allowed
        const int   nx,      // image width
        const int   ny,      // image height
        const int   warps,   // number of warpings per scale
        const bool  verbose,  // enable/disable the verbose mode
        float *u1,      // x component of the optical flow
        float *u2      // y component of the optical flow
        )
{

    const float l_t = lambda * theta;
    const int   size = nx * ny;
    const int n_d = DT_NEI;
    const int r = DT_R;
    PosNei *p = new PosNei[size];

    float *u1x    = new float[size];
    float *u1y    = new float[size];
    float *u2x    = new float[size];
    float *u2y    = new float[size];

    float *v1     = new float[size];
    float *v2     = new float[size];

    float *rho_c  = new float[size];
    float *grad   = new float[size];

    float *u1_    = new float[size];
    float *u2_    = new float[size];

    float *u1_tmp = new float[size];
    float *u2_tmp = new float[size];

    float *I1x    = new float[size];
    float *I1y    = new float[size];

    float *I1w    = new float[size];
    float *I1wx   = new float[size];
    float *I1wy   = new float[size];

    //Divergence
    float *div_xi1 = new float[size];
    float *div_xi2 = new float[size];

    float *u_N = new float[size];

    //Five point gradient of the right,left view. (1/12)*[-1 8 0 -8 1]
    centered_gradient(I1, I1x, I1y, nx, ny);
    initialize_pos_nei(nx, ny, n_d, r, p);

    for (int warpings = 0; warpings < warps; warpings++)
    {
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1,  u1, u2, I1w,  nx, ny, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);
        // #pragma omp parallel for
        for (int i = 0; i < size; i++)
        {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1(p + u))| (Warping image)
            grad[i] = hypot(Ix2 + Iy2,0.01);

            for (int j = 0; j < n_d; j++)
            {
                // std::printf("I:%d Iter:%d J:%d I:%d \n",i, j,  p[i].apj[j], p[i].api[j]);
                if (positive(p[i].api[j]) && positive(p[i].apj[j]))
                {
                    // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, p[i].apj[j]*nx + p[i].api[j], p[i].apj[j], p[i].api[j]);
                    assert(p[i].api[j] >= 0);
                    assert(p[i].apj[j] >= 0);
                    assert(p[i].apj[j]*nx + p[i].api[j] < nx*ny);
                    const int pos = p[i].apj[j]*nx + p[i].api[j];

                    p[i].b[j] = (I0[i] - I0[pos] - I1w[i] + I1w[pos] + I1wx[i] * u1[i]
                                 + I1wy[i] * u2[i])/grad[i];
                }
            }
        }

        for (int i = 0; i < nx*ny; i++){
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS)
        {
            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)
            // #pragma omp parallel for
            for (int i = 0; i < size; i++){
                int it = 0;
                for (int j = 0; j< n_d; j++)
                {
                    if (positive(p[i].api[j]) && positive(p[i].apj[j]))
                    {
                        // std::printf("J:%d I:%d\n",p[i].api[j],p[i].apj[j]);
                        p[i].ba[it] = -(p[i].b[j] -  (I1wx[i] * u1[i]
                                                      + I1wy[i] * u2[i])/grad[i]);
                        it++;
                    }
                }
                for (int j = 0; j < (p[i].n+1); j++)
                {
                    p[i].ba[it]= (p[i].n - 2*j)*l_t*grad[i];
                    it++;
                }

                std::sort(p[i].ba.begin(), p[i].ba.begin() + it);
                // v1[i] = u1[i] - l_t*I1wx[i]*p[i].ba[it/2+1]/grad[i];
                // v2[i] = u2[i] - l_t*I1wy[i]*p[i].ba[it/2+1]/grad[i];
                //TODO: Posible error en la minimizacion
                v1[i] = u1[i] - I1wx[i]*p[i].ba[it/2+1]/grad[i];
                v2[i] = u2[i] - I1wy[i]*p[i].ba[it/2+1]/grad[i];
            }
            //Data term

            //Dual variables
            forward_gradient(u1_,u1x,u1y,nx,ny);
            forward_gradient(u2_,u2x,u2y,nx,ny);
            tvcsad_getD(xi11,xi12,xi21,xi22,u1x,u1y,u2x,u2y,tau,size);

            //Primal variables
            divergence(xi11,xi12,div_xi1,nx,ny);
            divergence(xi21,xi22,div_xi2,nx,ny);

            //Almacenamos la iteracion anterior
            for (int i = 0; i < size; i++){
                u1_tmp[i] = u1[i];
                u2_tmp[i] = u2[i];
            }

            tvcsad_getP(u1,u2,v1,v2,div_xi1,div_xi2,u_N,theta,tau,size,&err_D);

            //(aceleration = 1);
            for (int i = 0; i < size; i++){
                u1_[i] = 2*u1[i] - u1_tmp[i];
                u2_[i] = 2*u2[i] - u2_tmp[i];

            }



        }
        if (verbose)
            fprintf(stderr, "Warping: %d,Iter: %d "
                            "Error: %f\n", warpings,n, err_D);
    }

    delete [] p;

    delete [] u1x;
    delete [] u1y;
    delete [] u2x;
    delete [] u2y;

    delete [] v1;
    delete [] v2;

    delete [] rho_c;
    delete [] grad;

    delete [] u1_;
    delete [] u2_;

    delete [] u1_tmp;
    delete [] u2_tmp;

    delete [] I1x;
    delete [] I1y;

    delete [] I1w;
    delete [] I1wx;
    delete [] I1wy;

    delete [] div_xi1;
    delete [] div_xi2;

    delete [] u_N;
    std::printf("Sale del nivel\n");

}

/////////////////////////////////////NLTV-CSAD//////////////////////////////////


void nltvcsad_PD(
        float *I0,           // source image
        float *I1,           // target image
        float *a,            // source image (color)
        int pd,               //Channel numer
        const float lambda,  // weight of the data term
        const float theta,   // weight of the data term
        const float tau,     // time step
        const float tol_OF,  // tol max allowed
        const int   w,      // image width
        const int   h,      // image height
        const int   warps,   // number of warpings per scale
        const bool  verbose, // enable/disable the verbose mode
        float *u1,       // x component of the optical flow
        float *u2      // y component of the optical flow
        )
{

    const int   size = w * h;
    const float l_t = lambda * theta;

    DualVariables_global *p = new DualVariables_global[size];
    DualVariables_global *q = new DualVariables_global[size];
    PosNei *pnei = new PosNei[size];
    float *v1     = new float[size];
    float *v2     = new float[size];
    float *rho_c  = new float[size];
    float *grad   = new float[size];
    float *u1_    = new float[size];
    float *u2_    = new float[size];
    float *u1_tmp = new float[size];
    float *u2_tmp = new float[size];
    float *I1x    = new float[size];
    float *I1y    = new float[size];
    float *I1w    = new float[size];
    float *I1wx   = new float[size];
    float *I1wy   = new float[size];
    float *div_p  = new float[size];
    float *div_q  = new float[size];

    int radius = MAX_BETA;
    int n_d = MAX_DUAL_VAR;

    const int ndt = DT_NEI;
    const int rdt = DT_R;


    //Initialization of the Dual variables.
    initialize_dual_variables(a, pd, w, h, n_d, radius, p, q);
    initialize_pos_nei(w, h, ndt, rdt, pnei);
    centered_gradient(I1, I1x, I1y, w, h);

    for (int warpings = 0; warpings < warps; warpings++)
    {
        // compute the warping of the Right image and its derivatives Ir(x + u1o), Irx (x + u1o) and Iry (x + u2o)
        bicubic_interpolation_warp(I1,  u1, u2, I1w,  w, h, true);
        bicubic_interpolation_warp(I1x, u1, u2, I1wx, w, h, true);
        bicubic_interpolation_warp(I1y, u1, u2, I1wy, w, h, true);
        //#pragma omp parallel for
        for (int i = 0; i < size; i++)
        {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1(p + u))| (Warping image)
            grad[i] = Ix2 + Iy2;
            if (grad[i] > GRAD_IS_ZERO)
            {
                for (int j = 0; j < ndt; j++)
                {
                    // std::printf("I:%d Iter:%d J:%d I:%d \n",i, j,  p[i].apj[j], p[i].api[j]);
                    if (positive(pnei[i].api[j]) && positive(pnei[i].apj[j]))
                    {
                        // std::printf("I:%d Iter:%d Pos: %d J:%d I:%d \n",i, j, p[i].apj[j]*nx + p[i].api[j], p[i].apj[j], p[i].api[j]);
                        assert(pnei[i].api[j] >= 0);
                        assert(pnei[i].apj[j] >= 0);
                        assert(pnei[i].apj[j]*w + pnei[i].api[j] < w*h);
                        const int pos = pnei[i].apj[j]*w + pnei[i].api[j];

                        pnei[i].b[j] = (I0[i] - I0[pos] - I1w[i] + I1w[pos] + I1wx[i] * u1[i]
                                        + I1wy[i] * u2[i])/sqrt(grad[i]);
                    }
                }
            }
        }

        for (int i = 0; i < size; i++)
        {
            u1_[i] = u1[i];
            u2_[i] = u2[i];
        }

        int n = 0;
        float err_D = INFINITY;
        // while (err_D > tol_OF*tol_OF && n < MAX_ITERATIONS)
        while (n < MAX_ITERATIONS)
        {
            n++;
            // estimate the values of the variable (v1, v2)
            //#pragma omp parallel for
            for (int i = 0; i < size; i++){
                v1[i] = u1[i];
                v2[i] = u2[i];
                if (grad[i] > GRAD_IS_ZERO)
                {
                    int it = 0;
                    for (int j = 0; j< ndt; j++)
                    {
                        if (positive(pnei[i].api[j]) && positive(pnei[i].apj[j]))
                        {
                            pnei[i].ba[it] = -(pnei[i].b[j] -  (I1wx[i] * u1[i]
                                                                + I1wy[i] * u2[i])/sqrt(grad[i]));
                            it++;
                        }
                    }
                    for (int j = 0; j < (pnei[i].n+1); j++)
                    {
                        pnei[i].ba[it]= (pnei[i].n - 2*j)*l_t*sqrt(grad[i]);
                        it++;
                    }

                    std::sort(pnei[i].ba.begin(), pnei[i].ba.begin() + it);
                    // v1[i] = u1[i] - l_t*I1wx[i]*pnei[i].ba[it/2+1]/grad[i];
                    // v2[i] = u2[i] - l_t*I1wy[i]*pnei[i].ba[it/2+1]/grad[i];
                    //TODO: Posible error en la minimizacion
                    v1[i] = u1[i] - I1wx[i]*pnei[i].ba[it/2+1]/sqrt(grad[i]);
                    v2[i] = u2[i] - I1wy[i]*pnei[i].ba[it/2+1]/sqrt(grad[i]);
                }
            }
            //Dual variables
            ofnltv_getD(u1_, u2_, size, n_d, tau, p, q);
            //Almacenamos la iteracion anterior
            for (int i = 0; i < size; i++)
            {
                u1_tmp[i] = u1[i];
                u2_tmp[i] = u2[i];
            }

            //Primal variables
            non_local_divergence(p, size, n_d, div_p);
            non_local_divergence(q, size, n_d, div_q);
            ofnltv_getP(v1, v2, div_p, div_q, theta, tau, size, u1, u2, &err_D);
            //(aceleration = 1);
            for (int i = 0; i < size; i++)
            {
                u1_[i] = 2*u1[i] - u1_tmp[i];
                u2_[i] = 2*u2[i] - u2_tmp[i];
            }

        }
        if (verbose)
            std::printf("Warping: %d,Iter: %d Error: %f\n", warpings,n, err_D);
    }

    delete [] p;
    delete [] q;
    delete [] pnei;

    delete [] v1;
    delete [] v2;

    delete [] rho_c;
    delete [] grad;

    delete [] u1_;
    delete [] u2_;

    delete [] u1_tmp;
    delete [] u2_tmp;

    delete [] I1x;
    delete [] I1y;

    delete [] I1w;
    delete [] I1wx;
    delete [] I1wy;

    delete [] div_p;
    delete [] div_q;
}












///////////////////////////////MAIN///////////


/**
 *
 *  Function to read images using the iio library
 *  It always returns an allocated the image.
 *
 */
static float *read_image(const char *filename, int *w, int *h) {
    float *f = iio_read_image_float(filename, w, h);
    if (!f)
        fprintf(stderr, "ERROR: could not read image from file "
                        "\"%s\"\n", filename);
    return f;
}

void rgb2gray(float *in, int w, int h, float *out) {
    int size = w*h;
    for (int i = 0; i < size; i++) {
        out[i] = .299*in[i] + .587*in[size + i] + .114*in[2*size + i];

    }

}


// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
static char *pick_option(int *c, char ***v, char *o, char *d) {
    int argc = *c;
    char **argv = *v;
    int id = d ? 1 : 0;
    for (int i = 0; i < argc - id; i++)
        if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o)) {
            char *r = argv[i+id]+1-id;
            *c -= id+1;
            for (int j = i; j < argc - id; j++)
                (*v)[j] = (*v)[j+id+1];
            return r;
        }
    return d;
}


/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -nprocs      number of threads to use (OpenMP library)
 *   -I0          first image
 *   -I1          second image
 *   -tau         time step in the numerical scheme
 *   -theta       attachment parameter between E_Data and E_Smooth
 *   -nscales     number of scales in the pyramidal structure
 *   -zfactor     downsampling factor for creating the scales
 *   -nwarps      number of warps per scales
 *   -out         name of the output flow field
 *   -verbose     switch on/off messages
 *
 */
int main(int argc, char *argv[]) {

    char *var_reg = pick_option(&argc, &argv, (char *)"m", (char *)"0");
    char *warps_val = pick_option(&argc, &argv, (char *)"w", (char *)"1");

    if (argc != 4) {
        fprintf(stderr, "Usage: %s  ims.txt in_flow.flo  out.flo"
                        "  [-m] val [-w] val"
                //                       0   1  2   3     4    5   6
                "\n", *argv);
        // 4
        return EXIT_FAILURE;
    }
    // O TV-l2 coupled 1 - ||Du + Du'||_{F}
    int val_method = atoi(var_reg);
    int nwarps = atoi(warps_val);


    //read the parameters

    //images
    int i = 1;
    //filename that contains all the images to use
    char* filename_images  = argv[i]; i++;
    char* image_flow_name = argv[i]; i++;
    char* outfile = argv[i]; i++;

    //filename of images
    char *filename_i_1;
    char *filename_i0;
    char *filename_i1;

    //Read txt file of images
    string line;
    ifstream infile;
    int num_files = 0;
    infile.open (filename_images);
    while(getline(infile, line)){

        ++num_files;
        if (num_files == 3){
            filename_i_1  = strdup(line.c_str());
        }else{
            if (num_files == 1){
                filename_i0  = strdup(line.c_str());
            }else{
                if (num_files == 2){
                    filename_i1  = strdup(line.c_str());
                }
            }
        }
    }
    infile.close();

    //read parameters
    float lambda  = (argc > i)? atof(argv[i]): PAR_DEFAULT_LAMBDA; i++;
    float theta   = (argc > i)? atof(argv[i]): PAR_DEFAULT_THETA;  i++;
    float tau     = (argc > i)? atof(argv[i]): PAR_DEFAULT_TAU;    i++;
    float tol_D   = (argc > i)? atof(argv[i]): PAR_DEFAULT_TOL_D;    i++;
    int   nproc    = (argc > i)? atoi(argv[i]): PAR_DEFAULT_NPROC;   i++;
    int   verbose  = (argc > i)? atoi(argv[i]): PAR_DEFAULT_VERBOSE; i++;


    //check parameters
    if (lambda <= 0) {
        lambda = PAR_DEFAULT_LAMBDA;
        if (verbose) fprintf(stderr, "warning: "
                                     "lambda changed to %g\n", lambda);
    }

    if (theta <= 0) {
        tau = PAR_DEFAULT_THETA;
        if (verbose) fprintf(stderr, "warning: "
                                     "theta changed to %g\n", theta);
    }

    if (tau <= 0 || tau > 0.25) {
        tau = PAR_DEFAULT_TAU;
        if (verbose) fprintf(stderr, "warning: "
                                     "tau changed to %g\n", tau);
    }

    if (tol_D <= 0) {
        tol_D = PAR_DEFAULT_TOL_D;
        if (verbose) fprintf(stderr, "warning: "
                                     "tol_D changed to %f\n", tol_D);
    }

    if (nproc < 0) {
        nproc = PAR_DEFAULT_NPROC;
        if (verbose) fprintf(stderr, "warning: "
                                     "nproc changed to %d\n", nproc);
    }



    // read the input images
    // nx es width y ny es height


    // open input images
    int w[4], h[4], pd[4];
    float *i_1;
    if (num_files == 3){
        i_1 = iio_read_image_float_split(filename_i_1, w + 3, h + 3, pd + 3);
    }else{
        i_1 = iio_read_image_float_split(filename_i1, w + 3, h + 3, pd + 3);
    }

    float *i0   = iio_read_image_float_split(filename_i0, w + 0, h + 0, pd + 0);
    float *i1   = iio_read_image_float_split(filename_i1, w + 1, h + 1, pd + 1);
    float *flow = iio_read_image_float_split(image_flow_name, w + 2, h + 2, pd + 2);
    //Ensure that dimensions match
    if (num_files == 3){
        if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
        if (w[0] != w[3] || h[0] != h[3] || pd[0] != pd[3])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
        if (w[1] != w[3] || h[1] != h[3] || pd[1] != pd[3])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
    }else{
        if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
            return fprintf(stderr, "ERROR: input images and flow size mismatch\n");
    }
    //Ensure dimensions match between flow and images
    if (w[0] != w[2] || h[0] != h[2] || pd[2] != 2)
        return fprintf(stderr, "ERROR: input flow field size mismatch\n");




    float *a;
    float *i0n;
    float *i1n;
    float *i_1n;
    float *xi11;
    float *xi12;
    float *xi21;
    float *xi22;

    int size = w[0]*h[0];
    // 0 - TVl2 coupled, otherwise Du
    if (val_method == M_NLTVL1 || val_method == M_NLTVL1_W || val_method == M_NLTVCSAD || val_method == M_NLTVCSAD_W ){
        printf("NL-TVL1 or NLTV-CSAD\n");
        std::printf("W:%d H:%d Pd:%d\n", w[0], h[0], pd[0]);
        a = new float[size*pd[0]];
        image_to_lab(i0, size, a);
    }

    i0n = new float[size];
    i1n = new float[size];

    if (pd[0] != 1){

        rgb2gray(i0, w[0], h[0], i0n);
        rgb2gray(i1, w[0], h[0], i1n);
        rgb2gray(i_1, w[0], h[0], i_1n);

    }else{

        memcpy(i0n, i0, size*sizeof(float));
        memcpy(i1n, i1, size*sizeof(float));
        memcpy(i_1n, i_1, size*sizeof(float));
    }
    image_normalization_3(i0n, i1n, i_1n, i0n, i1n, i_1n, size);
//    normalization(i0n, i1n, i0n, i1n, size);
    gaussian(i0n, w[0], h[0], PRESMOOTHING_SIGMA);
    gaussian(i1n, w[0], h[0], PRESMOOTHING_SIGMA);
    gaussian(i_1n, w[0], h[0], PRESMOOTHING_SIGMA);


    if (verbose)
        fprintf(stderr,"tau=%2.3f tol_D=%2.3f theta=%2.3f"
                       " lambda=%2.3f\n", tau, tol_D, theta, lambda);

    //allocate memory for the flow
    float *u = new float[size*2];
    float *v = u + size;


    //Initialize flow with flow from local faldoi
    for (int i = 0; i < size; i++){
        u[i] = flow[i];
        v[i] = flow[size + i];
    }


    //Initialize dual variables if necessary (TV)
    if (val_method == M_TVL1 || val_method == M_TVL1_W || val_method == M_TVCSAD || val_method == M_TVCSAD_W) {
        xi11 = new float[size];
        xi12 = new float[size];
        xi21 = new float[size];
        xi22 = new float[size];

        for (int i = 0; i < size; i++){
            xi11[i] = 0.0;
            xi12[i] = 0.0;
            xi21[i] = 0.0;
            xi22[i] = 0.0;
        }
    }

    // 0 - TVl2 coupled, otherwise Du
    if (val_method == M_TVL1 || val_method == M_TVL1_W){
        printf("TV-l2 coupled\n");
        tvl2OF(i0n, i1n, u, v, xi11, xi12, xi21, xi22,
               lambda, theta, tau, tol_D, w[0], h[0], nwarps, verbose);
    }else if (val_method == M_NLTVCSAD || val_method == M_NLTVCSAD_W){
        lambda = 0.85;
        theta  = 0.3;
        tau    = 0.1;
        printf("NLTV-CSAD\n");
        nltvcsad_PD(i0n, i1n, a, pd[0], lambda, theta, tau, tol_D,
                w[0], h[0], nwarps, verbose, u, v);
    }else if (val_method == M_NLTVL1 || val_method == M_NLTVL1_W){
        lambda = 2.0;
        theta  = 0.3;
        tau    = 0.1;
        printf("NLTV-L1\n");
        nltvl1_PD(i0n, i1n, a, pd[0], lambda, theta, tau, tol_D,
                w[0], h[0], nwarps, verbose, u, v);
    }else if (val_method == M_TVCSAD || val_method == M_TVCSAD_W){
        lambda = 0.85;
        theta  = 0.3;
        tau    = 0.125;
        printf("TV-CSAD\n");
        tvcsad_PD(i0n, i1n, xi11, xi12, xi21, xi22,
                  lambda, theta, tau, tol_D, w[0], h[0], nwarps, verbose, u,v);

    }else if (val_method == M_TVL1_OCC){
        lambda = 0.25;
        theta = 0.3;
        const float beta = 1;
        const float alpha = 0.01;
        const float tau_u = 0.125;
        const float tau_eta = 0.125;
        const float tau_chi = 0.125;
        tvl2OF_occ(i0n, i1n, i_1n, u, v, xi11, xi12, xi21, xi22,
               lambda, theta, tau_u, tau_eta, tau_chi, beta, alpha, tol_D, w[0], h[0], nwarps, verbose);

    }
    iio_save_image_float_split(outfile, u, w[0], h[0], 2);

    //delete allocated memory
    delete [] u;
    if (val_method == M_TVL1 || val_method == M_TVL1_W || val_method == M_TVCSAD || val_method == M_TVCSAD_W){
        delete [] xi11;
        delete [] xi12;
        delete [] xi21;
        delete [] xi22;
    }else{
        delete [] a;
    }

    delete [] i0n;
    delete [] i1n;

    return EXIT_SUCCESS;
}

#endif//GLOBAL_FALDOI
