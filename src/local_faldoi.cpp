// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// All rights reserved.
#ifndef LOCAL_FALDOI
#define LOCAL_FALDOI

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string.h>
#include <queue>
#include <random>
#include <future>
#include <algorithm>
#include <vector>
#include "energy_structures.h"
#include "aux_energy_model.h"
#include "energy_model.h"
#include "utils_preprocess.h"

extern "C" {
#include "iio.h"
#include "bicubic_interpolation.h"
#include "elap_recsep.h"
}
#include <omp.h>

#include <iostream>
#include <fstream>
#include <string>

#include "utils.h"

using namespace std;

//GLOBAL VARIABLES
// char GLOBAL_TMP_FILE[] = "/tmp/faldoy_XXXXXX"; // template for our file.

#define LOCAL_ITER 3
#define TU_TOL 0.01
#define FB_TOL 2

#define MAX_PATCH 50

#define MAX(x,y) ((x)>(y)?(x):(y))


// typedef  __gnu_pbds::priority_queue<SparseOF, CompareSparseOF,__gnu_pbds::pairing_heap_tag> pq_cand;

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////OUTLIERS FUNCTIONS/////////////////////////////
////////////////////////////////////////////////////////////////////////////////
static void delete_random(
        float tol,
        int w,
        int h,
        float *en_in0,
        float *en_in1
        ) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    int size = w*h;
    int n = 0;

    for (int i = 0; i < size; i++){
        float val = dis(gen);
        if (val > tol){
            en_in0[i] = INFINITY;
            en_in1[i] = INFINITY;
        }else{
            n++;
        }
    }
    std::printf("Too-Chosen: %f\n", (n*1.0)/size);
}

static void rand_local_patch_ini(
        float *u,
        int radius
        ) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-300.0, 300.0);
    int size = (2*radius + 1)*(2*radius + 1);
    for (int i = 0; i < size; i++){
        float val = dis(gen);
        assert(std::isfinite(val));
        u[i] = val;
    }
}

static void zero_local_patch_ini(
        float *u,
        int radius
        ) {
    int size = (2*radius + 1)*(2*radius + 1);
    for (int i = 0; i < size; i++){
        u[i] = 0.0;
    }
}

static float getsample_inf(float *x, int w, int h, int pd, int i, int j, int l) {
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return INFINITY;
    return x[(i + j*w)*pd + l];
}

static int too_uniform(float *a, float tol, int i, int j, int w, int h, int pd){
    float difference = 0;
    int neighborhood[4][2] = {
        {0,1}, {0,-1}, {1, 0}, {-1, 0}};
    for (int l = 0; l < pd; l++){
        float center = getsample_inf(a, w, h, pd, i, j , l);
        for (int k = 0; k < 4; k++){
            int px = i + neighborhood[k][0];
            int py = j + neighborhood[k][1];
            float neighborhood = getsample_inf(a, w, h, pd, px, py , l);
            if (std::isfinite(center) && std::isfinite(neighborhood)){
                float tmp = std::abs(neighborhood-center);
                //std::printf("Tmp: %f, Tol: %f neighborhood: %f Center: %f\n", tmp, difference, neighborhood, center);
                if (difference < tmp){
                    difference = tmp;
                }
            }
        }
    }

    if (difference < tol){
        return 1;
    }
    return 0;

    //return difference < (float) tol;
}

void too_uniform_areas(      
        float *a,
        float *b,
        float *in0,
        int  *trust_in0,
        int w,
        int h,
        float tol
        ){

    float *bw = new float[w*h];
    int size = w*h;
    int n = 0;

    bicubic_interpolation_warp(b, in0, in0 + size, bw, w, h, true);
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++){
            //If both areas present too uniform pixels, we remove the flow.
            if ((too_uniform(a, tol, i, j, w, h, 1) == 1 ) || (too_uniform(bw, tol, i, j, w, h, 1) == 1)){
                trust_in0[j*w + i] = 0;
            }else{
                trust_in0[j*w + i] = 1;
                n++;
            }
        }
    std::printf("Too-Chosen: %f\n", (n*1.0)/size);

    delete [] bw;
}

// Check forward-backward consistency check for of |u(x) + v(x+u(x))| < eps.
// Energy map related to that flows are put to INFINITY.
void fb_consistency_check(
        float *in0,
        float *in1,
        int *trust_in0,
        int w,
        int h,
        float epsilon
        ) {

    float *u1w = new float[w*h];
    float *u2w = new float[w*h];
    int size = w*h;
    int n = 0;

    bicubic_interpolation_warp(in1, in0, in0 + size, u1w, w, h, true);
    bicubic_interpolation_warp(in1 + w*h, in0, in0 + size, u2w, w, h, true);

    for (int i = 0; i < size; i++) {
        float tolerance = hypotf(in0[i] + u1w[i],in0[size + i] + u2w[i]);
        if (tolerance > epsilon){
            //(tolerance > epsilon) means Pixel Occluded
            trust_in0[i] = 0;
        }else{
            trust_in0[i] = 1;
            n++;
        }
    }
    std::printf("FB-Chosen: %f\n", (n*1.0)/size);
    delete [] u2w;
    delete [] u1w;
}

void pruning_method(
        OpticalFlowData *ofGo,
        OpticalFlowData *ofBa,
        float *a, //I0
        float *b, //I1
        int w, //width image
        int h, //height image
        float *tol, //tolerance too_uniform and f-b
        int *method, // if method[i]!=0, then
        int *trust_Go, //energy map of u
        float *go, //of to t, t+1
        int *trust_Ba, //energy map of v
        float *ba //of to t+1, t
        ) {
    int *go_fb_check   = new int[w*h];
    int *go_cons_check = new int[w*h];
    int *ba_fb_check   = new int[w*h];
    int *ba_cons_check = new int[w*h];

    for (int i = 0; i < w*h; i++){
        //0 - Invalid pixel 1 - Trustable pixel.
        trust_Go[i] = 1;
        trust_Ba[i] = 1;
    }


    //FB - consistency check
    if (method[0] == 1) {
        std::printf("FB-Consistency: %f\n", tol[0]);
        fb_consistency_check(go, ba, go_fb_check, w, h, tol[0]);
        fb_consistency_check(ba, go, ba_fb_check, w, h, tol[0]);
    }
    //Too-uniform consistency check
    if (method[1] == 1){
        std::printf("Too Uniform -Consistency: %f\n", tol[1]);
        too_uniform_areas(a, b, go, go_cons_check, w, h, tol[1]);
        too_uniform_areas(b, a, ba, ba_cons_check, w, h, tol[1]);
    }
    for (int i = 0; i < w*h; i++){
        if (method[0] == 1) {
            //FB-Consistency
            if (go_fb_check[i] == 0) {
                trust_Go[i] = 0;
            }
            if (ba_fb_check[i] == 0) {
                trust_Ba[i] = 0;
            }
        }
        //Too uniform -Consistency
        if (method[1] == 1) {
            if (go_cons_check[i] == 0) {
                trust_Go[i] = 0;
            }
            if (ba_cons_check[i] == 0) {
                trust_Ba[i] = 0;
            }
        }
    }

    delete [] go_fb_check;
    delete [] go_cons_check;
    delete [] ba_fb_check;
    delete [] ba_cons_check;
}

void delete_not_trustable_candidates(
        SpecificOFStuff *ofS,
        OpticalFlowData *ofD,
        float *in,
        float *ene_val
        ) {
    int *mask = ofD->trust_points;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    int w = ofD->params.w;
    int h = ofD->params.h;
    int n = 0;
    for (int i = 0; i < w*h; i++) {
        if (mask[i] == 0) {
            //printf("%f\n", ene_val[i]);
            if (ene_val[i] == 0.0){
                n++;
            }
            in[i]       = NAN;
            in[i + w*h] = NAN;
            u1[i]       = NAN;
            u2[i]       = NAN;
            ene_val[i]  = INFINITY;
        }
    }
    printf("Total_seeds: %d\n", n);
}


void compare_seeds(float *old_s, float *new_s, int w, int h) {
    int size = w*h;
    for (int i = 0; i < size; i++) {
        //Dejamos las semillas iniciales que pasan el pruning
        if (std::isfinite(old_s[i]) && std::isfinite(old_s[size + i]) &&
                !std::isfinite(new_s[i]) && !std::isfinite(new_s[size + i])) {
            old_s[i] = NAN;
            old_s[size + i] = NAN;
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
//////////////////LOCAL INITIALIZATION//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//Constant Interpolation
inline void interpolate_constant(
        OpticalFlowData *ofD,
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej, // end row
        float *v) {
    const int nx = ofD->params.w;
    int *mask = ofD->fixed_points;
    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++) {
            const int i = l*nx + k;
            //1 fixed - 0 not
            //if (mask[i] == 0){
            ofD->u1[i] = v[0];
            ofD->u2[i] = v[1];
            //}
        }
}



//Poisson Interpolation
void interpolate_poisson(
        OpticalFlowData *ofD,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej // end row
        ) {
    int w = ei - ii;
    int h = ej - ij;
    int *mask = ofD->fixed_points;
    int wR = ofD->params.w;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    float buf_in[2*MAX_PATCH*MAX_PATCH];
    float buf_out[2*MAX_PATCH*MAX_PATCH];
    assert(w * h < MAX_PATCH * MAX_PATCH);
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int x = i + ii;
            int y = j + ij;
            int xy = y * wR + x;
            //1 fixed - 0 not
            if (mask[xy] == 1) {
                buf_in[ j*w + i       ] = u1[ xy ];
                buf_in[ j*w + i + w*h ] = u2[ xy ];
            } else  { //not fixed
                buf_in[ j*w + i       ] = NAN;
                buf_in[ j*w + i + w*h ] = NAN;
            }
        }
    elap_recursive_separable(buf_out, buf_in, w, h, 2, 0.4, 3, 7);
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int x = i + ii;
            int y = j + ij;
            int xy = y * wR + x;
            u1[ xy ] = buf_out[ j*w + i       ];
            u2[ xy ] = buf_out[ j*w + i + w*h ];
        }
}

void nltv_regularization(
        OpticalFlowData *ofD,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej, // end row
        float tau, //timestep to regularize
        int niter //Iteration number
        ) {

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    int *mask = ofD->fixed_points;
    int n_d = NL_DUAL_VAR;

    //Columns and Rows
    const int w = ofD->params.w;

    BilateralWeight *p = ofD->weight;

    //Get the correct wt to force than the sum will be 1
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            float wt_tmp = 0.0;
            const int i = l*w + k;
            for (int j = 0; j < n_d; j++) {
                const int api = p[i].api[j];
                const int apj = p[i].apj[j];
                //The position should be the same
                const int ap = validate_ap_patch(ii, ij, ei, ej, w, api, apj);
                if (ap == 0) {
                    const float wp = p[i].wp[j];
                    assert(wp >= 0);
                    wt_tmp +=wp;
                }
            }
            p[i].wt = wt_tmp;
        }
    }
    //We perform the number of iterations
    int n = 0;
    while (n < niter) {
        n++;
        for (int l = ij; l < ej; l++){
            for (int k = ii; k < ei; k++){
                const int i = l*w + k;

                float g1 = 0.0;
                float g2 = 0.0;
                for (int j = 0; j < n_d; j++) {
                    const int api = p[i].api[j];
                    const int apj = p[i].apj[j];
                    //The position should be the same
                    const int ap = validate_ap_patch(ii, ij, ei, ej, w, api, apj);

                    if (ap == 0) {
                        const float wp = p[i].wp[j];
                        assert(wp >= 0);
                        g1 += (u1[i] - u1[apj*w + api])*wp;
                        g2 += (u2[i] - u2[apj*w + api])*wp;
                    }
                }
                g1 /= p[i].wt;
                g2 /= p[i].wt;
                //If the value is not fixed
                if (mask[i] == 0){
                    const float u1k = u1[i];
                    const float u2k = u2[i];
                    u1[i] = u1k  -tau*g1;
                    u2[i] = u2k  -tau*g2;
                }
            }
        }
    }

}


void bilateral_filter_regularization(
        OpticalFlowData *ofD,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej, // end row
        int niter //Iteration number
        ) {

    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    int *mask = ofD->fixed_points;
    int n_d = NL_DUAL_VAR;

    //Columns and Rows
    const int w = ofD->params.w;

    BilateralWeight *p = ofD->weight;

    //Get the correct wt to force than the sum will be 1
    for (int l = ij; l < ej; l++){
        for (int k = ii; k < ei; k++){
            float wt_tmp = 0.0;
            const int i = l*w + k;
            for (int j = 0; j < n_d; j++) {
                const int api = p[i].api[j];
                const int apj = p[i].apj[j];
                //The position should be the same
                const int ap = validate_ap_patch(ii, ij, ei, ej, w, api, apj);
                if (ap == 0) {
                    const float wp = p[i].wp[j];
                    assert(wp >= 0);
                    wt_tmp +=wp;
                }
            }
            p[i].wt = wt_tmp;
            assert(std::isfinite(wt_tmp));
        }
    }

    //We perform the number of iterations
    int n = 0;
    while (n < niter) {
        n++;
        for (int l = ij; l < ej; l++){
            for (int k = ii; k < ei; k++){
                const int i = l*w + k;

                float g1 = 0.0;
                float g2 = 0.0;
                for (int j = 0; j < n_d; j++) {
                    const int api = p[i].api[j];
                    const int apj = p[i].apj[j];
                    //The position should be the same
                    const int ap = validate_ap_patch(ii, ij, ei, ej, w, api, apj);

                    if (ap == 0) {
                        const float wp = p[i].wp[j];
                        assert(wp >= 0);
                        assert(std::isfinite(g1));
                        assert(std::isfinite(g2));
                        g1 += u1[apj*w + api]*wp;
                        g2 += u2[apj*w + api]*wp;
                    }
                }
                g1 /= p[i].wt;
                g2 /= p[i].wt;
                //If the value is not fixed
                if (mask[i]==0){
                    u1[i] = g1;
                    u2[i] = g2;
                }
            }
        }
    }

}

//Poisson Interpolation + nltv_regularization
void interpolate_poisson_nltv(
        OpticalFlowData *ofD,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej // end row
        ) {
    int w = ei - ii;
    int h = ej - ij;
    int *mask = ofD->fixed_points;
    int wR = ofD->params.w;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    float buf_in[2*MAX_PATCH*MAX_PATCH];
    float buf_out[2*MAX_PATCH*MAX_PATCH];
    assert(w * h < MAX_PATCH * MAX_PATCH);
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int x = i + ii;
            int y = j + ij;
            int xy = y * wR + x;
            //1 fixed - 0 not
            if (mask[xy] == 1) {
                buf_in[ j*w + i       ] = u1[ xy ];
                buf_in[ j*w + i + w*h ] = u2[ xy ];
            } else  { //not fixed
                buf_in[ j*w + i       ] = NAN;
                buf_in[ j*w + i + w*h ] = NAN;
            }
        }
    elap_recursive_separable(buf_out, buf_in, w, h, 2, 0.4, 3, 7);
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int x = i + ii;
            int y = j + ij;
            int xy = y * wR + x;
            u1[ xy ] = buf_out[ j*w + i       ];
            u2[ xy ] = buf_out[ j*w + i + w*h ];
        }
    //Here we use the inizialization of the previous poisson
    float tau = 0.25;
    int niter = 20;
    nltv_regularization(ofD, ii, ij, ei, ej, tau, niter);
}

//Insert 8-connected candidates into the priority queue with their energies.
void insert_candidates(
        pq_cand *queue,
        float *ene_val,
        OpticalFlowData *ofD,
        const int i,
        const int j,
        const float ener_N
        ) {

    int n_neigh = 4;
    int neighborhood[8][2] = {
        {0, 1}, {0, -1}, {1, 0}, {-1, 0},
        {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

    const int w = ofD->params.w;
    const int h = ofD->params.h;
    const float *sal = ofD->saliency;

    for (int k = 0; k < n_neigh; k++){
        int px = i + neighborhood[k][0];
        int py = j + neighborhood[k][1];

        if (px >= 0 && px < w && py >=0 && py < h){
            float new_ener = ener_N * sal[py*w + px];

            //std::printf("Ener_N: %f  Sim: %f \n", ener_N, ene_val[py*w + px]);
            if (!ofD->fixed_points[py*w + px] &&  new_ener < ene_val[py*w + px]){

                ene_val[py*w + px] = ener_N;
                SparseOF element;
                element.i = px; // column
                element.j = py; // row
                element.u = ofD->u1[py*w + px];
                element.v = ofD->u2[py*w + px];
                element.sim_node = new_ener;
                element.occluded = ofD->chi[py*w + px];
                queue->push(element);
            }
        }
    }

}


inline PatchIndexes get_index_patch(
        const int wr,
        const int w,
        const int h,
        const int i,
        const int j,
        const int factor
        ) {
    PatchIndexes index;
    //Points to begin and end. End is the previous value
    index.i = i;
    index.j = j;
    index.ii = ((i - factor * wr) < 0)? 0 : (i - factor * wr);
    index.ij = ((j - factor * wr) < 0)? 0 : (j - factor * wr);
    index.ei = ((i + 1 + factor * wr) > w)? w : (i + 1 + factor * wr);
    index.ej = ((j + 1 + factor * wr) > h)? h : (j + 1 + factor * wr);
    return index;

}

//TODO: Esto esta fatal. Si finalmenente funciona lo de los pesos arreglarlo para
//que faldoi sea independiente y este dentro de energy_model.cpp
inline void get_relative_index_weight(
        int *iiw, // initial column
        int *ijw, // initial row
        const int wr,
        const int i,
        const int j
        ) {

    (*iiw) = (((i -  wr)< 0)? -(i - wr) : 0);
    (*ijw) = (((j -  wr)< 0)? -(j - wr) : 0);
    assert(*iiw >= 0);
    assert(*ijw >= 0);
}
static void get_index_weight(
        int method,
        SpecificOFStuff *ofS,
        const int wr,
        int i,
        int j){

    int iiw, ijw;
    if (method == M_TVL1_W || method == M_NLTVCSAD_W || method == M_NLTVL1_W || method == M_TVCSAD_W) {
        get_relative_index_weight(&iiw, &ijw, wr, i, j);
    }
    switch (method) {
    case M_TVL1_W:
        ofS->tvl2w.iiw = iiw;
        ofS->tvl2w.ijw = ijw;
        break;
    case M_NLTVCSAD_W:
        ofS->nltvcsadw.iiw = iiw;
        ofS->nltvcsadw.ijw = ijw;
        break;
    case M_NLTVL1_W:
        ofS->nltvl1w.iiw = iiw;
        ofS->nltvl1w.ijw = ijw;
        break;
    case M_TVCSAD_W:
        ofS->tvcsadw.iiw = iiw;
        ofS->tvcsadw.ijw = ijw;
        break;
    default:
        break;
    }
}



//Copy over ofD->u1 and ofD->u2 the presented values in out.
inline void copy_fixed_coordinates(
        OpticalFlowData *ofD,
        float *out,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej  // end row
        ) {
    const int w = ofD->params.w;
    const int h = ofD->params.h;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    int *fixed = ofD->fixed_points;

    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++){
            //Copy only fixed values from the patch
            const int i = l*w + k;
            if (fixed[i] == 1){
                u1[i] = (float) out[i];
                u2[i] = (float) out[w*h + i];
                assert(std::isfinite(u1[i]));
                assert(std::isfinite(u2[i]));
            }
        }
}

static inline void update_fixed_coordinates(
        float *out,
        OpticalFlowData *ofD,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej  // end row
        ) {
    const int w = ofD->params.w;
    const int h = ofD->params.h;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    int *fixed = ofD->fixed_points;

    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++){
            //Copy only fixed values from the patch
            const int i = l*w + k;
            if (fixed[i] == 1){
                out[i]       = u1[i];
                out[w*h + i] = u2[i];
                assert(std::isfinite(out[i]));
                assert(std::isfinite(out[w*h + i]));
            }else{
                out[i] = NAN;
                out[w*h + i] = NAN;
            }
        }
}


//Poisson Interpolation
void copy_ini_patch(
        OpticalFlowData *ofD,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej // end row
        ) {
    int w = ei - ii;
    int h = ej - ij;
    int wR = ofD->params.w;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    float *u1_ini = ofD->u1_ini;
    float *u2_ini = ofD->u2_ini;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++){
            int x = i + ii;
            int y = j + ij;
            int xy = y * wR + x;
            assert(std::isfinite(u1_ini[j*w + i]));
            assert(std::isfinite(u2_ini[j*w + i]));
            assert(xy >= 0);
            assert(j*w + i >=0);
            u1[ xy ] = u1_ini[ j*w + i];
            u2[ xy ] = u2_ini[ j*w + i];
        }
}

//Poisson Interpolation
void copy_mix_patch(
        OpticalFlowData *ofD,
        float *out,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej // end row
        ) {
    int w = ei - ii;
    int h = ej - ij;
    int *fixed = ofD->fixed_points;
    int *trust = ofD->trust_points;
    int wR = ofD->params.w;
    float *u1 = ofD->u1;
    float *u2 = ofD->u2;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int x = i + ii;
            int y = j + ij;
            int xy = y * wR + x;
            assert(xy >= 0);
            if (fixed[xy] == 1) {
                u1[ xy ] = out[ xy ];
                u2[ xy ] = out[ xy ];
            }else if ((trust[xy] == 0 ) && (fixed[xy] == 0)){
                u1[ xy ] = NAN;
                u2[ xy ] = NAN;
            }
        }
}


//Check if there is at least one pixel that hasn't survived to the prunning. 
int check_trustable_patch(
        OpticalFlowData *ofD,
        int ii, // initial column
        int ij, // initial row
        int ei, // end column
        int ej // end row
        ) {

    const int w = ofD->params.w;

    int *fixed = ofD->trust_points;

    for (int l = ij; l < ej; l++)
        for (int k = ii; k < ei; k++){
            //Return 0 if it detects that at least one point it is not fixed
            const int i = l*w + k;
            if (fixed[i] == 0){
                //If the pixel it is not trustable.
                return 0;
            }
        }

    return 1;
}


static void add_neighbors(        
        const float *i0,
        const float *i1,
        const float *i_1,
        float *ene_val,
        OpticalFlowData *ofD,
        SpecificOFStuff *ofS,
        pq_cand *queue,
        const int i,
        const int j,
        const int mode,
        float *out,
        float *out_occ
        ) {

    const int w  = ofD->params.w;
    const int h  = ofD->params.h;
    const int wr = ofD->params.w_radio;
    float ener_N;


    const PatchIndexes index = get_index_patch(wr, w, h, i, j, 1);
    //    fprintf()
    int method = ofD->params.val_method;

    //TODO: Arreglar los de los pesos
    get_index_weight(method, ofS, wr, i, j);
    //////////////////////////////////////////
    //
    // FIRST STEP, ADD "POISSON" CANDIDATES
    //
    //////////////////////////////////////////
    //Poisson Interpolation (4wr x 4wr + 1)
    if (mode == 0) {
        //it > 0. Interpolate over the survivors of the pruning.
        copy_fixed_coordinates(ofD, out, index.ii, index.ij, index.ei, index.ej);
        interpolate_poisson(ofD, index.ii, index.ij, index.ei, index.ej);

    }else if (check_trustable_patch(ofD, index.ii, index.ij, index.ei, index.ej) == 0) {

        copy_fixed_coordinates(ofD, out, index.ii, index.ij, index.ei, index.ej);
        interpolate_poisson(ofD, index.ii, index.ij, index.ei, index.ej);
    }

    // Optical flow method on patch (2*wr x 2wr + 1)
    of_estimation(ofS, ofD, &ener_N, i0, i1, i_1, index);

    // update_fixed_coordinates(out, ofD, ii, ij, ei, ej);
    insert_candidates(queue, ene_val, ofD, i, j, (float) ener_N);

    //TODO:It is a strange step, if the energy over the patch is lower thant the
    //stored energy, we put the new one, if it's not, we put the old one.
    if (ene_val[j*w + i] > ener_N) {
        out[      j*w + i] = ofD->u1[j*w + i];
        out[w*h + j*w + i] = ofD->u2[j*w + i];
        ene_val[  j*w + i] = ener_N;
        out_occ[  j*w + i] = ofD->chi[j*w + i];

    }
}



int insert_initial_seeds(
        const float *i0,
        const float *i1,
        const float *i_1,
        float *in,
        pq_cand *queue,
        OpticalFlowData *ofD,
        SpecificOFStuff *ofS,
        int mode,
        float *ene_val,
        float *out,
        float *out_occ
        ) {
    const int w = ofD->params.w;
    const int h = ofD->params.h;
    const int wr = ofD->params.w_radio;
    int nfixed = 0;

    //Set to the initial conditions all the stuff
    for (int i = 0; i < w*h; i++){
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out[i] = NAN;
        out[w*h + i] = NAN;
        out_occ[i] = 0;
    }

    ofD->params.w_radio = 1;
    //Fix the initial seeds.
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++){

            //Indicates the initial seed in the similarity map
            if (std::isfinite(in[j*w +i]) && std::isfinite(in[w*h + j*w + i])){

                out[j*w + i] = in[j*w + i];
                out[w*h + j*w + i] = in[w*h + j*w +i];
                ofD->fixed_points[j*w + i] = 1;

                // add_neigbors 0 means that during the propagation interpolates the patch
                // based on the energy.
                add_neighbors(i0, i1, i_1, ene_val, ofD, ofS, queue, i, j, 0, out, out_occ);

                //These values may have been modified in the previous function
                out[j*w + i] = in[j*w + i];
                out[w*h + j*w + i] = in[w*h + j*w +i];
                ofD->fixed_points[j*w + i] = 1;
                ene_val[j*w + i] = 0.0;
            }
        }
    ofD->params.w_radio = wr;

    return nfixed;
}


//  Insert each pixel into the queue as possible candidate. Its related energy comes
// from the energy store at the moment that the pixel was fixed.
void insert_potential_candidates(
        const float *i0,
        const float *i1,
        const float *in,
        SpecificOFStuff *ofS,
        OpticalFlowData *ofD,
        pq_cand *queue,
        float *ene_val,
        float *out,
        float *out_occ
        ){
    //Note: in and out are the same pointer
    int w = ofD->params.w;
    int h = ofD->params.h;

    //Fixed the initial seeds.
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            //Indicates the initial seed in the similarity map
            if (std::isfinite(in[j*w + i]) && std::isfinite(in[w*h + j*w + i])) {

                SparseOF element;
                element.i = i; // column
                element.j = j; // row
                element.u = in[j*w +i];
                element.v = in[w*h + j*w +i];
                //Obs: Notice that ene_val contains (en)*saliency
                element.sim_node = ene_val[j*w + i];
                element.occluded = out_occ[j*w + i];
                assert(std::isfinite(ene_val[j*w + i]));
                queue->push(element);
            }
        }

    //Set to the initial conditions all the stuff
    for (int i = 0; i < w*h; i++) {
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out[i] = NAN;
        out[w*h + i] = NAN;
        out_occ[i] = 0;
    }
}

static void update_energy_map(
        const float *i0,
        const float *i1,
        const float *i_1,
        float *ene_val,
        SpecificOFStuff *ofS,
        OpticalFlowData *ofD,
        const float *out
        ) {

    const int w  = ofD->params.w;
    const int h  = ofD->params.h;
    const int wr = ofD->params.w_radio;
    const float *sal = ofD->saliency;
    float ener_N;

    //Copy the last version of the optical flow
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++){
            ofD->u1[j*w + i] = out[j*w + i];
            ofD->u2[j*w + i] = out[w*h + j*w + i];
        }

    //Update the energy map measuring the energy for each pixel.
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {

            //TODO:Arreglar los de los pesos
            int method = ofD->params.val_method;
            get_index_weight(method, ofS, wr, i, j);
            //////////////////////////////////////////
            //
            // FIRST STEP, ADD "POISSON" CANDIDATES
            //
            //////////////////////////////////////////
            PatchIndexes index = get_index_patch(wr, w, h, i, j, 1);
            eval_functional(ofS, ofD, &ener_N, i0, i1, i_1, index);
            ene_val[j*w + i] = ener_N * sal[j*w + i];
        }
}


//Initialize the data to prepare everything for the region growing
void prepare_data_for_growing(
        OpticalFlowData *ofD,
        SpecificOFStuff *ofS,
        float *ene_val,
        float *out
        ) {
    int w = ofD->params.w;
    int h = ofD->params.h;

    //Set to the initial conditions all the stuff
    for (int i = 0; i < w*h; i++) {
        ofD->fixed_points[i] = 0;
        ene_val[i] = INFINITY;
        out[i] = NAN;
        out[w*h + i] = NAN;
    }

}

void local_growing(        
        const float *i0,
        const float *i1,
        const float *i_1,
        pq_cand *queue,
        SpecificOFStuff *ofS,
        OpticalFlowData *ofD,
        int tm,
        const int nfixed,
        float *ene_val,
        float *out,
        float *out_occ
        ) {

    const int w = ofD->params.w;
    const int h = ofD->params.h;
    std::printf("Queue size at start = %d\n", (int)queue->size());
    while (! queue->empty()) {

        SparseOF element = queue->top();
        int i = element.i;
        int j = element.j;
        //While the queue is not empty, take an element to process
        queue->pop();

        if (!ofD->fixed_points[j*w + i]){
            assert(std::isfinite(element.sim_node));
            float u = element.u;
            float v = element.v;
            float energy = element.sim_node;

            if (!std::isfinite(u)){
                std::printf("U1 = %f\n", u);
            }
            if (!std::isfinite(v)){
                std::printf("U2 = %f\n", v);
            }

            ofD->fixed_points[j*w + i] = 1;

            //printf("%d\n", val);
            out[j*w + i] = u;
            out[w*h + j*w + i] = v;
            ene_val[j*w + i] = energy;
            out_occ[j*w + i] = u;
            // //TODO: Lo copiamos para que esos valores influyan en la minimizacion.
            // //MIRAR
            // ofD->u1[j*w + i] = u;
            // ofD->u2[j*w + i] = v;

            //tm stores if we made interpolation or not.
            add_neighbors(i0, i1, i_1, ene_val, ofD, ofS, queue, i, j, tm, out, out_occ);
        }
    }
}




void match_growing_variational(
        float *go,
        float *ba,
        float *i0,
        float *i1,
        float *i_1,
        float *i2,
        float *sal_go,
        float *sal_ba,
        Parameters params,
        float *ene_val,
        float *out_flow,
        float *out_occ
        ){
    int w = params.w;
    int h = params.h;
    std::printf("Initializing stuff\n");

    //Initialize all the stuff for optical flow computation
    //Optical flow t, t+1
    OpticalFlowData ofGo = init_Optical_Flow_Data(sal_go, params);
    float *oft0 = new float[w*h*2];
    float *ene_Go = new float[w*h];
    float *occ_Go = new float[w*h];


    //Optical flow t+1, t
    OpticalFlowData ofBa = init_Optical_Flow_Data(sal_ba, params);
    float *oft1 = new float[w*h*2];
    float *ene_Ba = new float[w*h];
    float *occ_Ba = new float[w*h];

    //Create queues
    int nfixed_go = 0;
    int nfixed_ba = 0;
    pq_cand queueGo;
    pq_cand queueBa;


    //Initialize all the auxiliar data.
    SpecificOFStuff stuffGo;
    SpecificOFStuff stuffBa;
    initialize_auxiliar_stuff(stuffGo, ofGo);
    initialize_auxiliar_stuff(stuffBa, ofBa);


    //i0n, i1n, i_1n, i2n are a gray and smooth version of i0, i1, i_1, i2
    float *i0n = nullptr;
    float *i1n = nullptr;
    float *i_1n = nullptr;
    float *i2n = nullptr;

    //Prepare data based on the functional chosen (energy_model.cpp)
    prepare_stuff(&stuffGo, &ofGo, &stuffBa, &ofBa, i0, i1, i_1, i2, params.pd, &i0n, &i1n, &i_1n, &i2n);
    std::printf("Finished initializing stuff\n");


    ////FIXED POINTS////
    //Insert initial seeds to queues
    std::printf("Inserting initial seeds\n");
    //#pragma omp parallel num_threads(2)
    {
        //#pragma omp sections
        {
            //#pragma omp section
//            nfixed_go =  insert_initial_seeds(i0n, i1n, i_1n, go, &queueGo, &ofGo, &stuffGo, 0, ene_Go, oft0, occ_Go);
            auto future_nfixed_go = std::async(std::launch::async,
                                               [&] { return insert_initial_seeds(i0n, i1n, i_1n, go, &queueGo, &ofGo, &stuffGo, 0, ene_Go, oft0, occ_Go); });
            //#pragma omp section
            nfixed_ba = insert_initial_seeds(i1n, i0n, i2n, ba, &queueBa, &ofBa, &stuffBa, 0, ene_Ba, oft1, occ_Ba);
            nfixed_go = future_nfixed_go.get();
        } /// End of sections
    } /// End of parallel section

    std::printf("Finished inserting initial seeds\n");

    const int iter = LOCAL_ITER;
    //Variables for pruning
    float tol[2] = {FB_TOL, TU_TOL};
    int   p[2] = {1 , 0};


    for (int i = 0; i < iter; i++){
        std::printf("Iteration: %d\n", i);

        //#pragma omp parallel num_threads(2)
        {
            //#pragma omp sections
            {
                //#pragma omp section
                //Estimate local minimization I0-I1
                auto growing_fwd = std::async(std::launch::async,
                                              [&] { local_growing(i0n, i1n, i_1n, &queueGo, &stuffGo, &ofGo, i, nfixed_go, ene_Go, oft0, occ_Go); });
                //#pragma omp section
                //Estimate local minimzation I1-I0
                local_growing(i1n, i0n, i2n, &queueBa, &stuffBa, &ofBa, i, nfixed_ba, ene_Ba, oft1, occ_Ba);
                growing_fwd.get();
            } /// End of sections
        } /// End of parallel section

        //Pruning method
        pruning_method(&ofGo, &ofBa, i0n, i1n, w, h, tol, p,
                       ofGo.trust_points, oft0, ofBa.trust_points, oft1);

        //Delete not trustable candidates based on the previous pruning
        delete_not_trustable_candidates(&stuffGo, &ofGo, oft0, ene_Go);
        delete_not_trustable_candidates(&stuffBa, &ofBa, oft1, ene_Ba);

        //Insert each pixel into the queue as possible candidate
        insert_potential_candidates(i0n, i1n, oft0, &stuffGo, &ofGo, &queueGo, ene_Go, oft0, occ_Go);
        insert_potential_candidates(i1n, i0n, oft1, &stuffBa, &ofBa, &queueBa, ene_Ba, oft1, occ_Ba);

        prepare_data_for_growing(&ofGo, &stuffGo, ene_Go, oft0);
        prepare_data_for_growing(&ofBa, &stuffBa, ene_Ba, oft1);

    }
    std::printf("Last growing\n");
    local_growing(i0n, i1n, i_1n, &queueGo, &stuffGo, &ofGo, 10, nfixed_go, ene_Go, oft0, occ_Go);

    //Copy the result t, t+1 as output.
    memcpy(out_flow, oft0, sizeof(float)*w*h*2);
    memcpy(ene_val, ene_Go, sizeof(float)*w*h);
    memcpy(out_occ, occ_Go, sizeof(float)*w*h);


    free_auxiliar_stuff(&stuffGo, &ofGo);
    free_auxiliar_stuff(&stuffBa, &ofBa);

    delete [] i1n;
    delete [] i2n;
    delete [] i0n;
    delete [] i_1n;

    delete [] ofGo.u1;
    delete [] ofBa.u1;

    delete [] ofGo.fixed_points;
    delete [] ofBa.fixed_points;

    delete [] ofGo.trust_points;
    delete [] ofBa.trust_points;

    delete [] ofGo.weight;
    delete [] ofBa.weight;

    delete [] ofGo.u1_ini;
    delete [] ofGo.u2_ini;
    delete [] ofBa.u1_ini;
    delete [] ofBa.u2_ini;


    delete [] oft0;
    delete [] oft1;

    delete [] ene_Go;
    delete [] ene_Ba;

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////MAIN/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void copy_additional_seeds(float *out, float *in, int w, int h){

    int size = w*h;

    for (int j = 0; j < h; j++){
        for (int i = 0; i < w; i++){
            //Copy non-NaN values and if there is colision keep the old value.
            if (std::isfinite(in[j*w +i]) && std::isfinite(in[size + j*w +i])
                    && !std::isfinite(out[j*w +i]) && !std::isfinite(out[size + j*w +i]) ){
                out[j*w + i] = in[j*w + i];
                out[w*h + j*w + i] = in[w*h + j*w + i];
            }
        }
    }
}



//
//Create an artificial optical flow. 
// Input: u(x)
// Output: v(x) = -u(x + u(x)).
static void reverse_optical_flow(float *in, int w, int h, float *out){
    int size = w*h;

    //Set the sparse output to NaN
    for (int i = 0; i < w*h; i++)
    {
        out[i] = NAN;
        out[w*h + i] = NAN;
    }

    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            //If a sparse optical flow.
            // NaN - not flow
            // Real valuables - flow.
            if (std::isfinite(in[j*w +i]) && std::isfinite(in[size + j*w +i])){
                int x = std::floor(i + in[j*w + i]);
                int y = std::floor(j + in[size + j*w + i]);
                //Check that is inside of the image domain.

                if (x >= 0 && x < w && y >=0 && y < h){
                    float val = 1;

                    //We check if there are two different flow for the same position
                    if (std::isfinite(out[y*w +x]) && std::isfinite(out[size + y*w +x])){
                        float pre,now;
                        pre = hypotf(out[y*w +x],out[size + y*w +x]);
                        now = hypotf(-in[y*w +x],-in[size + y*w +x]);
                        fprintf(stderr, "COLLISION: At least two flow for the same pixel\n");
                        //If there is a colision, we put the optical flow with higher norm.
                        if (now < pre){
                            val = 0;
                        }
                    }

                    if (val == 1){
                        out[y*w + x] = -in[j*w + i];
                        out[size + y*w + x] = -in[size + j*w + i];
                    }
                }else{
                    fprintf(stderr, "OUT: Position outside of the image domain \n");
                }
            }
        }

}

/**
 *
 *  Function to read images using the iio library
 *  It always returns an allocated the image.
 *
 */
static float *read_image(const char *filename, int *w, int *h){

    float *f = iio_read_image_float(filename, w, h);
    if (!f)
        fprintf(stderr, "ERROR: could not read image from file "
                        "\"%s\"\n", filename);
    return f;
}


static bool pick_option(std::vector<std::string>& args, const std::string& option){

    auto it = std::find(args.begin(), args.end(), "-" + option);

    bool found = it != args.end();
    if (found)
        args.erase(it);

    return found;
}

static std::string pick_option(std::vector<std::string>& args, const std::string& option, const std::string& default_value) {
    auto arg = "-" + option;

    for (auto it = args.begin(); it != args.end(); it++) {
        if (*it == arg) {
            auto next = it + 1;
            if (next == args.end())
                continue;
            auto result = *next;
            args.erase(it, it + 2);
            return result;
        }
    }
    return default_value;
}




//Main function that expands sparse flow
int main(int argc, char* argv[]){


    using namespace std::chrono;

    system_clock::time_point today = system_clock::now();
    time_t tt;

    tt = system_clock::to_time_t ( today );
    std::cerr << "today is: " << ctime(&tt);


    // process input
    std::vector<std::string> args(argv, argv + argc);
    auto windows_ratio = pick_option(args, "wr", "5"); // Warpings
    auto var_reg       = pick_option(args, "m",  "8"); // Methods
    auto file_params    = pick_option(args, "p",  ""); // File of parameters

    if (args.size() != 7 && args.size() != 9) {
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff occlusions.png"
                //                          0        1     2       3       4       5         6
                " [-m method_id] [-wr windows_radio]\n", args.size(), args[0].c_str());
        fprintf(stderr, "usage %lu :\n\t%s ims.txt in0.flo in1.flo out.flo sim_map.tiff occlusions.png sal0.tiff sal1.tiff"
                //                          0      1     2       3       4       5         6
                " [-m method_id] [-wr windows_radio]\n", args.size(), args[0].c_str());

        return 1;
    }

    //filename that contains all the images to use
    std::string filename_i_1, filename_i0, filename_i1, filename_i2;

    //Read txt file of images
    const std::string& filename_images = args[1];
    ifstream infile(filename_images);
    int num_files = 0;
    string line;
    while(getline(infile, line)){

        ++num_files;
        if (num_files == 3){
            filename_i_1  = line;
        }else{
            if (num_files == 1){
                filename_i0  = line;
            }else{
                if (num_files == 2){
                    filename_i1  = line;
                }else{
                    if (num_files == 4){
                        filename_i2  = line;
                    }
                }
            }
        }
    }
    infile.close();

    if (num_files == 3){
        fprintf(stderr, "ERROR: 3 images given as input\n");
        fprintf(stderr, "Usage: 4 images in the following order: I0, I1, I-1, I2\n");
        fprintf(stderr, "Usage: 2 images in the following order: I0, I1\n");
        return 1;
    }

    //Save other arguments
    const std::string& filename_go  = args[2];
    const std::string& filename_ba  = args[3];
    const std::string& filename_out = args[4];
    const std::string& filenme_sim  = args[5];
    const std::string& filename_occ = args[6];
    const char *filename_sal0 = nullptr;
    const char *filename_sal1 = nullptr;

    if (args.size() == 9){
        filename_sal0 = args[7].c_str();
        filename_sal1 = args[8].c_str();
    }

    //Optional arguments
    int w_radio = stoi(windows_ratio);
    int val_method = stoi(var_reg);

    // Open input images and .flo
    // pd: number of channels
    int w[8], h[8], pd[6];

    //Frame t-1 and t+2
    float *i_1 = nullptr;
    float *i2 = nullptr;
    if (num_files == 4){
        i_1 = iio_read_image_float_split(filename_i_1.c_str(), w + 6, h + 6, pd + 4);
        i2 = iio_read_image_float_split(filename_i2.c_str(), w + 7, h + 7, pd + 5);
    }else{
        i_1 = iio_read_image_float_split(filename_i0.c_str(), w + 6, h + 6, pd + 4);
        i2 = iio_read_image_float_split(filename_i1.c_str(), w + 7, h + 7, pd + 5);
    }

    //Frames t and t+1
    float *i0 = iio_read_image_float_split(filename_i1.c_str(), w + 0, h + 0, pd + 0);
    float *i1 = iio_read_image_float_split(filename_i2.c_str(), w + 1, h + 1, pd + 1);

    //Sparse Optical flow forward and backward
    float *go = iio_read_image_float_split(filename_go.c_str(), w + 2, h + 2, pd + 2);
    float *ba = iio_read_image_float_split(filename_ba.c_str(), w + 3, h + 3, pd + 3);


    //Ensure dimensions match in images
    if (num_files == 4){
        if (w[0] != w[1] || h[0] != h[1] || w[0] != w[6] || h[0] != h[6])
            return fprintf(stderr, "ERROR: input images size mismatch\n");

        if (w[0] != w[7] || h[0] != h[7] || w[1] != w[6] || h[1] != h[6])
            return fprintf(stderr, "ERROR: input images size mismatch\n");

        if (w[1] != w[7] || h[1] != h[7] || w[6] != w[7] || h[6] != h[7])
            return fprintf(stderr, "ERROR: input images size mismatch\n");

    }else{
        if (w[0] != w[1] || h[0] != h[1])
            return fprintf(stderr, "ERROR: input images size mismatch\n");
    }

    //Ensure dimensions match in flow
    if (w[2] != w[3] || h[2] != h[3] || pd[2] != 2 || pd[2] != pd[3])
        return fprintf(stderr, "ERROR: input flow field size mismatch\n");

    //Load or compute saliency
    float *sal0 = nullptr;
    float *sal1 = nullptr;
    if (args.size() == 9){
        sal0 = iio_read_image_float(filename_sal0, w + 4, h + 4);
        sal1 = iio_read_image_float(filename_sal1, w + 5, h + 5);
        fprintf(stderr, "Reading saliency values given\n");
    }else{
        fprintf(stderr, "Saliency values not given\n");
        sal0 = new float[w[0]*h[0]];
        sal1 = new float[w[0]*h[0]];
        for (int i = 0; i < w[0]*h[0]; i++){
            sal0[i] = 1.0;
            sal1[i] = 1.0;
        }
    }

    for (int i = 0; i < w[0]*h[0]; i++){
        sal0[i] = 1.0;
        sal1[i] = 1.0;
    }
    for (int i = 0; i < w[0]*h[0]; i++){
        assert(std::isfinite(sal0[i]));
        assert(std::isfinite(sal1[i]));
        if (sal0[i] < 0.0)
            fprintf(stderr, "cosa: %d\n", i);
        if (sal1[i] < 0.0)
            fprintf(stderr, "cosa: %d\n", i);
    }

    //Initialize output optical flow and energy
    float *out_flow = new float[w[0]*h[0]*2];
    float *out_occ = new float[w[0]*h[0]];
    float *ene_val = new float[w[0]*h[0]];

    for (int i = 0; i < w[0]*h[0]*2; i++){
        out_flow[i] = NAN;
    }


    // Print method used
    if (num_files == 2 && val_method == M_TVL1_OCC){
        //If only two images given for occ, something not working
        //TODO: When new methods with occlusions implemented, add here
        switch (val_method) {
        case M_TVL1_OCC:
            fprintf(stderr, "Since only two images given, method is changed to TV-l2 coupled\n");
            val_method  = M_TVL1;
            break;
        default:
            fprintf(stderr, "Method unknown\n");
            break;
        }

    }else{
        //If four images given for without occ, two not needed
        if(num_files == 4 && val_method >= 0 && val_method <= 7) {
            fprintf(stderr, "Only two of the four images given will be used, according to method selected\n");
            fprintf(stderr, "Method: ");
            switch(val_method){
            case M_NLTVL1: //NLTVL1
                fprintf(stderr, "NLTV-L1\n");
                break;
            case M_TVCSAD: //TV-CSAD
                fprintf(stderr, "TV-CSAD\n");
                break;
            case M_NLTVCSAD: //NLTV-CSAD
                fprintf(stderr, "NLTV-CSAD\n");
                break;
            case M_TVL1_W: //TV-l2 con pesos
                fprintf(stderr, "TV-l2 coupled Weights\n");
                break;
            case M_NLTVCSAD_W: //NLTV-CSAD con pesos
                fprintf(stderr, "NLTV-CSAD Weights\n");
                break;
            case M_NLTVL1_W: //NLTV-L1 con pesos
                fprintf(stderr," NLTV-L1 Weights\n");
                break;
            case M_TVCSAD_W: //TV-CSAD con pesos
                fprintf(stderr, "TV-CSAD Weights\n");
                break;
            default: //TV-l2 coupled
                fprintf(stderr, "TV-l2 coupled\n");
            }
        }else{
            fprintf(stderr, "Method: ");
            switch (val_method) {
            case M_TVL1_OCC: //TV-l2 with occlusion
                fprintf(stderr, "TV-l2 occlusions\n");
                break;
            default:
                break;
            }
        }
    }
    //Initialize parameters
    int step_alg = LOCAL_STEP;
    Parameters params =  init_params(file_params, step_alg);
    params.w = w[0];
    params.h = h[0];
    params.w_radio = w_radio;
    params.val_method = val_method;
    cerr << params;

    //Match growing algorithm
    match_growing_variational(go, ba, i0, i1, i_1, i2, sal0, sal1, params, ene_val, out_flow, out_occ);


    // Save results
    iio_save_image_float_split(filename_out.c_str(), out_flow, w[0], h[0], 2);
    iio_save_image_float(filenme_sim.c_str(), ene_val, w[0], h[0]);
    iio_save_image_float(filename_occ.c_str(), out_occ, w[0], h[0]);

    // cleanup and exit
    free(i_1);
    free(i0);
    free(i1);
    free(i2);

    free(go);
    free(ba);

    if (args.size() ==  8){ //c == 8
        free(sal0);
        free(sal1);
    }else{
        delete [] sal0;
        delete [] sal1;
    }

    delete [] out_flow;
    delete [] ene_val;
    delete [] out_occ;

    today = system_clock::now();

    tt = system_clock::to_time_t ( today );
    std::cerr << "today is: " << ctime(&tt);
    return 0;
}
#endif
