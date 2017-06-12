// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2014, Roberto P.Palomares <r.perezpalomares@gmail.com>
// All rights reserved.
#ifndef ENERGY_MODEL
#define ENERGY_MODEL

#include <cmath>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <algorithm>


extern "C" {
#include "bicubic_interpolation.h"
//#include "mask.h"
}

#include "utils.h"
#include "energy_structures.h"
#include "aux_energy_model.h"
#include "energy_model.h"

//Models
#include "tvl2_model.h"
#include "nltv_model.h"
#include "tvcsad_model.h"
#include "nltvcsad_model.h"

//Models with weights
#include "tvl2w_model.h"
#include "nltvw_model.h"
#include "nltvcsadw_model.h"
#include "tvcsadw_model.h"

//Models with occlusions
#include "tvl2_model_occ.h"

//IMAGE_PARAMETERS
#define PRESMOOTHING_SIGMA  0.90

//OPTICAL FLOW PARAMETERS
#define PAR_DEFAULT_LAMBDA  40//40
#define PAR_DEFAULT_THETA   0.3
#define PAR_DEFAULT_TAU     0.125 //0.25
#define PAR_DEFAULT_NWARPS  1  //5
#define PAR_DEFAULT_TOL_D   0.01
#define PAR_DEFAULT_VERBOSE 0  //0
#define PAR_DEFAULT_GAMMA 0.05  //0


void rgb_to_gray(const float *in, const int w, const int h, float *out){
    const int size = w*h;

    for (int i = 0; i < size; i++){

        out[i] = .299*in[i] + .587*in[size + i] + .114*in[2*size + i];

    }

}

float pow2( const float f ) {return f*f;}


//Asume que las imagenes no estan normalizadas
void rgb_to_lab(const float *in, int size, float *out){

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
        float correct_lab = exp(-color_attenuation*pow2(pow2(L/100) - 0.6));
        out[i] = L;
        out[i + size] = A*correct_lab;
        out[i + 2*size] = B*correct_lab;
    }
}



////AUXILIAR FUNCTIONS (NEUMANN BOUNDARY CONDITIONS)///////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AUXILIAR FUNCTIONS CSAD///////////////////

///////////////////////////////////////////////////////////

void initialize_auxiliar_stuff(
        SpecificOFStuff& ofStuff,
        OpticalFlowData& ofCore){
    switch(ofCore.method)
    {
    case M_NLTVL1: //NLTV-L1
        intialize_stuff_nltvl1(&ofStuff, &ofCore);
        break;
    case M_TVCSAD: //TV-CSAD
        intialize_stuff_tvcsad(&ofStuff, &ofCore);
        break;
    case M_NLTVCSAD: //NLTV-CSAD
        intialize_stuff_nltvcsad(&ofStuff, &ofCore);
        break;
    case M_TVL1_W: //TV-l2 coupled con pesos
        intialize_stuff_tvl2coupled_w(&ofStuff, &ofCore);
        break;
    case M_NLTVCSAD_W: //NLTV-CSAD con pesos
        intialize_stuff_nltvcsad_w(&ofStuff, &ofCore);
        break;
    case M_NLTVL1_W: //NLTV-L1 con pessos
        intialize_stuff_nltvl1_w(&ofStuff, &ofCore);
        break;
    case M_TVCSAD_W: //TV-CSAD con pesos
        intialize_stuff_tvcsad_w(&ofStuff, &ofCore);
        break;
    case M_TVL1_OCC: //TV-l2 with occlusion
        intialize_stuff_tvl2coupled_occ(ofStuff, ofCore);
        break;
    default: //TV-l2 coupled
        intialize_stuff_tvl2coupled(&ofStuff, &ofCore);
    }

}

void free_auxiliar_stuff(SpecificOFStuff *ofStuff, OpticalFlowData *ofCore)
{
    const int method = ofCore->method;

    switch(method)
    {
    case M_NLTVL1: //NLTVL1
        free_stuff_nltvl1(ofStuff);
        break;
    case M_TVCSAD: //TV-CSAD
        free_stuff_tvcsad(ofStuff);
        break;
    case M_NLTVCSAD: //NLTV-CSAD
        free_stuff_nltvcsad(ofStuff);
        break;
    case M_TVL1_W: //TV-l2 con pesos
        free_stuff_tvl2coupled_w(ofStuff);
        break;
    case M_NLTVCSAD_W: //NLTV-CSAD con pesos
        free_stuff_nltvcsad_w(ofStuff);
        break;
    case M_NLTVL1_W: //NLTV-L1 con pesos
        free_stuff_nltvl1_w(ofStuff);
        break;
    case M_TVCSAD_W: //TV-CSAD con pesos
        free_stuff_tvcsad_w(ofStuff);
        break;
    case M_TVL1_OCC: //TV-l2 with occlusion
        free_stuff_tvl2coupled_occ(ofStuff);
        break;
    default: //TV-l2 coupled
        free_stuff_tvl2coupled(ofStuff);
    }
}


void prepare_stuff(
        SpecificOFStuff *ofStuff1,
        OpticalFlowData *ofCore1,
        SpecificOFStuff *ofStuff2,
        OpticalFlowData *ofCore2,
        const float *i0,
        const float *i1,
        const float *i_1,
        const float *i2,
        const int pd,
        float **out_i0,
        float **out_i1,
        float **out_i_1,
        float **out_i2
        ) {

    const int method = ofCore1->method;
    const int w = ofCore1->w;
    const int h = ofCore1->h;

    switch(method){
    case M_NLTVL1: //NLTV-L1
    {
        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];

        float *alb = new float[w*h*pd];
        float *blb = new float[w*h*pd];

        const int n_d = NL_DUAL_VAR;
        const int radius = NL_BETA;

        rgb_to_lab(i0, w*h, alb);
        rgb_to_lab(i1, w*h, blb);
        // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
        nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                ofStuff1->nltvl1.p, ofStuff1->nltvl1.q);
        nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                ofStuff2->nltvl1.p, ofStuff2->nltvl1.q);
        if (pd != 1){

            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i1, w, h, b_tmp);

        }else{
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }
        // normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->nltvl1.I1x, ofStuff1->nltvl1.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->nltvl1.I1x, ofStuff2->nltvl1.I1y,
                          ofCore2->w, ofCore2->h);

        *out_i0 = a_tmp;
        *out_i1 = b_tmp;

        delete [] alb;
        delete [] blb;
    }
        break;
    case M_TVCSAD: //TVCSAD
    {
        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];
        const int rdt = DT_R;
        const int ndt = DT_NEI;
        std::printf("1 - Inicializado CSAD\n");
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->tvcsad.pnei);
        std::printf("2 - Inicializado CSAD\n");
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->tvcsad.pnei);
        if (pd != 1)
        {
            // std::printf("Numero canales:%d\n",pd);
            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i1, w, h, b_tmp);
        }
        else
        {
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }
        // normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->tvcsad.I1x, ofStuff1->tvcsad.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->tvcsad.I1x, ofStuff2->tvcsad.I1y,
                          ofCore2->w, ofCore2->h);
        *out_i0 = a_tmp;
        *out_i1 = b_tmp;
        std::printf("Salimos de CSAD\n");
    }
        break;
    case M_NLTVCSAD: //NLTV-CSAD
    {
        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];

        float *alb = new float[w*h*pd];
        float *blb = new float[w*h*pd];

        int n_d = NL_DUAL_VAR;
        int radius = NL_BETA;
        int rdt = DT_R;
        int ndt = DT_NEI;
        std::printf("Preparando CSAD\n");
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->nltvcsad.pnei);
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->nltvcsad.pnei);

        rgb_to_lab(i0, w*h, alb);
        rgb_to_lab(i1, w*h, blb);
        // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
        nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                ofStuff1->nltvcsad.p, ofStuff1->nltvcsad.q);
        nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                ofStuff2->nltvcsad.p, ofStuff2->nltvcsad.q);
        std::printf("Preparado NLTV\n");
        if (pd!=1)
        {
            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i0, w, h, b_tmp);
        }
        else
        {
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }
        // normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->nltvcsad.I1x, ofStuff1->nltvcsad.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->nltvcsad.I1x, ofStuff2->nltvcsad.I1y,
                          ofCore2->w, ofCore2->h);

        *out_i0 = a_tmp;
        *out_i1 = b_tmp;

        delete [] alb;
        delete [] blb;

        std::printf("Salimos NLTV_CSAD\n");

    }
        break;
    case M_TVL1_W://TV-l2 coupled con pesos
    {

        float *weight1 = ofStuff1->tvl2w.weight;
        float *weight2 = ofStuff2->tvl2w.weight;
        gaussian1Dweight(weight1, ofCore1->wr);
        gaussian1Dweight(weight2, ofCore2->wr);

        std::printf("Pesos\n");
        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];
        if (pd!=1)
        {
            // std::printf("Numero canales:%d\n",pd);
            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i1, w, h, b_tmp);
        }
        else
        {
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }
        // normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->tvl2w.I1x, ofStuff1->tvl2w.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->tvl2w.I1x, ofStuff2->tvl2w.I1y,
                          ofCore2->w, ofCore2->h);
        *out_i0 = a_tmp;
        *out_i1 = b_tmp;

    }
        break;
    case M_NLTVCSAD_W: //NLTV-CSAD con pesos
    {

        float *weight1 = ofStuff1->nltvcsadw.weight;
        float *weight2 = ofStuff2->nltvcsadw.weight;
        gaussian1Dweight(weight1, ofCore1->wr);
        gaussian1Dweight(weight2, ofCore2->wr);

        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];

        float *alb = new float[w*h*pd];
        float *blb = new float[w*h*pd];

        int n_d = NL_DUAL_VAR;
        int radius = NL_BETA;
        int rdt = DT_R;
        int ndt = DT_NEI;
        std::printf("Preparado CSAD\n");
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->nltvcsadw.pnei);
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->nltvcsadw.pnei);

        rgb_to_lab(i0, w*h, alb);
        rgb_to_lab(i1, w*h, blb);
        // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
        nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                ofStuff1->nltvcsadw.p, ofStuff1->nltvcsadw.q);
        nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                ofStuff2->nltvcsadw.p, ofStuff2->nltvcsadw.q);
        std::printf("Preparado NLTV\n");
        if (pd!=1)
        {
            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i1, w, h, b_tmp);
        }
        else
        {
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }
        // normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->nltvcsadw.I1x, ofStuff1->nltvcsadw.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->nltvcsadw.I1x, ofStuff2->nltvcsadw.I1y,
                          ofCore2->w, ofCore2->h);

        *out_i0 = a_tmp;
        *out_i1 = b_tmp;

        delete [] alb;
        delete [] blb;

        std::printf("Salimos NLTV_CSAD con pesos\n");

    }
        break;
    case M_NLTVL1_W: //NLTV-L1 con pesos
    {
        float *weight1 = ofStuff1->nltvl1w.weight;
        float *weight2 = ofStuff2->nltvl1w.weight;
        gaussian1Dweight(weight1, ofCore1->wr);
        gaussian1Dweight(weight2, ofCore2->wr);
        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];

        float *alb = new float[w*h*pd];
        float *blb = new float[w*h*pd];

        int n_d = NL_DUAL_VAR;
        int radius = NL_BETA;

        rgb_to_lab(i0, w*h, alb);
        rgb_to_lab(i1, w*h, blb);
        // std::printf("W:%d x H:%d\n Neir:%d, radius:%d\n",w,h,n_d,radius);
        nltv_ini_dual_variables(alb, pd, w, h, n_d, radius,
                                ofStuff1->nltvl1w.p, ofStuff1->nltvl1w.q);
        nltv_ini_dual_variables(blb, pd, w, h, n_d, radius,
                                ofStuff2->nltvl1w.p, ofStuff2->nltvl1w.q);
        if (pd!=1)
        {
            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i1, w, h, b_tmp);
        }
        else
        {
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }
        // normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->nltvl1w.I1x, ofStuff1->nltvl1w.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->nltvl1w.I1x, ofStuff2->nltvl1w.I1y,
                          ofCore2->w, ofCore2->h);

        *out_i0 = a_tmp;
        *out_i1 = b_tmp;

        delete [] alb;
        delete [] blb;
    }
        break;
    case M_TVCSAD_W: //TVCSAD con pesos
    {
        float *weight1 = ofStuff1->tvcsadw.weight;
        float *weight2 = ofStuff2->tvcsadw.weight;
        gaussian1Dweight(weight1, ofCore1->wr);
        gaussian1Dweight(weight2, ofCore2->wr);
        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];
        int rdt = DT_R;
        int ndt = DT_NEI;
        std::printf("1 - Inicializado CSAD\n");
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff1->tvcsadw.pnei);
        std::printf("2 - Inicializado CSAD\n");
        csad_ini_pos_nei(w, h, ndt, rdt, ofStuff2->tvcsadw.pnei);
        if (pd != 1) {
            // std::printf("Numero canales:%d\n",pd);
            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i1, w, h, b_tmp);
        }
        else
        {
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }
        // normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->tvcsadw.I1x, ofStuff1->tvcsadw.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->tvcsadw.I1x, ofStuff2->tvcsadw.I1y,
                          ofCore2->w, ofCore2->h);
        *out_i0 = a_tmp;
        *out_i1 = b_tmp;
        std::printf("Salimos de CSAD\n");
    }
        break;
    case M_TVL1_OCC:
    {
        float *i1_tmp = new float[w*h];
        float *i2_tmp = new float[w*h];
        float *i0_tmp = new float[w*h];
        float *i_1_tmp = new float[w*h];

        //Check if image is gray, otherwise transform
        if (pd != 1){
            rgb_to_gray(i0, w, h, i0_tmp);
            rgb_to_gray(i1, w, h, i1_tmp);
            rgb_to_gray(i_1, w, h, i_1_tmp);
            rgb_to_gray(i2, w, h, i2_tmp);

        }else{
            memcpy(i0_tmp, i0, w*h*sizeof(float));
            memcpy(i1_tmp, i1, w*h*sizeof(float));
            memcpy(i_1_tmp, i_1, w*h*sizeof(float));
            memcpy(i2_tmp, i2, w*h*sizeof(float));
        }

        // Normalize the images between 0 and 255
        image_normalization_4(i0_tmp, i1_tmp, i_1_tmp, i2_tmp, i0_tmp, i1_tmp, i_1_tmp, i2_tmp, w*h);

        //Make a little smoth of images
        gaussian(i1_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(i2_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(i0_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(i_1_tmp, w, h, PRESMOOTHING_SIGMA);

        //TODO: change the computation of derivatives
        centered_gradient(i1_tmp, ofStuff1->tvl2_occ.I1x, ofStuff1->tvl2_occ.I1y, w, h);
        centered_gradient(i_1_tmp, ofStuff1->tvl2_occ.I_1x, ofStuff1->tvl2_occ.I_1y, w, h);

        centered_gradient(i0_tmp, ofStuff2->tvl2_occ.I1x, ofStuff2->tvl2_occ.I1y, w, h);
        centered_gradient(i2_tmp, ofStuff2->tvl2_occ.I_1x, ofStuff2->tvl2_occ.I_1y, w, h);

        // Initialize g (weight)
        // The derivatives are taken from previous computation. Forward: I0, backward: I1
        init_weight(ofStuff1->tvl2_occ.g, ofStuff2->tvl2_occ.I1x, ofStuff2->tvl2_occ.I1y, w*h);
        init_weight(ofStuff2->tvl2_occ.g, ofStuff1->tvl2_occ.I1x, ofStuff1->tvl2_occ.I1y, w*h);

        //The following variables contain a gray and smooth version
        //of the corresponding image
        *out_i1 = i1_tmp;
        *out_i2 = i2_tmp;
        *out_i0 = i0_tmp;
        *out_i_1 = i_1_tmp;
    }
        break;
    default: //TV-l2 coupled
    {
        float *a_tmp = new float[w*h];
        float *b_tmp = new float[w*h];
        //Check if image is gray, otherwise transform
        if (pd != 1){
            // std::printf("Numero canales:%d\n",pd);
            rgb_to_gray(i0, w, h, a_tmp);
            rgb_to_gray(i1, w, h, b_tmp);
        }else{
            memcpy(a_tmp, i0, w*h*sizeof(float));
            memcpy(b_tmp, i1, w*h*sizeof(float));
        }

        // Normalize the images between 0 and 255
        image_normalization(a_tmp, b_tmp, a_tmp, b_tmp, w*h);
        gaussian(a_tmp, w, h, PRESMOOTHING_SIGMA);
        gaussian(b_tmp, w, h, PRESMOOTHING_SIGMA);
        centered_gradient(b_tmp, ofStuff1->tvl2.I1x, ofStuff1->tvl2.I1y,
                          ofCore1->w, ofCore1->h);
        centered_gradient(a_tmp, ofStuff2->tvl2.I1x, ofStuff2->tvl2.I1y,
                          ofCore2->w, ofCore2->h);
        *out_i0 = a_tmp;
        *out_i1 = b_tmp;

    }
    }
}


void of_estimation(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore,
        float *ener_N,
        const float *i0,  //first frame
        const float *i1,  //second frame
        const float *i_1,
        const PatchIndexes index // end row
        ) {

    //TODO: Each method should have its own set of parameters
    float lambda  = PAR_DEFAULT_LAMBDA;
    float theta   = PAR_DEFAULT_THETA;
    float tau     = PAR_DEFAULT_TAU;
    float tol_OF  = PAR_DEFAULT_TOL_D;
    int    verbose = PAR_DEFAULT_VERBOSE;
    int    warps   = PAR_DEFAULT_NWARPS;

    switch(ofCore->method) {
    case M_NLTVL1: //NLTV-L1
    {
        lambda = 2;
        theta = 0.3;
        tau = 0.1;
        //estimate_tvl1
        guided_nltvl1(i0, i1, ofCore, &(ofStuff->nltvl1), ener_N, index,
                      lambda, theta, tau, tol_OF, warps, verbose);

    }
        break;
    case M_TVCSAD: //TVCSAD
    {
        lambda = 0.85; //80/(49-2)
        theta = 0.3;
        tau = 0.1;
        //estimate_tvl1
        guided_tvcsad(i0, i1, ofCore, &(ofStuff->tvcsad), ener_N, index.ii, index.ij, index.ei, index.ej,
                      lambda, theta, tau, tol_OF, warps, verbose);
    }
        break;
    case M_NLTVCSAD: //NLTV-CSAD
    {
        lambda = 0.85; //80/(49-2)
        theta = 0.3;
        tau = 0.1;
        //estimate_tvl1
        guided_nltvcsad(i0, i1, ofCore, &(ofStuff->nltvcsad), ener_N, index.ii, index.ij, index.ei, index.ej,
                        lambda, theta, tau, tol_OF, warps, verbose);
    }
        break;
    case M_TVL1_W: //TV-l2 coupled con pesos
    {
        const float central = ofStuff->tvl2w.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        guided_tvl2coupled_w(i0, i1, ofCore, &(ofStuff->tvl2w), ener_N, index.ii, index.ij, index.ei, index.ej,
                             lambda, theta, tau, tol_OF, warps, verbose);
    }
        break;
    case M_NLTVCSAD_W: //NLTV-CSAD con pesos
    {
        lambda = 0.85; //80/(49-2)
        const float central = ofStuff->nltvcsadw.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        theta = 0.3;
        tau = 0.1;
        //estimate_tvl1
        guided_nltvcsad_w(i0, i1, ofCore, &(ofStuff->nltvcsadw), ener_N, index.ii, index.ij, index.ei, index.ej,
                          lambda, theta, tau, tol_OF, warps, verbose);
    }
        break;
    case M_NLTVL1_W: //NLTV-L1 cons pesos
    {
        lambda = 2;
        lambda = 0.85; //80/(49-2)
        const float central = ofStuff->nltvl1w.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        theta = 0.3;
        tau = 0.1;
        //estimate_tvl1
        guided_nltvl1_w(i0, i1, ofCore, &(ofStuff->nltvl1w), ener_N, index.ii, index.ij, index.ei, index.ej,
                        lambda, theta, tau, tol_OF, warps, verbose);

    }
        break;
    case M_TVCSAD_W: //TVCSAD con pesos
    {
        lambda = 0.85; //80/(49-2)
        const float central = ofStuff->tvcsadw.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        theta = 0.3;
        tau = 0.1;
        //estimate_tvl1
        guided_tvcsad_w(i0, i1, ofCore, &(ofStuff->tvcsadw), ener_N, index.ii, index.ij, index.ei, index.ej,
                        lambda, theta, tau, tol_OF, warps, verbose);
    }
        break;
    case M_TVL1_OCC:
    {
        lambda = 0.25;
        theta = 0.3;
        const float beta = 1;
        const float alpha = 0.01;
        const float tau_u = 0.125;
        const float tau_eta = 0.125;
        const float tau_chi = 0.125;
        //estimate_tvl2 with occlusions
        guided_tvl2coupled_occ(i0, i1, i_1, ofCore, &(ofStuff->tvl2_occ), ener_N, index,
                           lambda, theta, tau_u, tau_eta, tau_chi, beta, alpha, tol_OF, warps, verbose);
    }
        break;
    default: //TV-l2 coupled
        //estimate_tvl2
        guided_tvl2coupled(i0, i1, ofCore, &(ofStuff->tvl2), ener_N, index.ii, index.ij, index.ei, index.ej,
                           lambda, theta, tau, tol_OF, warps, verbose);
    }
}

///////////////////////////////////////////////////////////////////////



void eval_functional(SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore,
        float *ener_N,
        const float *i0,  //first frame
        const float *i1,  //second frame
        const float *i_1,
        const PatchIndexes index) {
    //TODO: Each method should have its own set of parameters

    float lambda  = PAR_DEFAULT_LAMBDA;
    float theta   = PAR_DEFAULT_THETA;

    switch(ofCore->method){
    case M_NLTVL1: //NLTV-L1
    {
        lambda  = 2.0;
        eval_nltvl1(i0, i1, ofCore, &(ofStuff->nltvl1), ener_N,
                    index, lambda, theta);
    }
        break;
    case M_TVCSAD: //TV-CSAD
    {
        lambda  = 0.85;
        eval_tvcsad(i0, i1, ofCore, &(ofStuff->tvcsad), ener_N,
                    index.ii, index.ij, index.ei, index.ej, lambda, theta);
    }
        break;
    case M_NLTVCSAD: //NLTV-CSAD
    {
        lambda  = 0.85;
        eval_nltvcsad(i0, i1, ofCore, &(ofStuff->nltvcsad), ener_N,
                      index.ii, index.ij, index.ei, index.ej, lambda, theta);
    }
        break;
    case M_TVL1_W: //TV-l2 con pesos
    {
        const float central = ofStuff->tvl2w.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        eval_tvl2coupled_w(i0, i1, ofCore, &(ofStuff->tvl2w), ener_N,
                           index.ii, index.ij, index.ei, index.ej, lambda, theta);
    }
        break;
    case M_NLTVCSAD_W: //NLTV-CSAD con pesos
    {
        lambda  = 0.85;
        const float central = ofStuff->nltvcsadw.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        eval_nltvcsad_w(i0, i1, ofCore, &(ofStuff->nltvcsadw), ener_N,
                        index.ii, index.ij, index.ei, index.ej, lambda, theta);
    }
        break;
    case M_NLTVL1_W: //NLTV-L1 con pesos
    {
        lambda  = 2.0;
        const float central = ofStuff->nltvl1w.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        eval_nltvl1_w(i0, i1, ofCore, &(ofStuff->nltvl1w), ener_N,
                      index.ii, index.ij, index.ei, index.ej, lambda, theta);
    }
        break;
    case M_TVCSAD_W: //TV-CSAD
    {
        lambda  = 0.85;
        const float central = ofStuff->tvcsadw.weight[ofCore->wr + 1];
        lambda  = lambda /(central*central);
        eval_tvcsad_w(i0, i1, ofCore, &(ofStuff->tvcsadw), ener_N,
                      index.ii, index.ij, index.ei, index.ej, lambda, theta);
    }
        break;
    case M_TVL1_OCC:
    {
        lambda = 0.25;
        theta = 0.3;
        const float beta = 1;
        const float alpha = 0.01;
        eval_tvl2coupled_occ(i0, i1, i_1, ofCore,
                &(ofStuff->tvl2_occ),
                ener_N,
                index,
                lambda,  // weight of the data term
                theta,
                alpha,
                beta);
    }
        break;
    default: //TV-l2 coupled
    {
        eval_tvl2coupled(i0, i1, ofCore, &(ofStuff->tvl2), ener_N,
                         index.ii, index.ij, index.ei, index.ej, lambda, theta);
    }
    }
}


/////////////////////////////////////////////////////////////////////////
///////////////////////////OCCLUSION ESTIMATION/////////////////////////
///////////////////////////////////////////////////////////////////////



void eval_functional_occ(
        SpecificOFStuff *ofStuff,
        OpticalFlowData *ofCore,
        float *ener_N,
        float *a,  //first frame
        float *b,  //second frame
        const int ii, // initial column
        const int ij, // initial row
        const int ei, // end column
        const int ej // end row,
        ) {
    //TODO: Each method should have its own set of parameters

    float lambda  = PAR_DEFAULT_LAMBDA;
    float theta   = PAR_DEFAULT_THETA;

    switch(ofCore->method){

    case M_TVL1_OCC:
    {
        eval_tvl2coupled(a, b, ofCore, &(ofStuff->tvl2), ener_N,
                         ii, ij, ei, ej, lambda, theta);
    }
        break;
    default:
        break;
    }
}

#endif //ENERGY_MODEL
