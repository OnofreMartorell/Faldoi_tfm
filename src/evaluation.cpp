// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license athis program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017, Onofre Martorell <onofremartorelln@gmail.com>
// All rights reserved.


#include <cmath>
#include <vector>
#include <cstring>
#include <iostream>
#include <string>

extern "C" {
#include "iio.h"
}

void  of_errors(float *flow1, float *flow2, int w, int h){

    //[ni, nj, nChannels] = size(flow1);
    auto image_ae = new float[w*h];
    auto image_ee = new float [w*h];
    float ae = 0;
    float ee = 0;
    int nPix = 0;

    for (int j = 0; j < h; j++){
        for (int i = 0; i < w; i++){

            //Exactly as in C code
            int p = i*h + j;
            float u_c = flow1[p];//(i,j,1);
            float v_c = flow1[p + w*h];//(i,j,2);

            float u_e = flow2[p];//(i,j,1);
            float v_e = flow2[p + w*h];//(i,j,2);

            if ((std::abs(u_c) < 10000) && (std::abs(v_c) < 10000)){
                nPix++;
                float esc_prod = u_c*u_e + v_c*v_e + 1.0;

                float term1 = (u_c * u_c + v_c * v_c + 1.0);
                float term2 = (u_e * u_e + v_e * v_e + 1.0);
                float term3 = std::sqrt(term1 * term2);

                float tmp = esc_prod / term3;


                tmp = std::fmin(tmp, 1.0);
                tmp = std::fmax(tmp, -1.0);

                ae = ae + acos(tmp);
                image_ae[p] = acos(tmp);
                ee = ee + std::sqrt((u_c - u_e)*(u_c - u_e) + (v_c - v_e)*(v_c - v_e));
                image_ee[p] = std::sqrt((u_c - u_e)*(u_c - u_e) + (v_c - v_e)*(v_c - v_e));
            }
        }
    }
    ae = ae/nPix;
    std::cout << "Angular error: "<< ae << "\n";
    ee = ee/nPix;
    std::cout << "Endpoint error: "<< ee << "\n";
}


int main(int argc, char* argv[]){

    std::vector<std::string> args(argv, argv + argc);
    if (args.size() != 3) {
        fprintf(stderr, "usage %lu :\n\t%s out.flo gt.flo", args.size(), args[0].c_str());

        return 1;
    }


    //Save other arguments
    const std::string& filename_of  = args[1];
    const std::string& filename_gt  = args[2];




    // Open input images and .flo
    // pd: number of channels
    int w[2], h[2], pd[2];



    //Sparse Optical flow forward and backward
    float *optical_flow = iio_read_image_float_split(filename_of.c_str(), w, h, pd);
    float *ground_truth = iio_read_image_float_split(filename_gt.c_str(), w + 1, h + 1, pd + 1);
    of_errors(optical_flow, ground_truth, w[0], h[0]);
    return 0;

}
