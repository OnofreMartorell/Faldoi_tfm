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
#include <fstream>
extern "C" {
#include "iio.h"
}


using namespace std;

void  of_error(
        const std::string& filename_of,
        const std::string& filename_gt,
        float *results
        ){


    // Open input images and .flo
    // pd: number of channels
    int w[2], h[2], pd[2];



    //Sparse Optical flow forward and backward
    float *flow1 = iio_read_image_float_split(filename_of.c_str(), w, h, pd);
    float *flow2 = iio_read_image_float_split(filename_gt.c_str(), w + 1, h + 1, pd + 1);



    //[ni, nj, nChannels] = size(flow1);
    auto image_ae = new float[w[0]*h[0]];
    auto image_ee = new float[w[0]*h[0]];
    float ae = 0;
    float ee = 0;
    int nPix = 0;

    for (int j = 0; j < h[0]; j++){
        for (int i = 0; i < w[0]; i++){

            //Exactly as in C code
            int p = i*h[0] + j;
            float u_c = flow1[p];//(i,j,1);
            float v_c = flow1[p + w[0]*h[0]];//(i,j,2);

            float u_e = flow2[p];//(i,j,1);
            float v_e = flow2[p + w[0]*h[0]];//(i,j,2);

            if ((abs(u_c) < 10000) && (abs(v_c) < 10000)){
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
    cout << "num pixels: " << nPix << "\n";
    ae = ae/nPix;
    std::cout << "Angular error: "<< ae << "\n";
    ee = ee/nPix;
    std::cout << "Endpoint error: "<< ee << "\n";
    if (nPix == 0){
        results[0] = 0;
        results[1] = 0;
    }else{
    results[0] = ae;
    results[1] = ee;
    }
}


int main(int argc, char* argv[]){

    std::vector<std::string> args(argv, argv + argc);
    if (args.size() != 3) {
        fprintf(stderr, "usage %lu :\n\t%s out.txt gt.txt", args.size(), args[0].c_str());

        return 1;
    }


    //Save other arguments
    //    const std::string& filename_of  = args[1];
    //    const std::string&  filename_gt = args[2];

    const string& filename_images_of = args[1];
    ifstream infile_of(filename_images_of);
    string line_of;

    const string& filename_images_gt = args[2];
    ifstream infile_gt(filename_images_gt);
    string line_gt;

    float global_ae = 0.0;
    float global_ee = 0.0;
    float num_files = 0.0;

    while(getline(infile_of, line_of)){

        num_files = num_files + 1.0;
        getline(infile_gt, line_gt);
        float *results = new float[2];
        const string& filename_of = line_of;
        const string& filename_gt = line_gt;

        of_error(filename_of, filename_gt, results);
        global_ae += results[0];
        global_ee += results[1];

    }

    infile_of.close();
    infile_gt.close();
    global_ae = global_ae/num_files;
    std::cout << "Global angular error: "<< global_ae << "\n";
    global_ee = global_ee/num_files;
    std::cout << "Global endpoint error: "<< global_ee << "\n";


    //    float *results = new float[2];

    //    of_error(filename_of, filename_gt, results);
    return 0;

}
