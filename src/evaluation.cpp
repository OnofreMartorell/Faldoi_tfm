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


#include "energy_structures.h"
#include "utils_preprocess.h"
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



    float *image_ee = new float[w[0]*h[0]];
    float *image_ae = new float[w[0]*h[0]];
    float ae = 0;
    float ee = 0;
    int nPix = 0;

    for (int j = 0; j < h[0]; j++){
        for (int i = 0; i < w[0]; i++){


            int p = i*h[0] + j;
            float u_c = flow1[p];
            float v_c = flow1[p + w[0]*h[0]];

            float u_e = flow2[p];
            float v_e = flow2[p + w[0]*h[0]];

            if ((abs(u_c) < 10000) && (abs(v_c) < 10000) && (abs(u_e) < 10000) && (abs(v_e) < 10000)){

                nPix++;
                float esc_prod = u_c*u_e + v_c*v_e + 1.0;

                float term1 = (u_c * u_c + v_c * v_c + 1.0);
                float term2 = (u_e * u_e + v_e * v_e + 1.0);
                float term3 = std::sqrt(term1 * term2);

                float tmp = esc_prod / term3;


                tmp = std::fmin(tmp, 1.0);
                tmp = std::fmax(tmp, -1.0);

                float ae_single = acos(tmp);
                ae = ae + ae_single;
                image_ae[p] = ae_single;

                float ee_single = std::sqrt((u_c - u_e)*(u_c - u_e) + (v_c - v_e)*(v_c - v_e));




                ee = ee + ee_single;
                image_ee[p] = ee_single;
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
    delete [] image_ae;
    delete [] image_ee;
    delete [] flow1;
    delete [] flow2;
}

void  of_error(
        const std::string& filename_of,
        const std::string& filename_gt,
        const std::string& filename_mask,
        float *results
        ){


    // Open input .flo
    int w[3], h[3], pd[3];


    //Sparse Optical flow forward and backward
    float *flow1 = iio_read_image_float_split(filename_of.c_str(), w, h, pd);
    float *flow2 = iio_read_image_float_split(filename_gt.c_str(), w + 1, h + 1, pd + 1);
    float *mask_occ = iio_read_image_float_split(filename_mask.c_str(), w + 2, h + 2, pd + 2);

    auto image_ae = new float[w[0]*h[0]];
    auto image_ee = new float[w[0]*h[0]];
    float ae = 0;
    float ee = 0;
    int nPix = 0;

    for (int j = 0; j < h[0]; j++){
        for (int i = 0; i < w[0]; i++){


            int p = i*h[0] + j;
            float u_c = flow1[p];
            float v_c = flow1[p + w[0]*h[0]];

            float u_e = flow2[p];
            float v_e = flow2[p + w[0]*h[0]];

            if ((abs(u_c) < 10000) && (abs(v_c) < 10000) && (mask_occ[p] == 0)){

                nPix++;
                float esc_prod = u_c*u_e + v_c*v_e + 1.0;

                float term1 = (u_c * u_c + v_c * v_c + 1.0);
                float term2 = (u_e * u_e + v_e * v_e + 1.0);
                float term3 = std::sqrt(term1 * term2);

                float tmp = esc_prod / term3;


                tmp = std::fmin(tmp, 1.0);
                tmp = std::fmax(tmp, -1.0);

                float ae_single = acos(tmp);
                ae = ae + ae_single;
                image_ae[p] = ae_single;

                float ee_single = std::sqrt((u_c - u_e)*(u_c - u_e) + (v_c - v_e)*(v_c - v_e));
                ee = ee + ee_single;
                image_ee[p] = ee_single;
            }
        }
    }
    cout << "num pixels: " << nPix << "\n";
    ae = ae/nPix;
    std::cout << "Angular error matched: "<< ae << "\n";
    ee = ee/nPix;
    std::cout << "Endpoint error matched: "<< ee << "\n";
    if (nPix == 0){
        results[0] = 0;
        results[1] = 0;
    }else{
        results[0] = ae;
        results[1] = ee;
    }
    delete [] image_ae;
    delete [] image_ee;
    delete [] flow1;
    delete [] flow2;
}

//compute s0-10, s10-40 and s40+
void  of_error(
        const std::string& filename_of,
        const std::string& filename_gt,
        float *results,
        float inf_limit,
        float sup_limit,
        const std::string& method_output
        ){


    // Open input images and .flo
    // pd: number of channels
    int w[2], h[2], pd[2];


    //Sparse Optical flow forward and backward
    float *flow1 = iio_read_image_float_split(filename_of.c_str(), w, h, pd);
    float *flow2 = iio_read_image_float_split(filename_gt.c_str(), w + 1, h + 1, pd + 1);


    float *image_ee = new float[w[0]*h[0]];
    float ee = 0;
    int nPix = 0;

    for (int j = 0; j < h[0]; j++){
        for (int i = 0; i < w[0]; i++){


            int p = i*h[0] + j;
            float u_c = flow1[p];
            float v_c = flow1[p + w[0]*h[0]];

            float u_e = flow2[p];
            float v_e = flow2[p + w[0]*h[0]];

            if ((abs(u_c) < 10000) && (abs(v_c) < 10000)){
                float norm_gt = sqrt((u_e)*(u_e) + (v_e)*(v_e));


                if (norm_gt >= inf_limit && norm_gt <= sup_limit){
                    nPix++;

                    float ee_single = sqrt((u_c - u_e)*(u_c - u_e) + (v_c - v_e)*(v_c - v_e));
                    ee = ee + ee_single;
                    image_ee[p] = ee_single;
                }
            }
        }
    }
    cout << "num pixels: " << nPix << "\n";
    ee = ee/nPix;
    std::cout << method_output << ": "<< ee << "\n";
    if (nPix == 0){
        results[0] = 0;
    }else{
        results[0] = ee;
    }
    delete [] image_ee;
    delete [] flow1;
    delete [] flow2;
}


void occ_error(const std::string& filename_of,
               const std::string& filename_gt,
               float *results){
    // Open input images and .flo
    // pd: number of channels
    int w[2], h[2], pd[2];

    cout << filename_of << "\n";
    //Sparse Optical flow forward and backward
    float *occ_results_float = iio_read_image_float_split(filename_of.c_str(), w, h, pd);
    auto *occ_gt_uint8 = iio_read_image_uint8(filename_gt.c_str(), w + 1, h + 1);


    auto occ_results = new int[w[0]*h[0]];
    auto occ_gt = new int[w[0]*h[0]];

    for (int k = 0; k < w[0]*h[0]; k++){
        occ_results[k] = occ_results_float[k];
        occ_gt[k] = occ_gt_uint8[k];
    }

    int pixels_gt = 0;
    int pixels_results = 0;

    for (int k = 0; k < w[0]*h[0]; k++){

        pixels_gt += occ_gt[k];
        pixels_results += occ_results[k];
    }

    cout << "Pixels gt: " << pixels_gt << "\n";
    cout << "Pixels results: " << pixels_results << "\n";


    int nPix = 0;

    float pixelTP = 0;
    float pixelFP = 0;
    float pixelFN = 0;
    float pixelTN = 0;

    for (int j = 0; j < h[0]; j++){
        for (int i = 0; i < w[0]; i++){
            nPix++;
            int p = i*h[0] + j;

            if (occ_results[p] == 1){
                if (occ_gt[p] == 1){
                    pixelTP++;
                }else{
                    pixelFP++;
                }
            }else{
                if (occ_gt[p] == 1){
                    pixelFN++;
                }else{
                    pixelTN++;
                }
            }
        }
    }

    cout << "pixelTP: " << pixelTP << "\n";
    cout << "pixelFP: " << pixelFP << "\n";
    cout << "pixelFN: " << pixelFN << "\n";
    cout << "pixelTN: " << pixelTN << "\n\n";


    float pixelPrecision = 0.0;


    if (pixelTP == 0.0 && pixelFP == 0.0){
        cout << "Precision changed to 0\n";
        pixelPrecision = 0.0;
    }else{
        pixelPrecision = pixelTP / (pixelTP + pixelFP);
    }

    float pixelSensitivity = 0.0;

    if (pixelTP == 0.0 && pixelFN == 0.0){
        cout << "Recall changed to 0\n";
        pixelSensitivity = 0.0;
    }else{
        pixelSensitivity = pixelTP / (pixelTP + pixelFN);
    }

    float pixelFMeasure = 0.0;

    if (pixelPrecision == 0.0 && pixelSensitivity == 0.0){
        cout << "Fmeasure hanged to 0\n";
        pixelFMeasure = 0.0;
    }else{
        pixelFMeasure = (2*pixelPrecision*pixelSensitivity)/(pixelPrecision + pixelSensitivity);
    }
    cout << "num pixels: " << nPix << "\n";

    std::cout << "Precision: "<< pixelPrecision << "\n";
    std::cout << "Recall: "<< pixelSensitivity << "\n";
    std::cout << "F-Measure: "<< pixelFMeasure << "\n\n";

    if (nPix == 0){
        results[0] = 0;
        results[1] = 0;
        results[2] = 0;
    }else{
        results[0] = pixelPrecision;
        results[1] = pixelSensitivity;
        results[2] = pixelFMeasure;
    }
    delete [] occ_gt;
    delete [] occ_results;
    delete [] occ_results_float;
    delete [] occ_gt_uint8;
    delete [] occ_results;
    delete [] occ_gt;
}

int main(int argc, char* argv[]){

    std::vector<std::string> args(argv, argv + argc);

    if (args.size() <= 1) {
        fprintf(stderr, "Algorithm needs a method and arguments, possible methods: \n"
                        "epe: computes EPE and AE for all flow fields (input is .flo)\n"
                        "epe_match: computes EPE-matched and AE-matched for all flow fields (input is .flo)\n"
                        "fmeasure: conputes precision, recall and Fmeasure for all occlusion masks (input is .png)\n");
        return 1;
    }

    const string& method = args[1];
    string filename_mask_occ = "";
    string method_output = "";
    float inf_limit = 0.0;
    float sup_limit = 0.0;

    if (method == "epe"){
        if (args.size() != 4) {
            fprintf(stderr, "Usage %lu :\n\t%s epe out_flow.txt gt_flow.txt", args.size(), args[0].c_str());

            return 1;
        }
        cout << "Computing EPE and AE\n";
    }else if(method == "epe_match"){
        if (args.size() != 5) {
            fprintf(stderr, "Usage %lu :\n\t%s epe_match out_flow.txt gt_flow.txt mask_gt_occ.txt", args.size(), args[0].c_str());

            return 1;
        }
        cout << "Computing EPE-match and AE-match\n";
        filename_mask_occ = args[4];
    }else if (method == "fmeasure") {
        if (args.size() != 4) {
            fprintf(stderr, "Usage %lu :\n\t%s fmeasure out_occ.txt gt_occ.txt", args.size(), args[0].c_str());

            return 1;
        }
        cout << "Computing precision, recall and F-measure\n";
    }else if (method == "s0_10") {
        if (args.size() != 4) {
            fprintf(stderr, "Usage %lu :\n\t%s s0_10 out_flow.txt gt_flow.txt", args.size(), args[0].c_str());

            return 1;
        }
        cout << "Computing s0-10\n";
        inf_limit = 0.0;
        sup_limit = 10.0;
        method_output = "s0-10";
    }else if (method == "s10_40") {
        if (args.size() != 4) {
            fprintf(stderr, "Usage %lu :\n\t%s s10_40 out_flow.txt gt_flow.txt", args.size(), args[0].c_str());

            return 1;
        }
        cout << "Computing s10-40\n";
        inf_limit = 10.0;
        sup_limit = 40.0;
        method_output = "s10-40";
    }else if (method == "s40+") {
        if (args.size() != 4) {
            fprintf(stderr, "Usage %lu :\n\t%s s40+ out_flow.txt gt_flow.txt", args.size(), args[0].c_str());

            return 1;
        }
        cout << "Computing s40+\n";
        inf_limit = 40.0;
        sup_limit = INFINITY;
        method_output = "s40+";
    }else{
        fprintf(stderr, "The first argument is not a known method, possible methods: \n"
                        "epe: computes EPE and AE for all flow fields (input is .flo)\n"
                        "epe_match: computes EPE-matched and AE-matched for all flow fields (input is .flo)\n"
                        "fmeasure: conputes precision, recall and Fmeasure for all occlusion masks (input is .png)\n");
        return 1;
    }



    const string& filename_images_results = args[2];
    ifstream infile_results(filename_images_results);
    string line_results;

    const string& filename_images_gt = args[3];
    ifstream infile_gt(filename_images_gt);
    string line_gt;

    ifstream infile_mask(filename_mask_occ);
    string line_mask;

    float global_ae = 0.0;
    float global_ee = 0.0;
    float global_precision = 0.0;
    float global_recall = 0.0;
    float global_fmeasure = 0.0;
    float num_files = 0.0;



    while(getline(infile_results, line_results)){

        num_files = num_files + 1.0;
        cout << num_files << "\n";
        getline(infile_gt, line_gt);


        const string& filename_of = line_results;
        const string& filename_gt = line_gt;

        if (method == "epe"){
            float *results = new float[2];
            if (num_files == 1){
                of_error(filename_of, filename_gt, results);
            }
            global_ae += results[0];
            global_ee += results[1];
        }else if(method == "epe_match"){

            float *results = new float[2];
            getline(infile_mask, line_mask);
            const string& filename_mask = line_mask;
            of_error(filename_of, filename_gt, filename_mask, results);

            global_ae += results[0];
            global_ee += results[1];
        }else if (method == "fmeasure") {

            float *results = new float[3];
            occ_error(filename_of, filename_gt, results);

            global_precision += results[0];
            global_recall += results[1];
            global_fmeasure += results[2];
        }else if(method == "s0_10" || method == "s10_40" || method == "s40+" ){
            float *results = new float[2];

            of_error(filename_of, filename_gt, results, inf_limit, sup_limit, method_output);

            global_ee += results[0];


        }
    }

    infile_results.close();
    infile_gt.close();

    if (method == "epe"){

        global_ae = global_ae/num_files;
        global_ee = global_ee/num_files;
        std::cout << "Mean angular error: "<< global_ae << "\n";
        std::cout << "Mean endpoint error: "<< global_ee << "\n";

    }else if(method == "epe_match"){

        global_ae = global_ae/num_files;
        global_ee = global_ee/num_files;
        std::cout << "Mean angular error matched: " << global_ae << "\n";
        std::cout << "Mean endpoint error matched: " << global_ee << "\n";

    }else if (method == "fmeasure") {

        global_precision = global_precision/num_files;
        global_recall = global_recall/num_files;
        global_fmeasure = global_fmeasure/num_files;
        std::cout << "Mean precision: " << global_precision << "\n";
        std::cout << "Mean recall: " << global_recall << "\n";
        std::cout << "Mean Fmeasure: " << global_fmeasure << "\n";
    }else if(method == "s0_10" || method == "s10_40" || method == "s40+" ){
        cout << global_ee << "\n";
        global_ae = global_ae/num_files;
        global_ee = global_ee/num_files;
        std::cout << "Mean " << method_output << ": "<< global_ee << "\n";
    }

    return 0;
}


