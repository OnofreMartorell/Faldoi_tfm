// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017, Onofre Martorell <onofremartorelln@gmail.com>
// All rights reserved.

#ifndef FALDOI
#define FALDOI


#include <stdlib.h>
#include "parameters.h"
#include "utils_preprocess.h"

int main(int argc, char *argv[]) {



    int method = M_TVL1_OCC;

    //#M_TVL1       0
    //#M_TVL1_W     1
    //#M_NLTVL1     2
    //#M_NLTVL1_W   3
    //#M_TVCSAD     4
    //#M_TVCSAD_W   5
    //#M_NLTVCSAD   6
    //#M_NLTVCSAD_W 7
    //#M_TVL1_OCC   8

    bool sparse_flow = True;
    bool local_of = True;
    bool global_of = True;

    using namespace std::chrono;
    system_clock::time_point today = system_clock::now();
    time_t tt;

    tt = system_clock::to_time_t ( today );
    std::cerr << "today is: " << ctime(&tt);

    // process input
    std::vector<std::string> args(argv, argv + argc);

    //#Local Wise Minimization
    //parser.add_argument("-wr", default = '5',
    //                    help = "Windows Radio Local patch"
    //                    "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)
    //#Global Mininization
    //parser.add_argument("-warps", default = '7',
    //                    help = "Number of warps finest scale")
    //#Energy model
    //parser.add_argument("-vm", default = str(method),
    //                    help = "Variational Method "
    //                    "(tv-l2 coupled: 0, tv-l2 coupled_w: 1, NLTVL1:2, NLTVL1:3...")
    //#Ground truth for evaluation
    //parser.add_argument("-ev", default = '',
    //                    help = "File with ground truth optical")
    //#Ground truth for evaluation
    //parser.add_argument("-p", default = '',
    //                    help = "File with parameters")

    auto windows_radio  = pick_option(args, "wr", "5");// Windows radius for local minimization
    auto warps_val      = pick_option(args, "warps", to_string(PAR_DEFAULT_NWARPS_GLOBAL));// Warpings global
    auto method_ev      = pick_option(args, "vm", std::to_string(method));

    auto file_gt        = pick_option(args, "ev",  ""); // Methods
    auto file_params    = pick_option(args, "p",  ""); // File of parameters

    if (args.size() != 1) {
        fprintf(stderr, "Usage: %lu  ims.txt in_flow.flo  out.flo [-m] val [-w] val \n", args.size());

        return EXIT_FAILURE;
    }

    string method_extension = "";
    if (method == M_TVL1){
            method_extension = "tvl1";
    }else{
        if (method == M_TVL1_OCC){
            method_extension = "tvl1_occ";
    }


    //with open(args.file_images, 'r') as file:
    //        # read a list of lines into data
    //        data = file.readlines()
    //for i in range(len(data)):
    //        data[i] = data[i][:-1]

    //sequence = data[0].split('.')[2].split('/')[-2]
    //core_name1 = data[0].split('.')[2].split('/')[-1]
    //core_name2 = data[1].split('.')[2].split('/')[-1]
    //var_m = args.vm
    //warps = args.warps
    //windows_radio = args.wr;
    //threshold = args.th;

    //#C++ program names

    std::string sparse_flow = "../build/sparse_flow";
    std::string match_propagation = "../build/local_faldoi";
    std::string of_var = "../build/global_faldoi";
    std::string evaluation = "../build/evaluation";


    //#Set the main directory that contains all the stuff
    //root_path = '%s/'%(os.getcwd())
    //#binary_path = root_path + "bin/"
    //binary_path = '../build/'
    //#f_path = root_path + "Results/"
    //f_path = '../Results/' + sequence + '/'
    //if not os.path.exists(f_path):
    //    os.makedirs(f_path)
    //filename_gt = args.ev
    //filename_params = args.p
    //if not filename_params == '':
    //        iteration_params = '_' + str(filename_params.split('/')[-1].split('.')[0].split('_')[1]).zfill(2)
    //else:
    //        iteration_params = ''


    //#Set the images input.
    //im_name0 = os.path.abspath(data[0])
    //im_name1 = os.path.abspath(data[1])
    //#Get the image size
    //from PIL import Image
    //with open(im_name1, 'r') as f:
    //    image = Image.open(f)
    //    width_im = image.size[0]
    //    height_im = image.size[1]

    //os.chdir(binary_path)

    //match_name_1 = '%s%s_exp_mt_1.txt'%(f_path, core_name1)
    //sparse_name_1 = '%s%s_exp_mt_1.flo'%(f_path, core_name1)

    //match_name_2 = '%s%s_exp_mt_2.txt'%(f_path, core_name2)
    //sparse_name_2 = '%s%s_exp_mt_2.flo'%(f_path, core_name2)

    //region_growing = '%s%s_rg_%s%s.flo'%(f_path, core_name1, method_extension, iteration_params)
    //sim_value = '%s%s_exp_sim_%s%s.tiff'%(f_path, core_name1, method_extension, iteration_params)
    //var_flow = '%s%s_exp_var_%s%s.flo'%(f_path, core_name1, method_extension, iteration_params)

    //occlusions_rg = '%s%s_rg_occ_%s%s.png'%(f_path, core_name1, method_extension, iteration_params)
    //occlusions_var = '%s%s_var_occ_%s%s.png'%(f_path, core_name1, method_extension, iteration_params)


    //#Create a sparse flow from the deepmatching matches.
    //print('Creating sparse from matches')

    //param = '%s %s %s %s\n'%(cut(delete(confi(im_name0, im_name1, match_name_1, f_path), threshold)), width_im, height_im, sparse_name_1)
    //command_line = '%s %s\n'%(sparse_flow, param)
    //if (sparse_flow){
    //        os.system(command_line)
    //}

    //param = '%s %s %s %s\n'%(cut(delete(confi(im_name1, im_name0, match_name_2, f_path), threshold)), width_im, height_im, sparse_name_2)
    //command_line = '%s %s\n'%(sparse_flow, param)

    //if (sparse_flow){
    //        os.system(command_line)
    //}

    //#Create a dense flow from a sparse set of initial seeds
    //if (!(filename_params == '')){
    //        options = '-m %s -wr %s -p %s'%(var_m, windows_radio, filename_params)
    //}else{
    //        options = '-m %s -wr %s'%(var_m, windows_radio)
    //}
    //param = '%s %s %s %s %s %s %s\n'%(args.file_images, sparse_name_1, sparse_name_2,
    //                            region_growing, sim_value, occlusions_rg, options)
    //#print param
    //command_line = '%s %s\n'%(match_propagation, param)

    //#print(command_line)
    //if local_of:
    //        print('Computing local faldoi')
    //        os.system(command_line)


    //#Put the dense flow as input for a variational method

    //# Tv-l2 coupled 0 Du 1
    //if not filename_params == '':
    //        options = '-m %s -w %s -p %s'%(var_m, warps, filename_params)
    //else:
    //        options = '-m %s -w %s'%(var_m, warps)
    //param = '%s %s %s %s %s %s\n'%(args.file_images,
    //                            region_growing, var_flow, occlusions_rg, occlusions_var, options)
    //command_line = '%s %s\n'%(of_var, param)
    //#print(command_line)
    //if global_of:
    //        print('Computing global faldoi')
    //        os.system(command_line)

    //#Evaluate results of method
    //    if  not filename_gt == '':
    //        command_line = evaluation + ' ' + var_flow + ' ' + filename_gt
    //    #print(command_line)
    //        os.system(command_line)
    //      today = system_clock::now();

    tt = system_clock::to_time_t ( today );
    std::cerr << "today is: " << ctime(&tt);
    return EXIT_SUCCESS;
}

#endif