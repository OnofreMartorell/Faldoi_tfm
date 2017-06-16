


void  of_errors(flow1, flow2, w, h, pd){

    //[ni, nj, nChannels] = size(flow1);
    image_ae = new float[w*h];
    image_ee = new float [w*h];
    float ae = 0;
    float ee = 0;
    int nPix = 0;
    for (int j = 1; j <= h; j++){
        for (int i = 1; i <= w; i++){

            //Exactly as in C code

            u_c = flow1[];//(i,j,1);
            v_c = flow1[];//(i,j,2);

            u_e = flow2[];//(i,j,1);
            v_e = flow2[];//(i,j,2);

            if ((abs(u_c) < 10000) && (abs(v_c) < 10000)){
                nPix = nPix + 1;
                esc_prod = u_c*u_e + v_c*v_e + 1.0;

                term1 = (u_c * u_c + v_c * v_c + 1.0);
                term2 = (u_e * u_e + v_e * v_e + 1.0);
                term3 = sqrt(term1 * term2);

                tmp = esc_prod / term3;


                tmp = min(tmp, 1.0);
                tmp = max(tmp, -1.0);

                ae = ae + acos(tmp);
                image_ae(i, j) = acos(tmp);
                ee = ee + sqrt((u_c - u_e)^2 + (v_c - v_e)^2);
                image_ee(i, j) = sqrt((u_c - u_e)^2 + (v_c - v_e)^2);
            }
        }
    }
    ae = ae/nPix;
    ee = ee/nPix;
}


int main(int argc, char* argv[]){

    std::vector<std::string> args(argv, argv + argc);
    if (args.size() != 2) {
        fprintf(stderr, "usage %d :\n\t%s out.flo gt.flo", args.size(), args[0].c_str());

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
    return 0;

}
