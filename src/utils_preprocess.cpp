#include "string"
#include "parameters.h"
#include "energy_structures.h"
#include <fstream>

Parameters init_params(const std::string& file_params, int step_alg){
    Parameters params;

    if (step_alg == LOCAL_STEP){
        params.warps = PAR_DEFAULT_NWARPS_LOCAL;
        params.iterations_of = MAX_ITERATIONS_LOCAL;
    }else{
        params.warps = PAR_DEFAULT_NWARPS_GLOBAL;
        params.iterations_of = MAX_ITERATIONS_GLOBAL;
    }
    if (file_params == ""){
        params.lambda = PAR_DEFAULT_LAMBDA;
        params.beta = PAR_DEFAULT_BETA;
        params.theta = PAR_DEFAULT_THETA;
        params.tau = PAR_DEFAULT_TAU;
        params.alpha = PAR_DEFAULT_ALPHA;
        params.tau_u = PAR_DEFAULT_TAU_U;
        params.tau_eta = PAR_DEFAULT_TAU_ETA;
        params.tau_chi = PAR_DEFAULT_TAU_CHI;

    }else{
        std::string::size_type sz;
        std::string line;
        std::ifstream infile;
        infile.open(file_params);
        getline(infile, line);

        params.lambda = std::stof(line, &sz); getline(infile, line);
        params.theta = std::stof(line, &sz); getline(infile, line);
        params.tau = std::stof(line, &sz); getline(infile, line);
        params.beta = std::stof(line, &sz); getline(infile, line);
        params.alpha = std::stof(line, &sz); getline(infile, line);
        params.tau_u = std::stof(line, &sz); getline(infile, line);
        params.tau_eta = std::stof(line, &sz); getline(infile, line);
        params.tau_chi = std::stof(line, &sz); getline(infile, line);

        infile.close();
    }
    params.tol_OF = PAR_DEFAULT_TOL_D;
    params.verbose = PAR_DEFAULT_VERBOSE;
    params.step_algorithm = step_alg;
    return params;
}
