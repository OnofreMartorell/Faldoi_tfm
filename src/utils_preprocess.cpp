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

    params.tol_OF = PAR_DEFAULT_TOL_D;
    params.verbose = PAR_DEFAULT_VERBOSE;
    params.step_algorithm = step_alg;

    if (file_params == ""){
        params.lambda = PAR_DEFAULT_LAMBDA;
        params.beta = PAR_DEFAULT_BETA;
        params.theta = PAR_DEFAULT_THETA;
        params.tau = PAR_DEFAULT_TAU;
        params.alpha = PAR_DEFAULT_ALPHA;
        params.tau_u = PAR_DEFAULT_TAU_U;
        params.tau_eta = PAR_DEFAULT_TAU_ETA;
        params.tau_chi = PAR_DEFAULT_TAU_CHI;
        params.mu = PAR_DEFAULT_MU;

    }else{
        std::string::size_type sz;
        std::string line;
        std::ifstream infile;
        infile.open(file_params);
        getline(infile, line);

        params.lambda = std::stof(line, &sz); getline(infile, line);
        if (params.lambda <= 0) {
            params.lambda = PAR_DEFAULT_LAMBDA;
            if (params.verbose) fprintf(stderr, "warning: lambda changed to %g\n", params.lambda);
        }
        params.theta = std::stof(line, &sz); getline(infile, line);
        if (params.theta <= 0) {
            params.theta = PAR_DEFAULT_THETA;
            if (params.verbose) fprintf(stderr, "warning: "
                                                "theta changed to %g\n", params.theta);
        }
        params.tau = std::stof(line, &sz); getline(infile, line);
        if (params.tau <= 0 || params.tau > 0.25) {
            params.tau = PAR_DEFAULT_TAU;
            if (params.verbose) fprintf(stderr, "warning: "
                                                "tau changed to %g\n", params.tau);
        }
        params.beta = std::stof(line, &sz); getline(infile, line);
        if (params.beta <= 0) {
            params.beta = PAR_DEFAULT_BETA;
            if (params.verbose) fprintf(stderr, "warning: "
                                         "beta changed to %f\n", params.beta);
        }
        params.alpha = std::stof(line, &sz); getline(infile, line);
        if (params.alpha <= 0) {
            params.alpha = PAR_DEFAULT_ALPHA;
            if (params.verbose) fprintf(stderr, "warning: "
                                         "alpha changed to %f\n", params.alpha);
        }
        params.tau_u = std::stof(line, &sz); getline(infile, line);
        if (params.tau_u <= 0 || params.tau_u > 0.25) {
            params.tau_u = PAR_DEFAULT_TAU_U;
            if (params.verbose) fprintf(stderr, "warning: "
                                                "tau_u changed to %g\n", params.tau_u);
        }
        params.tau_eta = std::stof(line, &sz); getline(infile, line);
        if (params.tau_eta <= 0 || params.tau_eta > 0.25) {
            params.tau_eta = PAR_DEFAULT_TAU_ETA;
            if (params.verbose) fprintf(stderr, "warning: "
                                                "tau_eta changed to %g\n", params.tau_eta);
        }
        params.tau_chi = std::stof(line, &sz); getline(infile, line);
        if (params.tau_chi <= 0 || params.tau_chi > 0.25) {
            params.tau_chi = PAR_DEFAULT_TAU_CHI;
            if (params.verbose) fprintf(stderr, "warning: "
                                                "tau_chi changed to %g\n", params.tau_chi);
        }
        params.mu = std::stof(line, &sz);
        if (params.mu <= 0) {
            params.mu = PAR_DEFAULT_MU;
            if (params.verbose) fprintf(stderr, "warning: "
                                         "mu changed to %f\n", params.mu);
        }
        infile.close();
    }
    return params;
}


