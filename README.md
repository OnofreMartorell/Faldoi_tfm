===============================================================================
== Compiling ==========================================================

To compile

Type

    make

in the directory where the Makefile is located. The compilation of the source
 code provides three executables:

1) sparse_flow - Convert the .txt files from  match_cli to a sparse flow

2) local_faldoi - Estimates an optical flow from a set of initial sparse 
                    matches.

3) global_faldoi - Estimates a global minimization using as initialization the 
                    output from local_faldoi.

===============================================================================
==  Executable - Usage ====================================================


Usage: ./sparse_flow sift_matches.txt colum row out.flo

Usage: ./local_faldoi i0 i1 in.flo out.flo sim_map.tiff [options...]
  options:
            -m (0)  Changes the functional (check aux_energy_model.h)
            -wr (5) Radius value 5 - patch 11x11 

Usage: ./global_faldoi I0 I1 input.flo out.flow
  options: 
            -m (0)  Changes the functional (check aux_energy_model.h)
            -w (5)  Number of warpings   
===============================================================================
==  Script - Usage ====================================================
To use the python script (fast_faldoi.py) that makes all the process:
    - Create the folders bin/ and Results/
    - Compile the executables through Makefile
    - Move the excutable files to bin/
    - Put the excutable deep matching from DeepMatching: Deep Convolutional Matching into 
bin (http://lear.inrialpes.fr/src/deepmatching/)

Usage: ./fast_faldoi.py i0 i1

==  Script - Usage ====================================================
To use the python script (fast_sift.py) that makes all the process:
    - Create the folders bin/ and Results/
    - Compile the executables through Makefile
    - Move the excutable files to bin/
    - Put the excutables sift_cli and match_cli from Anatomy of SIFT from IPOL into 
bin (http://www.ipol.im/pub/art/2014/82/)

Usage: ./fast_sift.py i0 i1
