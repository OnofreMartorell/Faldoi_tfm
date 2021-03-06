# Faldoi minimization strategy

This respository contains the source code of Faldoi algorithm described in [Faldoi](https://link.springer.com/article/10.1007/s10851-016-0688-y). In addition to all the functionals explained in the article this code has another functional implemented, the one from [Occlusion estimation](https://link.springer.com/chapter/10.1007/978-3-642-32717-9_4). 


## Compiling
To compile type
```
        mkdir build
        cd build
        cmake ../src -DCMAKE_BUILD_TYPE=RELEASE
        make -j4
```

in the directory where this file is located. The compilation of the source
 code provides three executables:

1) sparse_flow - Convert the .txt files from  match_cli to a sparse flow

2) local_faldoi - Estimates an optical flow from a set of initial sparse 
                    matches.

3) global_faldoi - Estimates a global minimization using as initialization the 
                    output from local_faldoi.

## Executable - Usage


Usage: ``./sparse_flow list_matches.txt colums rows out.flo``

Usage: ``./local_faldoi i0 i1 in.flo out.flo sim_map.tiff [options...]``


options:
-            -m (0)  Changes the functional (check aux_energy_model.h)
-            -wr (5) Radius value 5 - patch 11x11 

Usage: ``./global_faldoi I0 I1 input.flo out.flow``


options: 
-            -m (0)  Changes the functional (check aux_energy_model.h)
-            -w (5)  Number of warpings   
## Script - Usage
There exist some python scripts for executing all the code at once. To use them, go to scripts_python folder. In the folder example_data, there is an example of images to use.
### fast_faldoi
To use the python script (fast_faldoi.py) that makes all the process, first you need to do the following:
- Compile the executables through the commands above
- Put the excutable deep matching from DeepMatching: Deep Convolutional Matching into 
build directory from [DeepMatchings](http://lear.inrialpes.fr/src/deepmatching/)

Usage: ``python fast_faldoi.py i0 i1``

### fast_sift
To use the python script (fast_sift.py) that makes all the process, first you need to do the following:
- Compile the executables throughthe commands above
- Put the excutables sift_cli and match_cli from Anatomy of SIFT from IPOL into 
build directory from [SIFT](http://www.ipol.im/pub/art/2014/82/)

Usage: ``python fast_sift.py i0 i1``

### fast_faldoi_occ
To use the python script (fast_faldoi.py) that makes all the process, first you need to do the following:
- Compile the executables through the commands above
- Put the excutable deep matching from DeepMatching: Deep Convolutional Matching into 
build directory from [DeepMatchings](http://lear.inrialpes.fr/src/deepmatching/)
- Create a txt file with the path of the images in the following order: I_0, I_1, I_{-1}, I_2

Usage: ``python fast_faldoi_occ.py list_images.txt``    
