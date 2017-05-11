#! /usr/bin/python
"""
 Litle script for faldoy to execute the data from sift matches.

"""
import argparse
import os
import subprocess
import shlex
from auxiliar_faldoi_functions import cut_matching_list as cut
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from auxiliar_faldoi_functions import delete_outliers as delete
#Set the arguments to compute the images
parser = argparse.ArgumentParser(description='Faldoy Minimization')
parser.add_argument("i0", help="first image")
parser.add_argument("i1", help="second image")
#Energy model
parser.add_argument("-vm", default='0',
                    help="Variational Method "
                    "(tv-l2 coupled: 0, ||Du+Du'||: 1, NLTVL1:3")
#Local Wise Minimization 
parser.add_argument("-wr", default='5',
                    help="Windows Radio Local patch"
                    "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)
#Global Mininization
parser.add_argument("-warps",default='5',
                    help="Number of warps finest scale")
#Initial seeds (SIFT parameters)
parser.add_argument("-nsp", default='15',
                    help="Increase the sift matches")

parser.add_argument("-m",default ='0',
                    help="It uses the Gaussian weight over the Data Term");

args = parser.parse_args()
core_name1 = args.i0.split('.')[0].split('/')[-1]
core_name2 = args.i1.split('.')[0].split('/')[-1]
var_m = args.vm
warps = args.warps
windows_radio = args.wr;
gauss = args.m;
nsp = args.nsp;

param_sif = '-ss_nspo %s'%(nsp)

#C program names
feature_descriptor = "./sift_cli "
match_comparison = "./match_cli"
sparse_flow = "./sparse_flow"
match_propagation = "./local_faldoi"
of_var = "./global_faldoi"


#Set the main directory that contains all the stuff
root_path = '%s/'%(os.getcwd())
binary_path = root_path + "bin/"
f_path = root_path + "Results/"
#Set the folder where the binaries are.
#Set the images input.
im_name1 = os.path.abspath(args.i0)
im_name2 = os.path.abspath(args.i1)
#Get the image size
from PIL import Image
with open(im_name1, 'r') as f:
    image = Image.open(f)
    width_im = image.size[0]
    height_im = image.size[1]
os.chdir(binary_path)
desc_name_1 = '%s%s_sift_desc_1.txt'%(f_path, core_name1)
desc_name_2 = '%s%s_sift_desc_2.txt'%(f_path, core_name2)
match_name_1 = '%s%s_sift_mt_1.txt'%(f_path, core_name1)
sparse_name_1 = '%s%s_sift_mt_1.flo'%(f_path, core_name1)
match_name_2 = '%s%s_sift_mt_2.txt'%(f_path, core_name2)
sparse_name_2 = '%s%s_sift_mt_2.flo'%(f_path, core_name2)
region_growing = '%s%s_sift_rg.flo'%(f_path, core_name1)
sim_value = '%s%s_sift_sim.tiff'%(f_path, core_name1)
var_flow = '%s%s_sift_var.flo'%(f_path, core_name1)
# Obtain the matches' list for both (I0-I1 and I1-I0)
#Initial seeds (SIFT descriptors)
command_line = '%s %s %s\n'%(feature_descriptor, im_name1, param_sif)
print command_line
# exe_prog(command_line, desc_name_1)
command_line = '%s %s %s\n'%(feature_descriptor, im_name2, param_sif)
# exe_prog(command_line, desc_name_2)
#Obtain the matches' list
command_line = '%s %s %s\n'%(match_comparison, desc_name_1, desc_name_2)
exe_prog(command_line, match_name_1)
command_line = '%s %s %s\n'%(match_comparison, desc_name_2, desc_name_1)
# exe_prog(command_line, match_name_2)
#Create a sparse flow from the sift matches.
param = '%s %s %s %s\n'%(cut(match_name_1), width_im, height_im, sparse_name_1)
command_line = '%s %s\n'%(sparse_flow, param)
# os.system(command_line)
#Create a sparse flow from the sift matches.
param = '%s %s %s %s\n'%(cut(match_name_2), width_im, height_im, sparse_name_2)
command_line = '%s %s\n'%(sparse_flow, param)
# os.system(command_line)
#Create a dense flow from a sparse set of initial seeds
options = '-m %s -wr %s'%(var_m, windows_radio)
param = '%s %s %s %s %s %s %s\n'%(im_name1, im_name2, sparse_name_1,sparse_name_2, 
                            region_growing, sim_value, options)
#print param
command_line = '%s %s\n'%(match_propagation, param)
print command_line
os.system(command_line)
# Put the dense flow as input for a variational method
# Tv-l2 coupled 0 Du 1
options = '-m %s -w %s'%(var_m, warps)
param = '%s %s %s %s %s\n'%(im_name1, im_name2,
                            region_growing, var_flow, options)
command_line = '%s %s\n'%(of_var, param)
# os.system(command_line)
