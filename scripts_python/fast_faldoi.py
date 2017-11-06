#! /usr/bin/python
"""
 Litle script for faldoy to execute the data from Deepmatching matches.

"""
import argparse
import os
import subprocess
import shlex
import math
from rescore_prunning import  confidence_values as confi
from auxiliar_faldoi_functions import cut_deep_list as cut
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from auxiliar_faldoi_functions import delete_outliers as delete
#Set the arguments to compute the images
parser = argparse.ArgumentParser(description = 'Faldoy Minimization')
parser.add_argument("i0", help = "first image")
parser.add_argument("i1", help = "second image")

#Energy model
parser.add_argument("-vm", default = '0',
                    help = "Variational Method "
                    "(tv-l2 coupled: 0, tv-l2 coupled_w: 1, NLTVL1:2, NLTVL1:3...")
#define M_TVL1       0
#define M_TVL1_W     1
#define M_NLTVL1     2 
#define M_NLTVL1_W   3 
#define M_TVCSAD     4
#define M_TVCSAD_W   5
#define M_NLTVCSAD   6
#define M_NLTVCSAD_W 7

#Local Wise Minimization 
parser.add_argument("-wr", default = '5',
                    help = "Windows Radio Local patch"
                    "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)
#Global Mininization
parser.add_argument("-warps", default = '5',
                    help="Number of warps finest scale")
#Threshold for Deep Flow
parser.add_argument("-th", default = '0.45',
                    help = "Threshold to discard outliers from DeepFlow")


args = parser.parse_args()
core_name1 = args.i0.split('.')[2].split('/')[-1]
core_name2 = args.i1.split('.')[2].split('/')[-1]

var_m = args.vm
warps = args.warps
windows_radio = args.wr;
threshold = args.th;

#C program names
match_comparison = "../build/deepmatching"
sparse_flow = "../build/sparse_flow"
match_propagation = "../build/local_faldoi"
of_var = "../build/global_faldoi"


#Set the main directory that contains all the stuff
root_path = '%s/'%(os.getcwd())
binary_path = '../build/'
f_path = '../Results/'
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
match_name_1 = '%s%s_exp_mt_1.txt'%(f_path, core_name1)
sparse_name_1 = '%s%s_exp_mt_1.flo'%(f_path, core_name1)
match_name_2 = '%s%s_exp_mt_2.txt'%(f_path, core_name2)
sparse_name_2 = '%s%s_exp_mt_2.flo'%(f_path, core_name2)
region_growing = '%s%s_rg.flo'%(f_path, core_name1)
sim_value = '%s%s_exp_sim.tiff'%(f_path, core_name1)
var_flow = '%s%s_exp_var.flo'%(f_path, core_name1)
#Obtain the matches' list for both (I0-I1 and I1-I0)
print('Obtaining list of matches from DeepMatching')
max_scale = math.sqrt(2)
#I0-I1
#param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name1, im_name2, max_scale, match_name_1)
#command_line = '%s %s\n'%(match_comparison, param)
#os.system(command_line)
#I1-I0
#param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name2, im_name1, max_scale, match_name_2)
#command_line = '%s %s\n'%(match_comparison, param)
#os.system(command_line)
#Create a sparse flow from the deepmatching matches.
print('Creating sparse from matches')
param = '%s %s %s %s\n'%(cut(delete(confi(im_name1, im_name2, match_name_1, f_path), threshold)), width_im, height_im, sparse_name_1)
command_line = '%s %s\n'%(sparse_flow, param)
os.system(command_line)
param = '%s %s %s %s\n'%(cut(delete(confi(im_name2, im_name1, match_name_2, f_path), threshold)), width_im, height_im, sparse_name_2)
command_line = '%s %s\n'%(sparse_flow, param)
os.system(command_line)
#Create a dense flow from a sparse set of initial seeds
print('Computing local faldoi')
options = '-m %s -wr %s'%(var_m, windows_radio)
param = '%s %s %s %s %s %s %s\n'%(im_name1, im_name2, sparse_name_1,sparse_name_2, 
                            region_growing, sim_value, options)
#print param
command_line = '%s %s\n'%(match_propagation, param)
os.system(command_line)
#Put the dense flow as input for a variational method
print('Computing global faldoi')
# Tv-l2 coupled 0 Du 1
options = '-m %s -w %s'%(var_m, warps)
param = '%s %s %s %s %s\n'%(im_name1, im_name2,
                            region_growing, var_flow, options)
command_line = '%s %s\n'%(of_var, param)
os.system(command_line)
