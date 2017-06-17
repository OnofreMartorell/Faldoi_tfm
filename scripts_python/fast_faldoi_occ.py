#! /usr/bin/python
#"""
# Litle script for faldoy to execute the data from sift matches.

#"""
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
parser.add_argument("file_images", help = "File with images")

method = 8
matchings = False
local_of = False
global_of = True


#Energy model
parser.add_argument("-vm", default = str(method),
                    help = "Variational Method "
                    "(tv-l2 coupled: 0, tv-l2 coupled_w: 1, NLTVL1:2, NLTVL1:3...")
#M_TVL1       0
#M_TVL1_W     1
#M_NLTVL1     2 
#M_NLTVL1_W   3 
#M_TVCSAD     4
#M_TVCSAD_W   5
#M_NLTVCSAD   6
#M_NLTVCSAD_W 7
#M_TVL1_OCC   8       

if method == 0:
	method_extension = 'tvl1'
elif method == 8:
	method_extension = 'tvl1_occ'
else:
	method_extension = ''



#Local Wise Minimization 
parser.add_argument("-wr", default = '5',
                    help = "Windows Radio Local patch"
                    "1 -  3x3, 2 - 5x5,...") #(2*r +1) x (2*r+1)
#Global Mininization
parser.add_argument("-warps", default = '5',
                    help = "Number of warps finest scale")
#Threshold for Deep Flow
parser.add_argument("-th", default = '0.45',
                    help = "Threshold to discard outliers from DeepFlow")


args = parser.parse_args()
with open(args.file_images, 'r') as file:
	# read a list of lines into data
	data = file.readlines()
for i in range(len(data)):
	data[i] = data[i][:-1]

sequence = data[0].split('.')[2].split('/')[-2]
core_name1 = data[0].split('.')[2].split('/')[-1]
core_name2 = data[1].split('.')[2].split('/')[-1]
var_m = args.vm
warps = args.warps
windows_radio = args.wr;
threshold = args.th;

#C++ program names
match_comparison = "../build/deepmatching"
sparse_flow = "../build/sparse_flow"
match_propagation = "../build/local_faldoi"
of_var = "../build/global_faldoi"


#Set the main directory that contains all the stuff
root_path = '%s/'%(os.getcwd())
#binary_path = root_path + "bin/"
binary_path = '../build/'
#f_path = root_path + "Results/"
f_path = '../Results/' + sequence + '/'
if not os.path.exists(f_path):
    os.makedirs(f_path)

#Set the images input.
im_name0 = os.path.abspath(data[0])
im_name1 = os.path.abspath(data[1])
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

region_growing = '%s%s_rg_%s.flo'%(f_path, core_name1, method_extension)
sim_value = '%s%s_exp_sim_%s.tiff'%(f_path, core_name1, method_extension)
var_flow = '%s%s_exp_var_%s.flo'%(f_path, core_name1, method_extension)

occlusions_rg = '%s%s_rg_occ_%s.png'%(f_path, core_name1, method_extension)
occlusions_var = '%s%s_var_occ%s.png'%(f_path, core_name1, method_extension)

#Obtain the matches' list for both (I0-I1 and I1-I0)
print('Obtaining list of matches from DeepMatching')
max_scale = math.sqrt(2)

#I0-I1
param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name0, im_name1, max_scale, match_name_1)
command_line = '%s %s\n'%(match_comparison, param)
if matchings:
	os.system(command_line)

#I1-I0
param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name1, im_name0, max_scale, match_name_2)
command_line = '%s %s\n'%(match_comparison, param)
if matchings:
	os.system(command_line)


#Create a sparse flow from the deepmatching matches.
print('Creating sparse from matches')

param = '%s %s %s %s\n'%(cut(delete(confi(im_name0, im_name1, match_name_1, f_path), threshold)), width_im, height_im, sparse_name_1)
command_line = '%s %s\n'%(sparse_flow, param)
os.system(command_line)

param = '%s %s %s %s\n'%(cut(delete(confi(im_name1, im_name0, match_name_2, f_path), threshold)), width_im, height_im, sparse_name_2)
command_line = '%s %s\n'%(sparse_flow, param)
os.system(command_line)


#Create a dense flow from a sparse set of initial seeds

options = '-m %s -wr %s'%(var_m, windows_radio)
param = '%s %s %s %s %s %s %s\n'%(args.file_images, sparse_name_1, sparse_name_2, 
                            region_growing, sim_value, occlusions_rg, options)
#print param
command_line = '%s %s\n'%(match_propagation, param)

print(command_line)
if local_of:
	print('Computing local faldoi')
	os.system(command_line)


#Put the dense flow as input for a variational method

# Tv-l2 coupled 0 Du 1
options = '-m %s -w %s'%(var_m, warps)
param = '%s %s %s %s %s %s\n'%(args.file_images,
                            region_growing, var_flow, occlusions_rg, occlusions_var, options)
command_line = '%s %s\n'%(of_var, param)
print(command_line)
if global_of:
	print('Computing global faldoi')
	os.system(command_line)
