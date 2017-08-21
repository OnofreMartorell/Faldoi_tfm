#! /usr/bin/python
#"""
# Litle script for faldoy to execute the data from sift matches.

#"""

import argparse
import os
from os import listdir
from os.path import isfile, join
import subprocess
import shlex
import math
from rescore_prunning import  confidence_values as confi
from auxiliar_faldoi_functions import cut_deep_list as cut
from auxiliar_faldoi_functions import execute_shell_program as exe_prog
from auxiliar_faldoi_functions import delete_outliers as delete


max_scale = math.sqrt(2)

type_sequence = ['Sintel_final/', 'Sintel_clean/']
directories_images = ['../data/Sintel_final', '../data/Sintel_clean']


#Set the arguments to compute the images
parser = argparse.ArgumentParser(description = 'Deep_matchings')


#Threshold for Deep Flow
parser.add_argument("-th", default = '0.45',
                    help = "Threshold to discard outliers from DeepFlow")


args = parser.parse_args()
threshold = args.th;



#C++ program names
match_comparison = "../build/deepmatching"


#Set the main directory that contains all the stuff
root_path = '%s/'%(os.getcwd())

binary_path = '../build/'

os.chdir(binary_path)


#Obtain the matches' list for both (I0-I1 and I1-I0)


for idx_directory in range(len(directories_images)):
	directory = directories_images[idx_directory]
	folders_sequences = listdir(directory)
	for i in range(len(folders_sequences)):
		
		directory_sequence = join(directory, folders_sequences[i])

		images = sorted([f for f in listdir(directory_sequence) if isfile(join(directory_sequence, f)) and f.endswith(".png")])

		images_dir = [join(directory_sequence, f) for f in images]
		#print images_dir
		for j in range(len(images_dir) - 1):
			data = [images_dir[j], images_dir[j + 1]]
			sequence = data[0].split('.')[2].split('/')[-2]

			im_name0 = os.path.abspath(data[0])
			im_name1 = os.path.abspath(data[1])

			f_path = '../Results/' + type_sequence[idx_directory] + sequence + '/'
			if not os.path.exists(f_path):
    				os.makedirs(f_path)

			core_name1 = data[0].split('.')[2].split('/')[-1]
			core_name2 = data[1].split('.')[2].split('/')[-1]
			


			match_name_1 = '%s%s_exp_mt_1.txt'%(f_path, core_name1)
			match_name_2 = '%s%s_exp_mt_2.txt'%(f_path, core_name2)
			
			#I0-I1
			param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name0, im_name1, max_scale, match_name_1)
			command_line = '%s %s\n'%(match_comparison, param)

			
			print command_line
			#os.system(command_line)

			#I1-I0
			param = '%s %s -downscale 1 -max_scale %s -rot_range -45 +45 > %s'%(im_name1, im_name0, max_scale, match_name_2)
			command_line = '%s %s\n'%(match_comparison, param)
			
			print command_line
			#os.system(command_line)

			
























