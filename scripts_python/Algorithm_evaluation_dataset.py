

import os 
import argparse
import subprocess
import shlex
import math
from os import listdir
from os.path import isfile, join
from PIL import Image
from utils import list_images_dataset


method = 0

dataset = 'sintel'#middlebury sintel

prova = False

images_dataset = []	



if __name__ == '__main__':
	
	list_images = list_images_dataset(dataset)
	
	#Convert file of images into list
	
	#Evaluate all images of a dataset
	
	folder_out = '../Output_error_evaluation/' + dataset + '_method_' + str(method) + '/'
	if not os.path.exists(folder_out):
		os.makedirs(folder_out)
	folder_scripts = "../scripts_evaluation_datasets/" + dataset + '_method_' + str(method) + '/'
	if not os.path.exists(folder_scripts):
		os.makedirs(folder_scripts)
	if prova:
		range_images = [0]
	else:
		range_images = range(len(list_images))

	for i in range_images:
		image = list_images[i]
		#Create sh
			
		file_sh = folder_scripts + "Params_ev_im_" + str(i) + ".sh"
		with open(file_sh, 'w') as file:
			file.write("#!/bin/bash\n")
			cmd = 'Faldoi_evaluation_dataset.py ' + image + '  ' + dataset + ' -vm ' + str(method)

			f = 'time python ' + cmd + '  \n'
			file.write(f)

			#Create sub
		file_sub = folder_scripts + "Params_ev_im_" + str(i) + ".sub"
		with open(file_sub, 'w') as file:
			file.write("#!/bin/bash\n")
			file.write('#$ -N Params_ev_im_' + str(i) + '\n')
			file.write('#$ -cwd\n')
			file.write('#$ -o ' + folder_out + '\n')
			file.write('#$ -e ' + folder_out + '\n')
			file.write('#$ -pe smp 3' + '\n')
			file.write('#$ -q default.q \n')
			file.write('sh '+ file_sh + '\n')
				

		#Send job
		command = "qsub " + file_sub
		print command
		
		os.system(command)
	



