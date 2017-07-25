

import os 
import argparse
import subprocess
import shlex
import math
from PIL import Image

random_trials = 1
init_num = 0


method = 8

#create_sh = "for i in {1..18}; do echo 'time python fast_faldoi_occ.py ../data/Urban2/frames_9101112.txt -ev ../Ground_truth/Urban2/flow10.flo -p ../Parameters_files/params_$i.txt' > opti_urban_$i.sh; done"
#python -c 'import foo; print foo.hello()'
#fast_faldoi_occ(file_images, num_param_file, vm = str(method), wr = '5', warps = '7', th = '0.45', ev = '', p = ''):



#$ -N Params_optimization
#
# (2) We'll redirect the output files to
# our working directory (flags -cwd, -e, -o)
# ---------------------------------------
#$ -cwd
#$ -o ../Output_error/params_optimization
#$ -e ../Output_error/params_optimization
#
# (3)Finally, we call the script
# ------------------------------
#sh Params_opti.sh


prova = True


if __name__ == '__main__':
	two_image_files = {49, 75, 120, 122, 133, 137, 171}
	if prova:
		file_images_random = '../scripts_python/images.txt'
	else:
		file_images_random = '../scripts_python/Random_images.txt'
	#Convert file of images into list
	with open(file_images_random, 'r') as file:
		# read a list of lines into data
		list_images = file.readlines()
	for i in range(len(list_images)):
		list_images[i] = list_images[i][:-1]
	#For all random trials, evaluate all images
	for count in range(init_num, init_num + random_trials):
		folder_out = '../Output_error_optimization/params_optimization_trial_' + str(count)
		if not os.path.exists(folder_out):
			os.makedirs(folder_out)
		for i in range(len(list_images)):
			image = list_images[i]
			#Create sh
			file_sh = "../scripts_optimization/Params_opti_" + str(count) + "_im_" + str(i) + ".sh"
			with open(file_sh, 'w') as file:
				file.write("#!/bin/bash\n")
				if False:
					cmd = 'import Parameter_optimization as param; param.fast()'
				else:
					cmd = 'Faldoi_optimization.py ' + image + ' ' + str(count)

				f = 'time python ' + cmd + '  \n'
				file.write(f)

			#Create sub
			file_sub = "../scripts_optimization/Params_opti_" + str(count) + "_im_" + str(i) + ".sub"
			with open(file_sub, 'w') as file:
				file.write("#!/bin/bash\n")
				file.write('#$ -N Params_optimization_' + str(count) + "_im_" + str(i) + '\n')
				file.write('#$ -cwd\n')
				file.write('#$ -o ' + folder_out + '\n')
				file.write('#$ -e ' + folder_out + '\n')
				if i == len(list_images) - 1:
					
					file.write('#$ -M onofremartorelln@gmail.com\n')
					file.write('#$ -m e\n')
				file.write('sh '+ file_sh + '\n')
				

			#Send job
			command = "qsub " + file_sub
			#print command
			if i in two_image_files:
				print command
				#os.system(command)
	














