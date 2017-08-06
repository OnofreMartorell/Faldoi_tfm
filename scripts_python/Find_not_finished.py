

import os 
import argparse
import subprocess
import shlex
import math
from PIL import Image
random_trials = 6
init_num = 44

if __name__ == '__main__':
	
	file_images_random = '../scripts_python/Random_images.txt'
	#Convert file of images into list
	with open(file_images_random, 'r') as file:
		# read a list of lines into data
		list_images = file.readlines()
	for i in range(len(list_images)):
		list_images[i] = list_images[i][:-1]
	#For all random trials, evaluate all images
	not_done = 0
	for count in range(init_num, init_num + random_trials):

		folder_out = '../Output_error_optimization/params_optimization_trial_' + str(count)
		if not os.path.exists(folder_out):
			os.makedirs(folder_out)
		folder_lists = '../scripts_optimization/lists_images/'
		if not os.path.exists(folder_lists):
			os.makedirs(folder_lists)
		if True:
			for i in range(len(list_images)):
				image_txt = list_images[i]
				with open(image_txt, 'r') as file:
		
					list_images_png = file.readlines()
				for j in range(len(list_images_png)):
					list_images_png[j] = list_images_png[j][:-1]
				#Out flow file
				sequence = list_images_png[0].split('.')[2].split('/')[-2]
				core_name1 = list_images_png[0].split('.')[2].split('/')[-1]
				if 'Sintel_final' in list_images_png[0]:
					subset = 'Sintel_final/'
					
				else:
					subset = 'Sintel_clean/'

				f_path = '../Results/Experiment_' + str(count) + '/' + subset + sequence + '/'
				var_flow = '%s%s_exp_var.flo'%(f_path, core_name1)
				
				if not os.path.isfile(var_flow):
					
					#print var_flow 
					image = list_images[i]
					#Create sh
					file_sh = "../scripts_optimization/Params_opti_" + str(count) + "_im_" + str(i) + ".sh"
					with open(file_sh, 'w') as file:
						file.write("#!/bin/bash\n")
						
						cmd = 'Faldoi_optimization.py ' + image + ' ' + str(count)
		
						f = 'time python ' + cmd + '  \n'
						#print f
						file.write(f)

					#Create sub
					file_sub = "../scripts_optimization/Params_opti_" + str(count) + "_im_" + str(i) + ".sub"
					with open(file_sub, 'w') as file:
						file.write("#!/bin/bash\n")
						file.write('#$ -N Params_optimization_' + str(count) + "_im_" + str(i) + '\n')
						file.write('#$ -cwd\n')
						file.write('#$ -o ' + folder_out + '\n')
						file.write('#$ -e ' + folder_out + '\n')
						file.write('#$ -pe smp 2\n')
						file.write('sh '+ file_sh + '\n')
				

					#Send job
					command = "qsub " + file_sub
					print command
					
					not_done = not_done + 1
					#os.system(command)
	print not_done








