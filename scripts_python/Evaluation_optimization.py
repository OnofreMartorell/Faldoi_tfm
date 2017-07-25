

import os 
import argparse
import subprocess
import shlex
import math
from PIL import Image
random_trials = 10
init_num = 10

if __name__ == '__main__':
	
	file_images_random = '../scripts_python/Random_images.txt'
	#Convert file of images into list
	with open(file_images_random, 'r') as file:
		# read a list of lines into data
		list_images = file.readlines()
	for i in range(len(list_images)):
		list_images[i] = list_images[i][:-1]
	#For all random trials, evaluate all images
	for count in range(init_num, init_num + random_trials):

		folder_out = '../Output_error_optimization_results/params_optimization_trial_' + str(count)
		if not os.path.exists(folder_out):
			os.makedirs(folder_out)
		folder_lists = '../scripts_optimization/lists_images/'
		if not os.path.exists(folder_lists):
			os.makedirs(folder_lists)
		filename_gt = folder_lists + 'List_flow_gt_' + str(count) + '.txt'
		filename_of = folder_lists + 'List_flow_of_' + str(count) + '.txt'
		with open(filename_gt, 'w') as file_gt, open(filename_of, 'w') as file_of:
			for i in range(len(list_images)):
				image_txt = list_images[i]
				with open(image_txt, 'r') as file:
		
					list_images_png = file.readlines()
				for i in range(len(list_images_png)):
					list_images_png[i] = list_images_png[i][:-1]
				#Out flow file
				sequence = list_images_png[0].split('.')[2].split('/')[-2]
				core_name1 = list_images_png[0].split('.')[2].split('/')[-1]
				f_path = '../Results/Experiment_' + str(count) + '/' + sequence + '/'
				var_flow = '%s%s_exp_var.flo'%(f_path, core_name1)
				
				gt_flow = '../Ground_truth/flow/' + sequence + '/' + core_name1 + '.flo'	
				
				file_gt.write(gt_flow + '\n')
				file_of.write(var_flow + '\n')


		#Create sh
		file_sh = "../scripts_optimization/Evaluation_" + str(count) + ".sh"
		with open(file_sh, 'w') as file:
			file.write("#!/bin/bash\n")
				
			cmd = '../build/evaluation ' + filename_of + ' ' + filename_gt				
			file.write(cmd)

			#Create sub
		file_sub = "../scripts_optimization/Evaluation_" + str(count) + ".sub"
		with open(file_sub, 'w') as file:
			file.write("#!/bin/bash\n")
			file.write('#$ -N Params_optimization_results' + str(count) + '\n')
			file.write('#$ -cwd\n')
			file.write('#$ -o ' + folder_out + '\n')
			file.write('#$ -e ' + folder_out + '\n')
			file.write('sh '+ file_sh + '\n')
				

		#Send job
		command = "qsub " + file_sub
		#print command
		os.system(command)








