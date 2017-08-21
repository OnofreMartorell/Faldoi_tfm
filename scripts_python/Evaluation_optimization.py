

import os 
import argparse
import subprocess
import shlex
import math
from PIL import Image




random_trials = 44
init_num = 0
				
method_ev = 'fmeasure' #epe_match/fmeasure


if __name__ == '__main__':
	
	file_images_random = '../scripts_python/Random_images.txt'
	#Convert file of images into list
	with open(file_images_random, 'r') as file:
		# read a list of lines into data
		list_images = file.readlines()
	for i in range(len(list_images)):
		list_images[i] = list_images[i][:-1]

	folder_out = '../Output_error_optimization_results/' + method_ev
	if not os.path.exists(folder_out):
		os.makedirs(folder_out)
	folder_lists = '../scripts_optimization/lists_images/'
	if not os.path.exists(folder_lists):
		os.makedirs(folder_lists)

	#For all random trials, evaluate all images
	for count in range(init_num, init_num + random_trials):


		filename_gt = folder_lists + 'List_gt_' + method_ev + '_' + str(count) + '.txt'
		filename_out = folder_lists + 'List_out_' + method_ev + '_' + str(count) + '.txt'
		filename_mask = folder_lists + 'List_mask_' + method_ev + '_' + str(count) + '.txt'

		with open(filename_gt, 'w') as file_gt, open(filename_out, 'w') as file_out, open(filename_mask, 'w') as file_mask:
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

				if method_ev == 'epe' or method_ev == 'epe_match':

					f_path = '../Results/Experiment_' + str(count) + '/' + subset + sequence + '/'
					out = f_path + core_name1 + '_exp_var.flo'
				
					gt = '../Ground_truth/flow_sintel/' + sequence + '/' + core_name1 + '.flo'
					mask = '../Ground_truth/occlusions_sintel/' + sequence + '/' + core_name1 + '.png'

				elif method_ev == 'fmeasure':
					f_path = '../Results/Experiment_' + str(count) + '/' + subset + sequence + '/'
					out = f_path + core_name1 + '_var_occ.png'
					
					gt = '../Ground_truth/occlusions_sintel/' + sequence + '/' + core_name1 + '.png'
					mask = ''
			
				

				file_gt.write(gt + '\n')
				file_out.write(out + '\n')
				file_mask.write(mask + '\n')

		#Create sh
		file_sh = "../scripts_optimization/Evaluation_" + method_ev + '_'  + str(count) + ".sh"
		with open(file_sh, 'w') as file:
			file.write("#!/bin/bash\n")
				
			if method_ev == 'epe':
				cmd = '../build/evaluation epe ' + filename_out + ' ' + filename_gt
			elif method_ev == 'epe_match':
				cmd = '../build/evaluation epe_match ' + filename_out + ' ' + filename_gt + ' ' + filename_mask
			elif method_ev == 'fmeasure':
				cmd = '../build/evaluation fmeasure ' + filename_out + ' ' + filename_gt				
			file.write(cmd)
			
			#os.system(cmd)

		#Create sub
		file_sub = "../scripts_optimization/Evaluation_" + method_ev + '_' + str(count) + ".sub"
		with open(file_sub, 'w') as file:
			file.write("#!/bin/bash\n")
			file.write('#$ -N Params_optimization_results_' + str(count) + '\n')
			file.write('#$ -cwd\n')
			file.write('#$ -o ' + folder_out + '\n')
			file.write('#$ -e ' + folder_out + '\n')
			file.write('sh '+ file_sh + '\n')
				

		#Send job
		command = "qsub " + file_sub
		print command
		os.system(command)








