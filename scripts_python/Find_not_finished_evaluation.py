

import os 
import argparse
import subprocess
import shlex
import math
from PIL import Image
from utils import list_images_dataset


method = 0

dataset = 'sintel'#middlebury sintel

prova = False

images_dataset = []	

if __name__ == '__main__':
	list_images = list_images_dataset(dataset)
	var_m = method
	#Convert file of images into list
	
	#Evaluate all images of a dataset
	
	folder_out = '../Output_error_evaluation/' + dataset + '_method_' + str(method) + '/'
	if not os.path.exists(folder_out):
		os.makedirs(folder_out)
	folder_scripts = "../scripts_evaluation_datasets/" + dataset + '_method_' + str(method) + '/'
	if not os.path.exists(folder_scripts):
		os.makedirs(folder_scripts)
	range_images = range(len(list_images))


	#For all random trials, evaluate all images
	not_done = 0
	
		
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
		if 'sintel' in dataset:
			if 'Sintel_final' in image_txt:
				f_path = '../Results/Sintel_evaluation/' + 'Method_' + str(var_m) + '/Sintel_final/' + sequence + '/'
				
			else:
				f_path = '../Results/Sintel_evaluation/' + 'Method_' + str(var_m) + '/Sintel_clean/' + sequence + '/'
				
		else:
			f_path = '../Results/Middlebury_evaluation/' + 'Method_' + str(var_m) + '/' + sequence + '/'
			
		

		var_flow = f_path + core_name1 + '_exp_var.flo'
				
		if not os.path.isfile(var_flow):
					
			#print var_flow 
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
					
			not_done = not_done + 1
			#os.system(command)
	print not_done



