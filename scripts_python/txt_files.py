#! /usr/bin/python


from os import listdir
from os.path import isfile, join


directory_images_clean = '../data/Sintel_clean'
folders_sequences = listdir(directory_images_clean)
for i in range(len(folders_sequences)):
	directory_sequence = join(directory_images_clean, folders_sequences[i])
	images = sorted([f for f in listdir(directory_sequence) if isfile(join(directory_sequence, f)) and f.endswith(".png")])
	images_dir = [join(directory_sequence, f) for f in images]
	#print images_dir
	for j in range(2, len(images_dir) - 1):
		nums = str(j).zfill(2) + str(j + 1).zfill(2) + str(j - 1).zfill(2) + str(j + 2).zfill(2)
		
		#print directory_sequence + "/frames_" + nums + ".txt"
		with open(directory_sequence + "/frames_" + nums + ".txt", 'w') as file:
			file.write(str(images_dir[j - 1]) + "\n")
			file.write(str(images_dir[j]) + "\n")	
			file.write(str(images_dir[j - 2]) + "\n")			
			file.write(str(images_dir[j + 1]) + "\n")
	#Firts two images	
	nums = "0102"
	#print directory_sequence + "/frames_" + nums + ".txt"
	with open(directory_sequence + "/frames_" + nums + ".txt", 'w') as file:
			file.write(str(images_dir[0]) + "\n")
			file.write(str(images_dir[1]) + "\n")	
	nums = str(len(images_dir) - 1).zfill(2) + str(len(images_dir)).zfill(2)
	#Last two images
	#print directory_sequence + "/frames_" + nums + ".txt"
	with open(directory_sequence + "/frames_" + nums + ".txt", 'w') as file:
			file.write(str(images_dir[0]) + "\n")
			file.write(str(images_dir[1]) + "\n")



directory_images_clean = '../data/Sintel_final'
folders_sequences = listdir(directory_images_clean)
for i in range(len(folders_sequences)):
	directory_sequence = join(directory_images_clean, folders_sequences[i])
	images = sorted([f for f in listdir(directory_sequence) if isfile(join(directory_sequence, f)) and f.endswith(".png")])
	images_dir = [join(directory_sequence, f) for f in images]
	#print images_dir
	for j in range(2, len(images_dir) - 1):
		nums = str(j).zfill(2) + str(j + 1).zfill(2) + str(j - 1).zfill(2) + str(j + 2).zfill(2)
		
		#print directory_sequence + "/frames_" + nums + ".txt"
		with open(directory_sequence + "/frames_" + nums + ".txt", 'w') as file:
			file.write(str(images_dir[j - 1]) + "\n")
			file.write(str(images_dir[j]) + "\n")	
			file.write(str(images_dir[j - 2]) + "\n")			
			file.write(str(images_dir[j + 1]) + "\n")
	#Firts two images	
	nums = "0102"
	#print directory_sequence + "/frames_" + nums + ".txt"
	with open(directory_sequence + "/frames_" + nums + ".txt", 'w') as file:
			file.write(str(images_dir[0]) + "\n")
			file.write(str(images_dir[1]) + "\n")	
	nums = str(len(images_dir) - 1).zfill(2) + str(len(images_dir)).zfill(2)
	#Last two images
	#print directory_sequence + "/frames_" + nums + ".txt"
	with open(directory_sequence + "/frames_" + nums + ".txt", 'w') as file:
			file.write(str(images_dir[0]) + "\n")
			file.write(str(images_dir[1]) + "\n")
