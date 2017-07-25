#! /usr/bin/python


from os import listdir
from os.path import isfile, join
import numpy as np
images_dir = []
percentage_random = 10.0


directory_images_clean = '../data/Sintel_final'
folders_sequences = listdir(directory_images_clean)
print len(images_dir)
for i in range(len(folders_sequences)):
	directory_sequence = join(directory_images_clean, folders_sequences[i])
	images = sorted([f for f in listdir(directory_sequence) if isfile(join(directory_sequence, f)) and f.endswith(".txt")])
	images_dir = images_dir + [join(directory_sequence, f) for f in images]
	print len(images_dir)
#print images_dir
directory_images_clean = '../data/Sintel_clean'
folders_sequences = listdir(directory_images_clean)
for i in range(len(folders_sequences)):
	directory_sequence = join(directory_images_clean, folders_sequences[i])
	images = sorted([f for f in listdir(directory_sequence) if isfile(join(directory_sequence, f)) and f.endswith(".txt")])
	images_dir = images_dir + [join(directory_sequence, f) for f in images]
	print len(images_dir)


len_random = np.int(np.ceil(len(images_dir)*percentage_random/100.0))
print len_random
random_sample_images = np.random.choice(images_dir, size = len_random, replace = False)


with open("../scripts_python/Random_images.txt", 'w') as file:
	for image in random_sample_images:
		file.write(image + "\n")

