

import os

from os import listdir
from os.path import isfile, join


directory_images_clean = '../data/Sintel_final'
folders_sequences = listdir(directory_images_clean)
for i in range(len(folders_sequences)):
	directory_sequence = join(directory_images_clean, folders_sequences[i])
	images = sorted([f for f in listdir(directory_sequence) if isfile(join(directory_sequence, f))])
	images_dir = [join(directory_sequence, f) for f in images]

for file in listdir("/mydir"):
    if file.endswith(".txt"):
        print(join("/mydir", file))
