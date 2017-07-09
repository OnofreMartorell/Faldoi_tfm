#! /usr/bin/python




PAR_DEFAULT_LAMBDA = 40 #40
PAR_DEFAULT_THETA = 0.3
PAR_DEFAULT_TAU = 0.125 #0.25
PAR_DEFAULT_BETA =   0.5 #1
PAR_DEFAULT_ALPHA =  0.01
PAR_DEFAULT_TAU_U =  0.125
PAR_DEFAULT_TAU_ETA =  0.125
PAR_DEFAULT_TAU_CHI = 0.1
PAR_DEFAULT_MU = 1

import os

lambdas = [35, 50, 55]
thetas = [0.1, 0.2, 0.4, 0.5]
tau = [0.125]
betas = [0.125, 0.25, 1]
alphas = [0.001, 0.005]
taus_u = [0.1]
taus_eta = [0.1, 0.075]
taus_chi = [0.125, 0.075]
mus = [0.5, 1.5]
directory = "../Parameters_files/"
name_file = directory + "params_"
if not os.path.exists(directory):
	os.makedirs(directory)

def write_list(list_parameters, filename, count):
	print(filename + str(count) + ".txt")
	with open(filename + str(count) + ".txt", 'w') as file:
		for i in range(len(list_parameters)):
			file.write(str(list_parameters[i]) + "\n")


def list_params(directory, list_param, defaults, idx_param, count, filename):
	parameters = list(defaults)
	#print(parameters)
	for i in range(len(list_param)):
		parameters[idx_param] = list_param[i]
		print(parameters)
		count = count + 1
		write_list(parameters, filename, count)
	#print(count)	
	return count

parameters = [PAR_DEFAULT_LAMBDA, PAR_DEFAULT_THETA, PAR_DEFAULT_TAU, PAR_DEFAULT_BETA, PAR_DEFAULT_ALPHA, PAR_DEFAULT_TAU_U, PAR_DEFAULT_TAU_ETA, PAR_DEFAULT_TAU_CHI, PAR_DEFAULT_MU]
count = 0
count = list_params(directory, lambdas, parameters, 0, count, name_file)
count = list_params(directory, thetas, parameters, 1, count, name_file)
count = list_params(directory, tau, parameters, 2, count, name_file)
count = list_params(directory, betas, parameters, 3, count, name_file)
count = list_params(directory, alphas, parameters, 4, count, name_file)
count = list_params(directory, taus_u, parameters, 5, count, name_file)
count = list_params(directory, taus_eta, parameters, 6, count, name_file)
count = list_params(directory, taus_chi, parameters, 7, count, name_file)
count = list_params(directory, mus, parameters, 8, count, name_file)

