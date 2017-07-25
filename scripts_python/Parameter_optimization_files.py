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
from numpy import random

random_trials = 50
init_num = 0


lambda_range = [10, 100]
theta_range = [0.05, 1]
tau = [0.125]
beta_range = [0.125, 0.25, 1]
alpha_range = [0.001, 0.005]
tau_u_range = [0.001, 0.125]
tau_eta_range = [0.001, 0.125]
tau_chi_range = [0.001, 0.125]
mu_range = [0.05, 2]

directory = "../Parameters_optimization_files/"
name_file = directory + "params_"
if not os.path.exists(directory):
	os.makedirs(directory)

def write_list(list_parameters, filename, count):
	print(name_file + str(count) + ".txt")
	with open(name_file + str(count) + ".txt", 'w') as file:
		for i in range(len(list_parameters)):
			file.write(str(list_parameters[i]) + "\n")



for count in range(init_num, init_num + random_trials):
	
	parameters = [PAR_DEFAULT_LAMBDA, PAR_DEFAULT_THETA, PAR_DEFAULT_TAU, PAR_DEFAULT_BETA, PAR_DEFAULT_ALPHA, PAR_DEFAULT_TAU_U, PAR_DEFAULT_TAU_ETA, PAR_DEFAULT_TAU_CHI, PAR_DEFAULT_MU]
	parameters[0] = (lambda_range[1] - lambda_range[0]) * random.random_sample() + lambda_range[0]
	parameters[1] = (theta_range[1] - theta_range[0]) * random.random_sample() + theta_range[0]
	parameters[2] = PAR_DEFAULT_TAU
	parameters[3] = (beta_range[1] - beta_range[0]) * random.random_sample() + beta_range[0]
	parameters[4] = (alpha_range[1] - alpha_range[0]) * random.random_sample() + alpha_range[0]
	parameters[5] = (tau_u_range[1] - tau_u_range[0]) * random.random_sample() + tau_u_range[0]
	parameters[6] = (tau_eta_range[1] - tau_eta_range[0]) * random.random_sample() + tau_eta_range[0]
	parameters[7] = (tau_chi_range[1] - tau_chi_range[0]) * random.random_sample() + tau_chi_range[0]
	parameters[8] = (mu_range[1] - mu_range[0]) * random.random_sample() + mu_range[0]
	
	write_list(parameters, name_file, count)






