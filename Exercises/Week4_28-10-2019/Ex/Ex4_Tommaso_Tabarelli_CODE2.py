#!/usr/bin/python

import os
import subprocess as sub
import sys
import numpy as np

file_list = ["results_m.txt", "results_T.txt", "results_F.txt"]
existing_files = []
labels = []

path = os.getcwd()
temp_name = path+"/results.txt"

for filename in file_list:
	if os.path.isfile(filename):
		existing_files.append(filename)
		# Creating a label for every different file
		labels.append(filename[filename.find(".txt")-2:filename.find(".txt")])

# Executing task
for file_,label in zip(existing_files,labels):
	
	file_ = path+"/results_m.txt"
	# Renaming file to use Gnuplot script
	sub.run(str("mv "+file_+" "+temp_name), shell=True)
	os.system(str('gnuplot "Gnu_fit.gnu"'))
	
	# Resetting original name
	sub.run(str("mv "+temp_name+" "+file_), shell=True)
	
	# Renaming files created by Gnuplot
	sub.run(str("mv "+path+"/fit_time_dim.log "+path+"/fit_time_dim"+label+".log"), shell=True)
	sub.run(str("mv "+path+"/quadratic_fit_coef.txt "+path+"/quadratic_fit_coef"+label+".txt"), shell=True)
	sub.run(str("mv "+path+"/cubic_fit_coef.txt "+path+"/cubic_fit_coef"+label+".txt"), shell=True)
	sub.run(str("mv "+path+"/Fit.png "+path+"/Fit"+label+".png"), shell=True)

os.system(str('gnuplot "plot_res.gnu"'))
