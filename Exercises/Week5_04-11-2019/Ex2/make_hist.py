import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import subprocess as sub
import sys

# ---------- LOADING FILES USING ARGUMENTS AND APPENDING THEM TO A UNIQUE ARRAY ----------

# Taking as input:
#  dir_	: the directory from which take the data files
#			it should be only 1 at a time
#
# Directory name should be:
#	"TYPE_dim****_num****__CURRENTDATE"
#	TYPE		: "real" or "comp"
#	CURRENTDATE	: automatically saved by script when the directory was created
#	****		: numbers, any length is allowed
#	dim, num	: fixed pieces of directory name

dir_ = str(sys.argv[1])
dim_ = int(sys.argv[2])


# Getting type of matrix (comp or real) from directory name
data_type = dir_[:4]

# Getting number of files from directory name
num_files = int(dir_[dir_.find("num")+3:dir_.find("__")])

# Current dir
path = os.getcwd()


# Reading spacings to eventually store other results
intervals = []

ofile = open("spacings.txt", "r")
for line in ofile:
	intervals.append(int(line))
ofile.close()

intervals = np.array(intervals)

data_s = []
# Creating proper arrays to store results
for ii in range(len(intervals)+1):	# The +1 counts also the "standard" non-local spacings
	data_s.append(np.array([]))


# Reading files and concatenating them

# Standard spacings
for jj in range(num_files):
	current_data = np.loadtxt(str(dir_+"/s_i_distr"+str(jj)+".txt"))
	data_s[0] = np.append(data_s[0], current_data)

# Other spacings
for ii in range(len(intervals)):
	for jj in range(num_files):
		current_data = np.loadtxt(str(dir_+"/s_i_distr_interv"+str("{:04d}".format(intervals[ii]))+"_"+str(jj)+".txt"))
		data_s[ii+1] = np.append(data_s[ii+1], current_data)

# This check was correct last time ( should give (dim-1)*num_files )
#print(data_s.shape)


# ---------- BULDING HISTOGRAM ----------

my_hist = []
bin_edges = []


for ii in range(len(intervals)+1):

	my_hist, bin_edges = np.histogram(data_s[ii], density=True, bins='auto')
	# Evaluating bin centers (summing bin edges and dividing by 2)
	bin_centres = (bin_edges[1:]+bin_edges[:-1])/2

	# Building a dataset with:
	#	"x"	: bin centers
	#	"y"	: bin height
	#  so that one can use it as an approximation of a distribution

	graph_ = pd.DataFrame([bin_centres, my_hist]).T
	
	# Saving histogram results in "hist.txt" file in the same folder
	if (ii==0):
	
		graph_.to_csv(str(path+"/"+dir_+"/hist_standard.txt"), header=None, index=None, sep="\t")
		print("\n\tSaved \"hist._standard.txt\" at \""+dir_+"\" folder.")
	else:
		graph_.to_csv(str(path+"/"+dir_+"/hist_interv"+str("{:04d}".format(intervals[ii-1]))+".txt"), header=None, index=None, sep="\t")
		print("\n\tSaved \"hist_interv"+str("{:04d}".format(intervals[ii-1]))+".txt\" at \""+dir_+"\" folder.")

print()
