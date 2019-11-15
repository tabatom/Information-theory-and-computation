#!/usr/bin/python

import os
import subprocess as sub
import sys
import numpy as np

N_min = int(sys.argv[1])
N_max = int(sys.argv[2])
Step = int(sys.argv[3])

# Checking reasonable dimensions
if (N_min<=0 or N_max<=0):
	print("\n\tN_min/N_max souhld be greater than 0...\n\tStopping program.\n")
	sys.exit()

# Checking N_min <= N_max
if (N_min>N_max):
	print("\n\tN_min sould be smaller than N_max...\n\tStopping program.\n")
	sys.exit()

# Creating a logarithmic spacing (it spans more magnitudes than linear spacing)
# dim_list = np.logspace( np.log10(N_min), np.log10(N_max), 10, dtype=np.int16)
dim_list = np.arange( N_min, N_max, Step, dtype=np.int16)

# Cleaning old files
file_list = ["results_m.txt", "results_T.txt", "results_F.txt"]

for filename in file_list:
	if os.path.isfile(filename):
		os.remove(filename)

# Executing task
for ii in dim_list:
	ofile = open("dim_file.txt", "w+")
	ofile.write(str(ii))
	ofile.close()
	sub.run(["./Ex4_Tommaso_Tabarelli_CODE.exe"])
