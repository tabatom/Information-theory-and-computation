	#!/usr/bin/python

import os
import subprocess as sub
import sys
import numpy as np

# Taking as input:
#  L_min	: the min value of the interval width
#  L_max	: the max value of the interval width
#			interval is: [-L, L]
#  L_steps	: the number of linear(?) steps to go from L_min to L_max
#  k_min	: the min value of the matrix dimension
#  k_max	: the max value of the matrix dimension
#  k_steps	: the number of linear(?) step to go from k_min to k_max
L_min   = int(sys.argv[1])
L_max   = int(sys.argv[2])
L_steps = int(sys.argv[3])
k_min   = int(sys.argv[4])
k_max   = int(sys.argv[5])
k_steps = int(sys.argv[6])

L_list = np.linspace(L_min, L_max, L_steps)
k_list = np.logspace(np.log10(k_min), np.log10(k_max), k_steps, dtype=int)


# Initializing files
# writing first element of L_list to proper file
ofile = open("L.txt", "w+")
ofile.write(str(L_list[0]))
ofile.close()

# writing first element of k_list to proper file
ofile = open("k_split.txt", "w+")
ofile.write(str(k_list[0]))
ofile.close()

# Getting local path
path = os.getcwd()

# Getting current date and hour
current_date = str(sub.check_output(["date", "+'%d-%m-%y_%H:%M'"]))[3:-4].replace("/", "-")

print(current_date)

# Checking reasonable interval
if (np.any(L_list<0)):
	print("\n\tParameters \"L\" souhld be greater than 0...\n\tExiting program.\n")
	sys.exit()

# Checking reasonable dimensions
if (np.any(k_list<0)):
	print("\n\tNumber of intervals souhld be greater than 0...\n\tExiting program.\n")
	sys.exit()

# Creating directories names
for ll in L_list:
	for kk in k_list:
		dir_name = str(current_date+"__Res_L"+str(ll)+"_"+"k"+str(kk))
		test_dir = sub.run(str("test -d "+path+"/"+dir_name), shell=True, stdout=sub.PIPE)
		# If directory does not exist the returncode is 1, which turns into a TRUE boolean in python
		#  so reversing the result
		test_dir = not(bool(int(test_dir.returncode)))
		# Creating directory if not existing
		if (not(test_dir)):
			sub.run(str("mkdir "+path+"/"+dir_name), shell=True)
		else:
			print("The directories alredy exist.")
			print("If you want to overwrite previous results type 'y', otherwise type 'n' and wait about 60 seconds.")
			answer = input()
			while (answer.lower()!='n' and answer.lower()!='y'):
				print("Please, type only 'y' to overwrite results or 'n' to skip.")
				answer = input()
			if (answer=='n'):
				print("Exiting...")
				sys.exit()
		
		# ---------- DOING TASK ----------
		
		print("\n\tEvaluating L="+str(ll)+" k="+str(kk)+"...")
		
		# Writing to files the parameters that fortran script will read
		# writing first element of L_list to proper file
		ofile = open("L.txt", "w+")
		ofile.write(str(ll))
		ofile.close()

		# writing first element of k_list to proper file
		ofile = open("k_split.txt", "w+")
		ofile.write(str(kk))
		ofile.close()

		# Calling fortran executable
		sub.run(["./Ex6_Tommaso_Tabarelli_CODE.exe"])
		
		#break
		
		# Plotting the eigenvalues
		os.system("gnuplot plot_eig.gnu")
		
		# Moving image
		sub.run(str("mv "+path+"/eig_plot.png "+path+"/"+dir_name+"/eig_plot.png"), shell=True)
		
		# Evaluating the first 10 eigenfunctions
		for index_ in range(1,11):
			os.system("gnuplot -e \"L="+str(ll)+"; dim="+str(kk+1)+"; n_eigv="+str(index_)+"; img_name=\'eigf_"+str(index_)+"\'\" eigvect_comp.gnu")
			sub.run(str("mv "+path+"/eigf_"+str(index_)+" "+path+"/"+dir_name+"/eigf"+str(index_-1)+"_L"+str(ll)+"_k"+str(kk)+".png"), shell=True)

		# Renaming and moving results to proper directory
		sub.run(str("mv "+path+"/eval_exp.txt "+path+"/"+dir_name+"/eval_exp.txt"), shell=True)
		sub.run(str("mv "+path+"/eval_theo.txt "+path+"/"+dir_name+"/eval_theo.txt"), shell=True)
		sub.run(str("mv "+path+"/evect_exp.txt "+path+"/"+dir_name+"/evect_exp.txt"), shell=True)
		sub.run(str("mv "+path+"/evect_theo.txt "+path+"/"+dir_name+"/evect_theo.txt"), shell=True)
		

sys.exit()
