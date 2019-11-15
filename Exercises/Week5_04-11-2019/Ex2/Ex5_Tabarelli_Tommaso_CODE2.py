#!/usr/bin/python

import os
import subprocess as sub
import sys
import numpy as np

# Taking as input:
#  dim_	: the dimension of the matrix
#  num_mat	: the number of matrices to be evaluated
dim_ = int(sys.argv[1])
num_mat = int(sys.argv[2])

# Writing dim_ to proper file
ofile = open("dim_file.txt", "w+")
ofile.write(str(dim_))
ofile.close()

path = os.getcwd()

current_date = str(sub.check_output(["date", "+'%d-%m-%y_%H:%M'"]))[3:-4].replace("/", "-")

print(current_date)

# Creating directories
dir_real = str("real_dim"+str(dim_)+"_"+"num"+str(num_mat)+"__"+current_date)
dir_comp = str("comp_dim"+str(dim_)+"_"+"num"+str(num_mat)+"__"+current_date)

sub.run(str("mkdir "+path+"/"+dir_comp), shell=True)
sub.run(str("mkdir "+path+"/"+dir_real), shell=True)

# Checking reasonable dimensions
if (dim_<=0):
	print("\n\tDimensions souhld be greater than 0...\n\tStopping program.\n")
	sys.exit()

intervals = [dim_/200, dim_/100, dim_/50, dim_/25, dim_/10, dim_/5, dim_/2]

intervals = np.array(intervals, dtype=int)

# Writing intervals to file
ofile = open("spacing_intervals.txt", "w")
for v in intervals:
	ofile.write(str("{:04d}".format(v)+"\n"))
ofile.close()


print("\nIntervals are: ", intervals)

# Executing task: evaluating matrices and storing results
for ii in range(num_mat):

	# Info on process
	print("\nDoing operation: "+str(ii+1)+"/"+str(num_mat), end="\n\n")
	
	# Operating on hermitian matrixa
	sub.run(["./Ex5_Tommaso_Tabarelli_CODE1.exe"])
	sub.run(str("mv "+path+"/s_i_distr_comp.txt "+path+"/"+dir_comp+"/s_i_distr"+str(ii)+".txt"), shell=True)
	
	# Saving results for other spacings
	for jj in range(len(intervals)):
		interv_str = str("{:04d}".format(intervals[jj]))
		# str("{:04d}".format(num)) is to force numbers to have 3 digits (Fortran is hard to
		#  deal with for what concerns TRIMming stings and stuff like that...
		sub.run(str("mv "+path+"/s_i_distr_comp_local_"+interv_str+".txt "+path+"/"+dir_comp+"/s_i_distr_interv"+interv_str+"_"+str(ii)+".txt"), shell=True)
	
	# Operating on "real matrix" (only vector diagonal)
	sub.run(["./Ex5_Tommaso_Tabarelli_CODE2.exe"])
	sub.run(str("mv "+path+"/s_i_distr_real.txt "+path+"/"+dir_real+"/s_i_distr"+str(ii)+".txt"), shell=True)
	# Saving results for other spacings
	for jj in range(len(intervals)):
		interv_str = str("{:04d}".format(intervals[jj]))
		# str("{:03d}".format(num)) is to force numbers to have 3 digits (Fortran is hard to
		#  deal with for what concerns TRIMming stings and stuff like that...
		sub.run(str("mv "+path+"/s_i_distr_real_local_"+interv_str+".txt "+path+"/"+dir_real+"/s_i_distr_interv"+interv_str+"_"+str(ii)+".txt"), shell=True)


# ---------- Calling other scripts: ----------

# Making histograms for both "comp" and "real"
#	(parameters are passed to .sh script and in their turn
#		passed to python scripts called by bash)
print(path)
#sub.run(str(path+"/python make_hist.py "+"./"+dir_real), shell=True)
#sub.run(str("python "+path+"/python make_hist.py "+"./"+dir_comp), shell=True)
#sub.run(str("source "+path+"/make_hist.sh "+"./"+dir_real+" ./"+dir_comp), shell=True)
#sub.call([path+"/make_hist.sh", "./"+dir_comp, "./"+dir_real])
os.system("sh make_hist.sh "+dir_comp+" "+dir_real+" "+str(dim_))

# Calling gnuplot scripts.
# INPUT VARIABLES:
#	filename
#	fit_log
#	img_name
#	data_title (to be used in legend in image)
# gnuplot -e "filename='...'; fit_log='...'; img_name='...'; data_title='...'" fit_hist.gnu
info_real = dir_real
info_comp = dir_comp

# FITTING REAL DIAGONAL MATRIX RESULTS
print("\n\tFITTING REAL DIAGONAL MATRIX RESULTS\n")

# Standard spacings
print("Standard approach")
os.system(str('gnuplot -e "filename=\''+dir_real+'/hist_standard.txt\'; fit_log=\''+dir_real+'/fit_real_standard.log\'; img_name=\''+dir_real+'/hist_standard_fit\'; fit_save=\''+dir_real+'/fit_real_standard_res.txt\'; data_title=\'real'+str(dim_)+' standard\'" "gnu_fit.gnu"'))

# Other spacings
for ii in range(len(intervals)):
	print("Intervals: " + str(intervals[ii]))
	interv_str = str("{:04d}".format(intervals[ii]))
	os.system(str('gnuplot -e "filename=\''+dir_real+'/hist_interv'+interv_str+'.txt\'; fit_log=\''+dir_real+'/fit_real_interv'+interv_str+'.log\'; img_name=\''+dir_real+'/hist_interval_'+interv_str+'_fit\'; fit_save=\''+dir_real+'/fit_real_interv'+interv_str+'_res.txt\'; data_title=\'real'+str(dim_)+' interval'+interv_str+'\'" "gnu_fit.gnu"'))


# FITTING HERMITIAN MATRIX RESULTS
print("\n\tFITTING HERMITIAN MATRIX RESULTS\n")

# Standard spacings
print("Standard approach")
os.system(str('gnuplot -e "filename=\''+dir_comp+'/hist_standard.txt\'; fit_log=\''+dir_comp+'/fit_comp_standard.log\'; img_name=\''+dir_comp+'/hist_standard_fit\'; fit_save=\''+dir_comp+'/fit_comp_standard_res.txt\'; data_title=\'comp'+str(dim_)+'_standard\'" "gnu_fit.gnu"'))

# Other spacings
for ii in range(len(intervals)):
	print("Intervals: " + str(intervals[ii]))
	interv_str = str("{:04d}".format(intervals[ii]))
	os.system(str('gnuplot -e "filename=\''+dir_comp+'/hist_interv'+interv_str+'.txt\'; fit_log=\''+dir_comp+'/fit_comp_interv'+interv_str+'.log\'; img_name=\''+dir_comp+'/hist_interval_'+interv_str+'_fit\'; fit_save=\''+dir_comp+'/fit_comp_interv'+interv_str+'_res.txt\'; data_title=\'comp'+str(dim_)+' interval '+interv_str+'\'" "gnu_fit.gnu"'))

