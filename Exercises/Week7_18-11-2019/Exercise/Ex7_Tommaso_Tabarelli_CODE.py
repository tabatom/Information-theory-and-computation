#!/usr/bin/python

import os
import subprocess as sub
import sys
import numpy as np

# Taking as input:
#  L_lim		: the SPACE grid is [-L_lim, L_lim]
#  k_split		: the number of intervals of the SPACE grid
#  T_tot		: the TIME grid is [0, T_tot]
#  time_split	: the number of intervals of the TIME grid
#  gnuplot_steps	: the number of frames to generate the ".gif" file
L_lim		= float(sys.argv[1])
k_split	= int(sys.argv[2])
T_tot		= float(sys.argv[3])
time_split	= int(sys.argv[4])

gnuplot_steps = int(sys.argv[5])

delta_x	= float(2*L_lim/k_split)
delta_t	= float(T_tot/time_split)


# ---------- DOING "REASONABLE" CHECKS ----------

# Checking reasonable SPACE interval
if (L_lim<0):
	print("\n\tParameters \"L\" souhld be greater than 0...\n\tExiting program.\n")
	sys.exit()

# Checking reasonable SPACE splitting
if (k_split<0):
	print("\n\tNumber of SPACE intervals souhld be greater than 0...\n\tExiting program.\n")
	sys.exit()

# Checking reasonable TIME interval
#	No main reason to consider only positive times, but conventionally taking them > 0
if (T_tot<0):
	print("\n\tParameters \"T_tot\" souhld be greater than 0...\n\tExiting program.\n")
	sys.exit()

# Checking reasonable TIME splitting
if (time_split<0):
	print("\n\tNumber of TIME intervals souhld be greater than 0...\n\tExiting program.\n")
	sys.exit()


# ---------- Printing useful informations ----------
print("\n------ PARAMETERS ----------\n")
print("Space interval: [-"+str(L_lim)+","+str(L_lim)+"]")
print("Number of SPACE intervals:",k_split)
print("Δx:",delta_x)
print("Time interval: [0,"+str(T_tot)+"]")
print("Number of TIME intervals:",time_split)
print("Δt:",delta_t)
print()


# ---------- Initializing files ----------

# writing L_lim to proper file
ofile = open("L.txt", "w+")
ofile.write(str(L_lim))
ofile.close()

# writing k_split to proper file
ofile = open("k_split.txt", "w+")
ofile.write(str(k_split))
ofile.close()

# writing T_tot to proper file
ofile = open("T_tot.txt", "w+")
ofile.write(str(T_tot))
ofile.close()

# writing time_split to proper file
ofile = open("time_split.txt", "w+")
ofile.write(str(time_split))
ofile.close()


# ---------- PREPARING DIRECTORY ----------

# Getting local path
path = os.getcwd()

# Getting current date and hour
current_date = str(sub.check_output(["date", "+'%y-%m-%d_%H:%M'"]))[3:-4].replace("/", "-")

print("Current date:", current_date, end="\n\n")

# Creating directory name
dir_name = str(current_date+"__L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t))

test_dir = sub.run(str("test -d "+path+"/"+dir_name), shell=True, stdout=sub.PIPE)

# If directory does not exist the returncode is 1, which turns into a TRUE boolean in python
#  so here the code read the reversed result
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

print("\nEvaluating:"+"\n"+"L\t= "+str(L_lim)+"\nΔx\t= "+str(delta_x)+"\nT_tot\t= "+str(T_tot)+"\nΔt\t= "+str(delta_t))

# Calling fortran executable
sub.run(["./Ex7_Tommaso_Tabarelli_CODE.exe"])

print("\nDone...\n")

print("Plotting results in gif files...\n")

os.system("gnuplot -e \"tot_time_steps="+str(gnuplot_steps)+"\" gif_maker.gp")


# ---------- Moving results to the proper folder ----------

# Moving data files
sub.run(str("mv "+path+"/pot_evolut.txt "+path+"/"+dir_name+"/pot_evolut__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".txt"), shell=True)
sub.run(str("mv "+path+"/temp_evolut.txt "+path+"/"+dir_name+"/temp_evolut__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".txt"), shell=True)
sub.run(str("mv "+path+"/temp_evolut_mod.txt "+path+"/"+dir_name+"/temp_evolut_mod__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".txt"), shell=True)
sub.run(str("mv "+path+"/temp_evolut_mod_theo.txt "+path+"/"+dir_name+"/temp_evolut_mod_theo__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".txt"), shell=True)

# Moving gif files
sub.run(str("mv "+path+"/Wave_func_evo.gif "+path+"/"+dir_name+"/Wave_func_evo__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".gif"), shell=True)
sub.run(str("mv "+path+"/Wave_func_comp.gif "+path+"/"+dir_name+"/Wave_func_comp__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".gif"), shell=True)
sub.run(str("mv "+path+"/Wave_func_track.gif "+path+"/"+dir_name+"/Wave_func_track__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".gif"), shell=True)

# Moving image file
sub.run(str("mv "+path+"/first_time_steps.png "+path+"/"+dir_name+"/col_t_steps__"+"L"+str(L_lim)+"_deltax"+str(delta_x)+"_Ttot"+str(T_tot)+"_deltat"+str(delta_t)+".png"), shell=True)

print("\nDone\n")

sys.exit()
