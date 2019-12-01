# Plotting time evolution of wave function in a gif format
# Input parameters:
#	tot_time_steps

reset

set term gif animate size 1920,1280

set xtics (-4,-3,-2,-1,0,1,2,3,4)

tot_time = tot_time_steps + 1 + 1

set grid

# ---------- Plotting the module time evolution ----------

set output "Wave_func_evo.gif"

set tics font ",40"

set ylabel "|ψ(t)|²" font ",40"
set xlabel 'x [a.u.]' font ",40"

time = 1

load "comparison.gnuplot"
#960,640
#1920,1280

# ---------- Plotting the time evolution of module, real part, immaginary part ----------

set term gif animate size 960,640

set tics font ",20"

set ylabel "|ψ(t)|²" font ",20"
set xlabel 'x [a.u.]' font ",20"

set output "Wave_func_comp.gif"

time = 1

load "4_plots.gnuplot"

# ---------- Plotting the track of time evolution (using colors) ----------

set term gif animate size 960,640

set tics font ",20"

set ylabel "|ψ(t)|²" font ",20"
set xlabel 'x [a.u.]' font ",20"

set output "Wave_func_track.gif"

time = 1

load "wave_track.gnuplot"

# ---------- Plotting the time evolution using colors ----------

# Loading script to make single image
load "time_color.gnuplot"
