time = time + 1

set yrange[-0.1:1.1]

load "track_step.gnuplot" 

if (time<tot_time) reread
