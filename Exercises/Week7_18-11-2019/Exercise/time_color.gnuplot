#set palette defined (0 '#00aa00', 0.5 '#555555', 1 '#ff0000')
set palette defined (0 '#0000ff',  1 '#dd7700')
set terminal png size 1920,1280
set output 'first_time_steps.png'
set xlabel 'x [a.u.]' font ",20"
set ylabel '|ψ(t)|²' font ",20"

set yrange[-0.1:1.1]
set xrange[-5:5]

set xtics (-4,-3,-2,-1,0,1,2,3,4)
set grid
set title 'Time evolution of eigenfunction n=0' font ",30" offset 0,-1
set key right top
set key font ",20"
set xtics font "10,25"
set ytics font "10,30"


set cbrange [0:100]
set cblabel "Time steps" font ",30"
set cbtics font ",20"
#plot for [i=1:100] 'temp_evolut_mod.txt' using 1:(column(i+1)) w l lw 3 lt palette frac 1.0*i/100.0 title ""
plot for [i=1:25] 'temp_evolut_mod.txt' using 1:(column(4*i+1)) w l lw 3 lt palette frac 1.0*4*i/100.0 title ""
