set palette defined (0 '#0000ff',  1 '#dd7700')
set cbrange [0:100]
set cbtics font ",10"
set cblabel "Time steps"

set tics font ",20"

plot 'pot_evolut.txt' u 1:time w l lc rgb "#00dd00" lt 1 lw 3 title sprintf("Potential time evolution: t=%i", time-1), for [i=1:int(time/4)] 'temp_evolut_mod.txt' using 1:(column(4*i+1)) w l lw 3 lt palette frac 1.0*4*i/100.0 title "", 'temp_evolut_mod.txt' u 1:time w l lc rgb "#ee0000" lt 1 lw 3 title sprintf("Single step time evolution: t= %i", time-1)

#plot for [i=1:100] 'temp_evolut_mod.txt' using 1:(column(i+1)) w l lw 3 lt palette frac 1.0*i/100.0 title ""
