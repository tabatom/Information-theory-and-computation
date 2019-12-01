time = time + 1

set yrange[-0.1:1.1]
plot 'temp_evolut_mod.txt' u 1:time w l lc rgb "#ee0000" lt 1 lw 3 title sprintf("Single step time evolution: t=%i", time-1), 'pot_evolut.txt' u 1:time w l lc rgb "#00dd00" lt 1 lw 3 title sprintf("Potential time evolution: t=%i", time-1)

if (time<tot_time) reread
