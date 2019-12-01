time = time + 1

### Start multiplot (2x2 layout)
set multiplot layout 2,2 rowsfirst

# --- GRAPH a
set yrange[-0.1:1.1]
#set label 1 'a' at graph 0.92,0.9 font ',8'
plot 'temp_evolut_mod.txt' u 1:time w l ls 1 lc rgb "#aa0000" lw 2 title sprintf("Single step time evolution: t= %i", time-1), 'pot_evolut.txt' u 1:time w l lc rgb "#00dd00" lt 1 lw 2  title sprintf("Potential time evolution: t=%i", time-1)
# --- GRAPH b
set yrange[-0.1:1.1]
#set label 1 'b' at graph 0.92,0.9 font ',8'
plot 'temp_evolut_mod_theo.txt' u 1:time w l ls 1 lc rgb "#aa0000" lw 2 title sprintf("Multi step time evolution:  t= %i", time-1), 'pot_evolut.txt' u 1:time w l lc rgb "#00dd00" lt 1 lw 2  title sprintf("Potential time evolution: t=%i", time-1)
# --- GRAPH c
set yrange[-1.1:1.1]
#set label 1 'c' at graph 0.92,0.9 font ',8'
plot 'temp_evolut.txt' u 1:(column(2*time-2)) w l ls 1 lc rgb "#aa0000" lw 2 title sprintf("Real part s.s.t.evo.:  t= %i", time-1)
# --- GRAPH d
set yrange[-1.1:1.1]
#set label 1 'd' at graph 0.92,0.9 font ',8'
plot 'temp_evolut.txt' u 1:(column(2*time-1)) w l ls 1 lc rgb "#aa0000" lw 2 title sprintf("Immaginary part s.s.t.evo.:  t= %i", time-1)

unset multiplot
### End multiplot

if (time<tot_time) reread
