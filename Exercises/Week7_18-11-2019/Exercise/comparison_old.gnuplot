time = time + 1

set for [i=0:2] ytics (sprintf("%1.1f", 0) 12./10.*i)
set for [i=0:2] ytics add (sprintf("%1.1f", 1) 12./10.*i+1.)
set for [i=0:2] ytics add (sprintf("%1.1f", 0.5) 12./10.*i+0.5)

set yrange[-0.1:3.1]
plot 'temp_evolut_mod.txt' u 1:time w l lc rgb "#ff9900" lt 1 lw 2 title sprintf("Single step time evolution: t= %i", time-1), 'temp_evolut_mod_theo.txt' u 1:time w l lc rgb "#0000ff" lt 0 lw 2 title sprintf("Multi step time evolution:  t= %i", time-1), 'temp_evolut_mod.txt' u 1:(column(time)+1.2) w l lc rgb "#ff9900" lt 1 lw 2 title sprintf("Shifted s.s.t.evo.: t= %i", time-1), 'temp_evolut_mod_theo.txt' u 1:(column(time)+2.4) w l lc rgb "#0000ff" lt 1 lw 2 title sprintf("Shifted m.s.t.evo. t= %i", time-1)

if (time<tot_time) reread
