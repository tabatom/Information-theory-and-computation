set terminal pdf size 16,9 linewidth 1
set output "eigval_vs_lambda.pdf"

set xrange [0:3.1]
set yrange [-3.5:0]

set lmargin at screen 0.05

set ylabel "Energy/particle [a.u.]" font ",25" offset -3.5,0
set xlabel "Lambda [a.u.]" font ",25"

set xtics font ",25"
set ytics font ",25"

set key font ",25"

plot "Results.txt" u 1:2 title "N=2" w lp lw 2 pt 7, "Results.txt" u 1:3 title "N=3" w lp lw 2 pt 7 lc rgb "#2222ee", "Results.txt" u 1:4 title "N=4" w lp lw 2 pt 7 lc rgb "#00dddd", "Results.txt" u 1:5 title "N=5" w lp lw 2 pt 7, [0:2] (-1-x**2/4) w l lw 2 lc rgb "#000000" title "Mean field result", [2:4] -x w l lw 2 lc rgb "#000000" title ""
