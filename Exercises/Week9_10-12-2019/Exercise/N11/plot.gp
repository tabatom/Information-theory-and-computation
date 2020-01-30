# N	: input; is the number of particles

set terminal pdf size 16,9 linewidth 1
set output "4_eigval_vs_lambda.pdf"

set xrange [0:3.1]

set lmargin at screen 0.05

set ylabel "Energy [a.u.]" font ",25" offset -3.5,0
set xlabel "Lambda [a.u.]" font ",25"

set xtics font ",25"
set ytics font ",25"

set key font ",25"

plot "eig_val_results.txt" u 1:2 title "1st eigenvalue" w lp lw 2 pt 7, "eig_val_results.txt" u 1:3 title "2nd eigenvalue" w lp lw 2 pt 7 lc rgb "#2222ee", "eig_val_results.txt" u 1:4 title "3rd eigenvalue" w lp lw 2 pt 7 lc rgb "#00dddd", "eig_val_results.txt" u 1:5 title "4th eigenvalue" w lp lw 2 pt 7, [0:2] N*(-1-x**2/4) w l lw 2 lc rgb "#000000" title "Theoritical result", [2:4] -N*x w l lw 2 lc rgb "#000000" title ""
