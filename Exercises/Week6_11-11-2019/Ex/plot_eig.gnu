set key top left
set xlabel "n"
set ylabel "Eigenvalues"
set terminal png size 960,640
set output "eig_plot.png"
plot "ev_exp.txt" w p pt 7 ps 1.3 lc rgb "#9999aa" title "Approximated eigenvalues", "ev_theor.txt" w p pt 7 ps 1.3 lc rgb "#0000aa" title "Theoretical eigenvalues"
