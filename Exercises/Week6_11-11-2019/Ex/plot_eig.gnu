set key top left
set xlabel "n"
set ylabel "Eigenvalues"
set terminal png size 960,640
set output "eig_plot.png"
plot "eval_exp.txt" w p pt 7 ps 1 lc rgb "#9999aa" title "Approximated eigenvalues", "eval_theo.txt" w p pt 7 ps 1 lc rgb "#0000aa" title "Theoretical eigenvalues"
