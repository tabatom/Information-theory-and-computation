set key top left
set fit logfile "fit_time_dim_T.log"
cubic(x) = a3 + b3*x + c3*x**2 + d3*x**3
fit cubic(x) "results_T.txt" u 1:2 via a3,b3,c3,d3
set print "fit_coef_T.txt"
print a3
print b3
print c3
print d3
set terminal png size 1920,1280
set output "Fit_T.png"
plot cubic(x) w l title "Function fit", "results_m.txt" w p pt 7 ps 3 lc rgb "blue" title "Manual transposed data"
