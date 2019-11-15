set key top left
set fit logfile "fit_time_dim.log"
parab(x) = a2 + b2*x + c2*x**2
cubic(x) = a3 + b3*x + c3*x**2 + d3*x**3
fit parab(x) "results.txt" u 1:2 via a2,b2,c2
fit cubic(x) "results.txt" u 1:2 via a3,b3,c3,d3
set print "quadratic_fit_coef.txt"
print a2
print b2
print c2
set print "cubic_fit_coef.txt"
print a3
print b3
print c3
print d3
set terminal png size 1920,1280
set output "Fit.png"
plot parab(x) w l title "Function quadratic fit" lw 2 lt rgb "#dd8822", cubic(x) w l title "Function cubic fit" lw 2 lt rgb "#000000", "results.txt" w p pt 7 ps 3 lc rgb "blue" title "Data"
