set key top left
set fit logfile "fit_time_dim_m.log"
cubicm(x) = a3m + b3m*x + c3m*x**2 + d3m*x**3
fit cubicm(x) "results_m.txt" u 1:2 via a3m,b3m,c3m,d3m
set print "fit_coef_m.txt"
print a3m
print b3m
print c3m
print d3m
set terminal png size 1920,1280
set output "Fit_m.png"
plot cubicm(x) w l title "Function fit" lw 2 lt rgb "#dd8822", "results_m.txt" w p pt 7 ps 3 lc rgb "blue" title "Manual data"


set key top left
set fit logfile "fit_time_dim_T.log"
cubicT(x) = a3T + b3T*x + c3T*x**2 + d3T*x**3
fit cubicT(x) "results_T.txt" u 1:2 via a3T,b3T,c3T,d3T
set print "fit_coef_T.txt"
print a3T
print b3T
print c3T
print d3T
set terminal png size 1920,1280
set output "Fit_T.png"
plot cubicT(x) w l title "Function fit", "results_T.txt" w p pt 7 ps 3 lc rgb "blue" title "Manual transposed data"

set key top left
set fit logfile "fit_time_dim_F.log"
cubicF(x) = a3F + b3F*x + c3F*x**2 + d3F*x**3
fit cubicF(x) "results_F.txt" u 1:2 via a3F,b3F,c3F,d3F
set print "fit_coef_F.txt"
print a3F
print b3F
print c3F
print d3F
set terminal png size 1920,1280
set output "Fit_F.png"
plot cubicF(x) w l title "Function fit", "results_F.txt" w p pt 7 ps 3 lc rgb "blue" title "Matmul data"

set terminal png size 1920,1280
set output "Fit_comp.png"
plot cubicm(x) w l lw 2 lt rgb "#ff00ff" title "Function fit", "results_m.txt" w p pt 7 ps 3 lc rgb "blue" title "Manual data", cubicT(x) w l lw 2 lt rgb "#aa0000" title "Function fit", "results_T.txt" w p pt 7 ps 3 lc rgb "black" title "Manual transposed data", cubicF(x) w l lw 2 lt rgb "#0000aa" title "Function fit", "results_F.txt" w p pt 7 ps 3 lc rgb "green" title "Matmul data"
