set key top left
set terminal png size 1920,1280
set output "Results.png"
plot "results_m.txt" title "Manual" w lp, "results_T.txt" title "Transposed" w lp, "results_F.txt" title "Matmul" w lp

parabm(x) = a2m + b2m*x + c2m*x**2
cubicm(x) = a3m + b3m*x + c3m*x**2 + d3m*x**3
fit parabm(x) "results_m.txt" u 1:2 via a2m,b2m,c2m
fit cubicm(x) "results_m.txt" u 1:2 via a3m,b3m,c3m,d3m


parabT(x) = a2T + b2T*x + c2T*x**2
cubicT(x) = a3T + b3T*x + c3T*x**2 + d3T*x**3
fit parabT(x) "results_T.txt" u 1:2 via a2T,b2T,c2T
fit cubicT(x) "results_T.txt" u 1:2 via a3T,b3T,c3T,d3T


parabF(x) = a2F + b2F*x + c2F*x**2
cubicF(x) = a3F + b3F*x + c3F*x**2 + d3F*x**3
fit parabF(x) "results_F.txt" u 1:2 via a2F,b2F,c2F
fit cubicF(x) "results_F.txt" u 1:2 via a3F,b3F,c3F,d3F

set terminal png size 1920,1280
set output "Fit_comp.png"
plot cubicm(x) w l lw 2 lt rgb "#ff00ff" title "Function fit", "results_m.txt" w p pt 7 ps 3 lc rgb "blue" title "Manual data", cubicT(x) w l lw 2 lt rgb "#aa0000" title "Function fit", "results_T.txt" w p pt 7 ps 3 lc rgb "black" title "Manual transposed data", cubicF(x) w l lw 2 lt rgb "#0000aa" title "Function fit", "results_F.txt" w p pt 7 ps 3 lc rgb "green" title "Matmul data"
