# filename, fit_log, img_name, data_title, fit_save ARE ALL INPUTS VARIABLE
#
# To pass input variables to gnuplot scripts use following command from
#	shell:
#	gnuplot -e "filename='example.txt'; fit_log='example2.log'; " fit_hist.gnu

#if (!exists("filename"))
#	{
#	print("\n\tERROR: file not found")
#	print(filename)
#	print("\n\tQuiting...\n")
#	exit
#	}

set key top right
# Shutting down real-time log for fit
set fit quiet
set fit logfile fit_log

P(x) = a*(x**alpha)*exp(-b*(x**beta))

a = 7
alpha = 2
b = 2
beta = 1

fit P(x) filename u 1:2 via a, alpha, b, beta

save fit fit_save

set terminal png size 960,640
set output img_name
plot P(x) w l title "Fit function" lw 2 lt rgb "#dd8822", filename u 1:2 w p pt 7 ps 1 lc rgb "blue" title data_title
