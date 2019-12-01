# Arguments are:
#	L		: interval is [-L,L]
#	dim		: the dimension of the matrix (i.e. the number of points used
#				to evaluate the eigenfunction approximation)
#	n_eigv	: the n-th eigenvector (should be < 10)
#	img_name	: the name of the image

set xrange[0:dim]
set yrange[-1:1]

array tics_[11]

do for [i=0:10] {tics_[i+1]=-L+i*(2*L/10)}
set key top left
set xlabel "Interval position"
set ylabel "y"
set terminal png size 960,640
set output img_name
n_eigv = int(n_eigv)
set for [i=0:10] xtics (sprintf("%1.2f", tics_[i+1]) i*dim/10)
plot "evect_exp.txt" u n_eigv w p pt 7 ps 0.7 lc rgb "#0000aa" title "Approximated eigenfunctions", "evect_exp.txt" using (column(int(n_eigv)))*(-1) w p pt 7 ps 0.7 lc rgb "#aa0000" title "Approximated eigenfunctions", "evect_theo.txt" u n_eigv w l lw 2 lc rgb "#22aa33" title "Theoretical eigenfunction"
