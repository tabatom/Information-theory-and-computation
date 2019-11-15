# Arguments are:
#	L		: interval is [-L,L]
#	dim		: the dimension of the matrix (i.e. the number of points used
#				to evaluate the eigenfunction approximation)
#	n_eigv	: the n-th eigenvector (should be < 10)
#	img_name	: the name of the image
#	y_label_	: the title of y label

Lm = -L

set xrange[0:dim]

set key top left
set xlabel "Eigenfunction"
set ylabel y_label_
set terminal png size 960,640
set output img_name
set xtics Lm,10,L
plot "evect_theo.txt" u n_eigv w l lw 2 lc rgb "#9999aa" title "Theoretical eigenfunction", "evect_exp.txt" u n_eigv w p pt 7 ps 1.3 lc rgb "#0000aa" title "Approximated eigenfunctions"
