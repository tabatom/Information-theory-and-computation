cd ../Exercise/N2
gnuplot -e "N=2" plot.gp
cd ../N5
gnuplot -e "N=5" plot.gp
cd ../N10
gnuplot -e "N=10" plot.gp
cd ../N11
gnuplot -e "N=11" plot.gp
cd ../N12
gnuplot -e "N=12" plot.gp
cd ../N13
gnuplot -e "N=13" plot.gp

cd ../../Report

mv ../Exercise/N2/4_eigval_vs_lambda.pdf ./4_eigval_vs_lambda_N2.pdf
mv ../Exercise/N5/4_eigval_vs_lambda.pdf ./4_eigval_vs_lambda_N5.pdf
mv ../Exercise/N10/4_eigval_vs_lambda.pdf ./4_eigval_vs_lambda_N10.pdf
mv ../Exercise/N11/4_eigval_vs_lambda.pdf ./4_eigval_vs_lambda_N11.pdf
mv ../Exercise/N12/4_eigval_vs_lambda.pdf ./4_eigval_vs_lambda_N12.pdf
mv ../Exercise/N13/4_eigval_vs_lambda.pdf ./4_eigval_vs_lambda_N13.pdf
