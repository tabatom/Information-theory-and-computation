for T in 0.1 1 10 100;
	do
	for tsteps in 10 100 500;
		do python Ex7_Tommaso_Tabarelli_CODE.py 5 100 $T $tsteps $tsteps;
		done;
done
