
	1. In the first exercise I just wrote a program that write in terminal "Hello world".

	2.	a) In the second exercise, after writing the first version program, the compiler can recognize that there is an "out of range" error and it does not compile. Newer versions suggest this:
			Error: Arithmetic overflow converting INTEGER(4) to INTEGER(2) at (1). This check can be disabled with the option ‘-fno-range-check’
		Using that flag the program runs but gives the wrong result.
		
		b) Using INT*4 the compiler gives no errors.

	3. The program is set to work with 2 matrices that have dimensions 500*500 (arbitrary choice). To compile there is a bash script that compiles the file using:
		- no optimization
		- the option -O1
		- the option -O2
		- the option -O3
	Then, there is another script that runs the executable and write results on a file called "my_res.txt"
