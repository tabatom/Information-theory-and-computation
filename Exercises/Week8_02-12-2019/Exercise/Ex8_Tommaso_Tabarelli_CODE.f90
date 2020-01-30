	PROGRAM Ex8
	! Program that studies a system made by N generic subsystems: both interacting and non-interacting cases are analyzed
	! Variables are:
	!	real_temp		: variable to initialize real parts of matrix elements
	!	comp_temp		: variable to initialize imaginary parts of matrix elements
	!	num_par		: it is the number of particles. Parameters read from file "num_par.txt"
	!	dim_			: it is the dimension of every subsystem. Parameters read from file "dim_.txt"
	!	ii,jj			: simply some numbers to do loops
	!	ii_min, ii_max	: indexes to help storing psi_sep_temp elements
	!	res_t			: temporary result when dividing and "decoding" the position to create psi_sep_temp
	!	rem_t			: temporary remainder using when "decoding" the position to create psi_sep_temp
	!	trace			: variable used to check the trace of the density matrix (it should be 1)
	!	psi_sep_temp	: the "compressed" result when using separable states: it is a matrix dim_*num_par
	!	psi_sep		: the explicit wave function in the separable case
	!	psi_not_sep		: the explicit wave function in the non separable case
	!	rho_not_sep		: density matrix in the non separable case
	!	rho_sep		: density matrix in the separable case
	!	rho_squ		: squared rho_sep to check trace=1
	!	rho_red		: the reduced density matrix (this program only can evaluate it when num_par=2 and for generic dim_)
	!	red_dim		: the dimension of the reduced matrix
	!	particle		: a number used to select the particle to be traced out: should be 1 ≤ particle ≤ num_par
	!					N.B.: as aforementioned, the evaluation of the reduced matrix only works with num_par = 2
	!	ii_red		: index to evaluate the trace of the submatrices when evaluating the reduced matrix
	!	index_i,index_j	: variables to compres information about the indexes when evaluating the reduced matrix
	!	Debug_flag		: logical value to enable runtime error warnings and checkpointing
	!	Herm_flag		: logical value to test if the densiti matrices are hermitian
	!	counter		: integer variable to be used to count checkpoints steps

	USE DEBUG

	IMPLICIT NONE

	REAL*8 :: real_temp, comp_temp, trace
	INTEGER*8 :: num_par, dim_, ii, jj, ii_min, ii_max, res_t, rem_t
	COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi_sep_temp, psi_sep, psi_not_sep
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: rho_not_sep, rho_sep, rho_squ
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: rho_red
	INTEGER*8 :: red_dim, particle, ii_red, index_i, index_j
	LOGICAL :: Debug_flag, Herm_flag
	
	INTEGER*8 :: counter
	
	COMMON counter
	counter = 0

	! Debugging
	Debug_flag = .FALSE.


	! Opening file to read number of space intervals
	OPEN(11, file="num_par.txt", status='old', action='READ')
	! Reading number of space intervals from file
	READ(11,*) num_par
	! Closing file
	CLOSE(11)
	
	! Opening file to read total time
	OPEN(11, file="dim_.txt", status='old', action='READ')
	! Reading total time from file
	READ(11,*) dim_
	! Closing file
	CLOSE(11)
	
	
	! Checking dimensions greater than 0
	IF (Debug_flag.AND.(dim_.LE.0)) THEN
		CALL Print_error(Debug_flag, "Dim<=0")
	END IF
	
	! Checking number of particles
	IF (Debug_flag.AND.(dim_.LE.0)) THEN
		CALL Print_error(Debug_flag, "Num_partic<=0")
	END IF


	
! ---------- ALLOCATING VARIABLES ----------

	ALLOCATE(psi_sep_temp(dim_*num_par))
	
	ALLOCATE(psi_sep(dim_**num_par))
	ALLOCATE(psi_not_sep(dim_**num_par))
	
	ALLOCATE(rho_not_sep(dim_**num_par,dim_**num_par))
	ALLOCATE(rho_sep(dim_**num_par,dim_**num_par))
	
	ALLOCATE(rho_squ(dim_**num_par,dim_**num_par))



! ---------- DOING TASK ----------

! ---------- Wave functions ----------

	! Initializing NOT SEPARABLE wave function vector
	CALL Checkpoint(Debug_flag, 'Read', "psi_not_sep", psi_not_sep, dim_**num_par)
	DO ii=1,dim_**num_par
		CALL RANDOM_NUMBER(real_temp)
		CALL RANDOM_NUMBER(comp_temp)
		psi_not_sep(ii) = CMPLX(2*real_temp-1.,2*comp_temp-1)
	END DO
	
	! Normalizing
	psi_not_sep(:) = psi_not_sep(:) / SQRT(SUM(ZABS(psi_not_sep)**2))
	
	CALL Checkpoint(Debug_flag, 'Write', "psi_not_sep", psi_not_sep, dim_**num_par)

	
	! Checking wave function normalization
	WRITE(*,*) "Norm of psi (general case): ", SQRT(SUM(ZABS(psi_not_sep)**2))



	! Initializing SEPARABLE wave function vector (initializing only "single subsystem wave functions")
	CALL Checkpoint(Debug_flag, 'Read', "psi_sep_temp", psi_sep_temp, dim_*num_par)
	
	WRITE(*,*) "SEPARABLE CASE"
	WRITE(*,*) "Checking every wave function normalization:"
	
	DO ii=0,num_par-1
		DO jj=1,dim_
			CALL RANDOM_NUMBER(real_temp)
			CALL RANDOM_NUMBER(comp_temp)
			psi_sep_temp(ii*dim_+jj) = CMPLX(2*real_temp-1.,2*comp_temp-1)
		END DO
		ii_min = ii*dim_+1
		ii_max = ii*num_par+jj
		! Normalizing every single particle wave function
		psi_sep_temp(ii_min:ii_max) = psi_sep_temp(ii_min:ii_max) / SQRT(SUM(ZABS(psi_sep_temp(ii_min:ii_max)**2)))
		
		WRITE(*,*) SQRT(SUM(ZABS(psi_sep_temp(ii_min:ii_max)**2)))
	END DO
	CALL Checkpoint(Debug_flag, 'Write', "psi_sep_temp", psi_sep_temp, dim_*num_par)


! ---------- Density matrices ----------

! NOT SEPARABLE case

	! NOT SEPARABLE density matrix
	DO ii=1,dim_**num_par
		DO jj=1,dim_**num_par
			rho_not_sep(ii,jj) = psi_not_sep(ii)*CONJG(psi_not_sep(jj))
		END DO
	END DO


	! Checking trace of rho^2
	trace = 0
	rho_squ = MATMUL(rho_not_sep, rho_not_sep)
	
	DO ii=1,dim_**num_par
		trace = trace + rho_squ(ii,ii)
	END DO
	
	WRITE(*,*) "Checking trace of rho (general case)", trace


	! Checking the density matrix is hermitian
	Herm_flag = .TRUE.	
	DO ii=1,dim_**num_par
		DO jj=ii,dim_**num_par
			IF (rho_not_sep(ii,jj).NE.CONJG(rho_not_sep(jj,ii))) THEN
				Herm_flag = .FALSE.
			END IF
		END DO
	END DO

	WRITE(*,*) "Checking matrix is hermitian", Herm_flag



! SEPARABLE case

	! Bulding psi (the complete one)
	DO ii=1,dim_**num_par
	
		! The "code" should give the result "0" (all indexes to 0) at beginning, not at the end
		!	so here one has to start from 0, not from 1
		res_t = ii-1
		psi_sep(ii) = (1.,0.)
		
		
		DO jj=1,num_par
			rem_t = MOD(res_t,dim_)
			
			! Evaluating the product using this trick
			psi_sep(ii) = psi_sep(ii) * psi_sep_temp(dim_*(num_par - jj) + rem_t + 1)
			! Updating the number to divide
			res_t = (res_t-rem_t)/dim_
		END DO

	END DO


	! SEPARABLE density matrix
	DO ii=1,dim_**num_par
		DO jj=1,dim_**num_par
			rho_sep(ii,jj) = psi_sep(jj)*CONJG(psi_sep(ii))
		END DO
	END DO


	! Checking the trace of rho^2
	trace = 0
	rho_squ = MATMUL(rho_sep, rho_sep)
	
	DO ii=1,dim_**num_par
		trace = trace + rho_squ(ii,ii)
	END DO
	
	WRITE(*,*) "Checking trace of rho (separable case)", trace


	! Checking the density matrix is hermitian
	Herm_flag = .TRUE.	
	DO ii=1,dim_**num_par
		DO jj=ii,dim_**num_par
			IF (rho_sep(ii,jj).NE.CONJG(rho_sep(jj,ii))) THEN
				Herm_flag = .FALSE.
			END IF
		END DO
	END DO

	WRITE(*,*) "Checking matrix is hermitian", Herm_flag

! ---------- Evaluating the partial Trace ----------

! For case N=2, evaluating the right and the left trace
!	For N=2 (dim_ = d, with d generic integer >= 0)
!	Matrices are d*d, so the system density matrix is (d^2)^N
!
!	When tracing out a subsystem (one at a time) the dimension becomes: (d^2)^(N-1)

	red_dim = dim_**(num_par-1)

	ALLOCATE(rho_red(red_dim,red_dim))
	
	! Initialize rho_red to 0s
	DO ii=1,red_dim
		DO jj=1,red_dim
			rho_red(ii,jj) = CMPLX(0d0,0d0)
		END DO
	END DO

	! 1 <= Particle <= num_par: selecting 1 particle
	!	In our case, the last particle in psi_sep_temp is that changing factor at every "step" in the matrix
	particle = 1

	
	! Looping over the reduced matrix, building it
	DO ii=0,dim_**(num_par-1)-1
		DO jj=0,dim_**(num_par-1)-1
			! Building the reduced matrix directly
			DO ii_red=0,dim_-1
				index_i = ii*dim_**(particle-1)+ii_red*dim_**(num_par-particle) + 1
				index_j = jj*dim_**(particle-1)+ii_red*dim_**(num_par-particle) + 1
				rho_red(ii+1,jj+1) = rho_red(ii+1,jj+1) + rho_sep(index_i,index_j)
			END DO
		END DO
	END DO
	
	WRITE(*,*) "Printing evaluated reduced matrix"
	DO ii=1,dim_**(num_par-1)
		DO jj=1,dim_**(num_par-1)
			WRITE(*,*) rho_red(ii,jj)
		END DO
	END DO

	! Writing expected rho_red (in N=2 case it corresponds to the density matrix of the other particle)
	WRITE(*,*) "Original reduced"
	DO ii=1,dim_
		DO jj=1,dim_
			WRITE(*,*) psi_sep_temp((num_par-particle)*dim_+jj)*CONJG(psi_sep_temp((num_par-particle)*dim_+ii))
		END DO
	END DO

! ---------- WRITING RESULTS TO FILE ----------

	! Writing rho_sep to file
	OPEN(12, file="rho_sep.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dim_**num_par
		WRITE(12,"(F11.8,F13.8)",advance="No") rho_sep(ii,1)

		DO jj=2,(dim_**(num_par-1)-1)
			WRITE(12,"(F15.8,F13.8)",advance="No") rho_sep(ii, jj)
		END DO
		
		WRITE(12,"(F15.8,F13.8)",advance="Yes") rho_sep(ii,dim_**(num_par-1))
	END DO
	
	CLOSE(12)


	! Writing evaluated rho_red to file
	OPEN(12, file="rho_red_eval.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dim_**(num_par-1)
		WRITE(12,"(F11.8,F13.8)",advance="No") rho_red(ii,1)

		DO jj=2,(dim_**(num_par-1)-1)
			WRITE(12,"(F15.8,F13.8)",advance="No") rho_red(ii, jj)
		END DO
		
		WRITE(12,"(F15.8,F13.8)",advance="Yes") rho_red(ii,dim_**(num_par-1))
	END DO
	
	CLOSE(12)


	! Writing expected rho_red to file (in N=2 case it corresponds to the density matrix of the other particle)
	OPEN(12, file="rho_red_expec.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dim_**(num_par-1)
		WRITE(12,"(F11.8,F13.8)",advance="No") psi_sep_temp((num_par-particle)*dim_+1)*CONJG(psi_sep_temp((num_par-particle)*dim_+ii))

		DO jj=2,(dim_**(num_par-1)-1)
			WRITE(12,"(F15.8,F13.8)",advance="No") psi_sep_temp((num_par-particle)*dim_+jj)*CONJG(psi_sep_temp((num_par-particle)*dim_+ii))
		END DO
		
		WRITE(12,"(F15.8,F13.8)",advance="Yes") psi_sep_temp((num_par-particle)*dim_+dim_)*CONJG(psi_sep_temp((num_par-particle)*dim_+ii))
	END DO
	
	CLOSE(12)
! ---------- DEALLOCATING ----------
	
	DEALLOCATE(psi_sep)
	DEALLOCATE(psi_not_sep)
	DEALLOCATE(psi_sep_temp)
	DEALLOCATE(rho_not_sep)
	DEALLOCATE(rho_sep)
	DEALLOCATE(rho_red)
	
	STOP
	
	END PROGRAM Ex8

