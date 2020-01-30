	PROGRAM Ex8
	! Program that studies random hermitian matrices
	! Variables are:
	! 

	USE DEBUG

	IMPLICIT NONE

	REAL*8 :: coeff, real_temp, comp_temp, trace
	COMPLEX*16 :: const, trace_app
	COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi_sep_temp, psi_sep, psi_not_sep
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: rho_not_sep, trace_rho
	COMPLEX*16, DIMENSION(:), ALLOCATABLE :: app_red, rho_red, rho_sep
	INTEGER*8 :: num_par, dim_, ii, jj, ii_min, ii_max, temp_r, steps, my_r
	INTEGER*8 :: red_dim, particle, ii_red, jj_red, jump, max_index, step_par
	INTEGER*8 :: num_oper, app_index, n_iter
	LOGICAL :: Debug_flag
	!REAL*8 :: h_bar, mass, omega, pi_
	
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
	ALLOCATE(rho_sep((dim_**num_par)*(dim_**num_par)))
	
	ALLOCATE(trace_rho(dim_**num_par,dim_**num_par))



! ---------- DOING TASK ----------

! ---------- Wave functions ----------

	! Initializing NOT SEPARABLE wave function vector
	CALL Checkpoint(Debug_flag, 'Read', "psi_not_sep", psi_not_sep, dim_**num_par)
	DO ii=1,dim_**num_par
		CALL RANDOM_NUMBER(real_temp)
		CALL RANDOM_NUMBER(comp_temp)
		psi_not_sep(ii) = CMPLX(real_temp,comp_temp)
	END DO
	
	! Normalizing
	psi_not_sep(:) = psi_not_sep(:) / SQRT(SUM(ZABS(psi_not_sep)**2))
	
	CALL Checkpoint(Debug_flag, 'Write', "psi_not_sep", psi_not_sep, dim_**num_par)

	
	! Initializing SEPARABLE wave function vector (initializing only "single body wave functions")
	CALL Checkpoint(Debug_flag, 'Read', "psi_sep_temp", psi_sep_temp, dim_*num_par)
	DO ii=0,num_par-1
		DO jj=1,dim_
			CALL RANDOM_NUMBER(real_temp)
			CALL RANDOM_NUMBER(comp_temp)
			psi_sep_temp(ii*dim_+jj) = CMPLX(real_temp,comp_temp)
		END DO
		ii_min = ii*dim_+1
		ii_max = ii*num_par+jj
		! Normalizing every single particle wave function
		psi_sep_temp(ii_min:ii_max) = psi_sep_temp(ii_min:ii_max) / SQRT(SUM(ZABS(psi_sep_temp(ii_min:ii_max)**2)))
	END DO
	CALL Checkpoint(Debug_flag, 'Write', "psi_sep_temp", psi_sep_temp, dim_*num_par)

! ---------- Density matrices ----------

	! NOT SEPARABLE density matrix
	DO ii=1,dim_**num_par
		DO jj=1,dim_**num_par
			rho_not_sep(ii,jj) = psi_not_sep(ii)*CONJG(psi_not_sep(jj))
		END DO
	END DO



	
	! Separable case
	! Bulding psi (the complete one)
	DO ii=1,dim_**num_par
		steps = 0
		! The "code" should give the result "0" (all indexes to 0) at beginning, not at the end
		!	so here one has to start from 0, not from 1
		temp_r = ii-1
		const = (1.,0.)
		
		
		DO jj=1,num_par
			my_r = MOD(temp_r,dim_)
			
			! Evaluating the product using this trick
			const = const * psi_sep_temp(dim_*(num_par - jj) + my_r + 1)
			! Updating the number to divide
			temp_r = (temp_r-my_r)/dim_
		END DO
		
		psi_sep(ii) = const

	END DO
	
	
	! SEPARABLE density matrix
	DO ii=1,dim_**num_par
		DO jj=1,dim_**num_par
			rho_sep((ii-1)*dim_**num_par+jj) = psi_sep(jj)*CONJG(psi_sep(ii))
		END DO
	END DO
	
	trace = 0
	!trace_rho = MATMUL(rho_sep, rho_sep)
	
	!DO ii=1,dim_**num_par
	!	trace = trace + trace_rho(ii,ii)
	!END DO
	
	!WRITE(*,*) "Checking track of rho", trace


! ---------- Evaluating the partial Trace ----------

! For case N=2, evaluating the right and the left trace
!	For N=2 (dim_ = d, with d generic integer >= 0)
!	Matrices are d*d, so the system density matrix is (d^2)^N
!
!	When tracing out a subsystem (one at a time) the dimension becomes: (d^2)^(N-1)

	red_dim = dim_**(num_par-1)
	
	ALLOCATE(app_red(dim_*dim_))
	ALLOCATE(rho_red(red_dim*red_dim))

	! 1 <= Particle <= num_par: selecting 1 particle
	!	In our case, the last particle in psi_sep_temp is that changing factor at every "step" in the matrix
	particle = 1
	
	! The max value of the index to explore the matrix
	max_index = dim_**(2*(num_par + 1 - particle))
	
	! Evaluating the step to find out the matrices for a particle
	step_par = max_index / dim_**2
	
	! Evaluating the number of operations to do
	num_oper = dim_**(2*num_par) / max_index
	
	WRITE(*,*) "DUMPING RHO_SEP"
	WRITE(*,*) "Printing product"
	DO ii=1,dim_**num_par
		DO jj=1,dim_**num_par
			WRITE(*,*) ii,jj,psi_sep(jj)*CONJG(psi_sep(ii))
		END DO
	END DO
	
! ---------- WRONG ALGORITHM ---------- (see notes for details)
	
	! Looping over the number of iterations (have to evaluate dim_**(2*(num_par-1)) terms to insert in reduced matrix)
	DO n_iter = 1, dim_**(2*(num_par-1))
		! Getting the submatrices (ii index in this way should do only dim_*dim_ steps, so it fits in the "app_red matrix")
		app_index = 1
		
		! When a cycle is done, go to next
		!	When n_iter has reached step_par, then the process jumps to next interval (if any)
		IF (((n_iter-1).GT.0).AND.(MOD((n_iter-1),step_par).EQ.0)) THEN
			jump = jump + 1
		END IF
		
		
		DO ii=max_index*(jump)+n_iter,max_index*(1+jump),step_par
			app_red(app_index) = rho_sep(ii)
			app_index = app_index + 1
			WRITE(*,*) "Element", app_red(app_index-1)
		END DO

	
		! Evaluating trace and storing it
		trace_app = (0.,0.)
		DO jj = 1,dim_*dim_,(dim_+1)
			trace_app = trace_app + app_red(jj)
		END DO
		
		WRITE(*,*) "Temporary trace", trace_app
			
		rho_red(n_iter) = trace_app
	END DO
	
	WRITE(*,*) "Evaluated reduced"
	DO ii=1,red_dim*red_dim
		WRITE(*,*) rho_red(ii)
	END DO
	
	WRITE(*,*) "Original reduced"
	DO ii=1,dim_
		DO jj=1,dim_
			WRITE(*,*) psi_sep_temp((particle-1)*dim_+jj)*CONJG(psi_sep_temp((particle-1)*dim_+ii))
		END DO
	END DO
	
	WRITE(*,*) "Original reduced"
	DO ii=1,dim_
		DO jj=1,dim_
			WRITE(*,*) psi_sep_temp((particle)*dim_+jj)*CONJG(psi_sep_temp((particle)*dim_+ii))
		END DO
	END DO
	
	WRITE(*,*) "Printing product"
	DO ii=1,dim_**num_par
		DO jj=1,dim_**num_par
			WRITE(*,*) ii,jj,psi_sep(jj)*CONJG(psi_sep(ii))
		END DO
	END DO


! ---------- DEALLOCATING ----------
	
	DEALLOCATE(psi_sep)
	DEALLOCATE(psi_not_sep)
	
	DEALLOCATE(psi_sep_temp)
	
	DEALLOCATE(rho_not_sep)
	DEALLOCATE(rho_sep)
	
	STOP
	
	END PROGRAM Ex8

