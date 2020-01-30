	PROGRAM Ex7
	! Program that helps studying the development of a wave function
	! Variables are:
	! Hamilt	:	hamiltonian matrix
	! L_interv	:	the interval is [-L_interv, L_interv]
	! k_split	:	number of intervals used to split the considered interval
	! dims	:	= k_split+1, is the matrix dimension.
	! h_spac	:	= 2L_interv/k_split ---> it is the width of the interval used to discretize
	! x_pos	:	the positions in which calculous are done (x_i = -L + i*(2*L/k))
	! eig_fun_eval :	the evaluation of the theoretical eigenfunctions on the interval points
	! info_eig	:	is an output variable to be used in eigenvalues subroutine; if info_eig==0,
	!				then eigenvalues are stored in ev_exp in ascending order
	! ev_exp	:	vector of calculated eigenvalues of Hamiltonian matrix
	! ev_theor	:	vector of theoretical eigenvalues of Hamiltonian matrix
	! ii,jj	:	indeces/numbers to be used when needed
	! dims_	:	helpful variable to store both dimension matrix in a single variable
	!				and to use debug subroutine
	! Debug_flag:	logical value to enable runtime error warnings and checkpoints
	! counter	:	integer variable to be used to count checkpoints steps
	! h_bar	:	Plank constant
	! mass	:	mass of the quantistic particle
	! omega	:	=sqrt(k/mass), where k is the elastic constant (not considered in this program).
	!				omega is the classic frequence of the oscillation
	! eigf_const:	the constant to be evaluated for each eigenfunction
	! pi_		:	is the constant "greek pie" ~ 3.1415...
	! fact_const:	is the factorial normalizing value: different for every eigenfunction

	USE DEBUG
	USE HERMITE
	USE MY_MATH

	USE, INTRINSIC :: iso_c_binding 
	INCLUDE 'fftw3.f03'

	REAL*8 :: L_interv, h_spac, T_tot, delta_t, max_p, max_Vx
	REAL*8, DIMENSION(:), ALLOCATABLE :: x_pos, p_moments
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: psi0_eval, psi0_eval_t
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: psi0_squ, psi0_squ_t
	COMPLEX*16, DIMENSION(:), ALLOCATABLE :: app1, fft1, app2, fft2, app3
	INTEGER*8 :: ii, jj, k_split, time_split, dim_spac, dim_time
	LOGICAL :: Debug_flag
	REAL*8 :: h_bar, mass, omega, eigf_const, pi_, fact_const

	! Declearing variables to be used in "fft" subroutines
	TYPE(C_PTR) :: plan
	!DOUBLE COMPLEX, DIMENSION(:) :: IN_, OUT_
	
	COMMON counter
	counter = 0

	! Debugging
	Debug_flag = .FALSE.

	! Evaluating pi
	pi_ = 2.*ASIN(1.)

	! Opening file to read interval dimensions
	OPEN(11, file="L.txt", status='old', action='READ')
	! Reading interval length from file
	READ(11,*) L_interv
	! Closing file
	CLOSE(11)

	! Opening file to read number of space intervals
	OPEN(11, file="k_split.txt", status='old', action='READ')
	! Reading number of space intervals from file
	READ(11,*) k_split
	! Closing file
	CLOSE(11)
	
	! Opening file to read total time
	OPEN(11, file="T_tot.txt", status='old', action='READ')
	! Reading total time from file
	READ(11,*) T_tot
	! Closing file
	CLOSE(11)
	
	! Opening file to read number of time intervals
	OPEN(11, file="time_split.txt", status='old', action='READ')
	! Reading number of time intervals from file
	READ(11,*) time_split
	! Closing file
	CLOSE(11)
	
	
	! Getting space intervals width
	h_spac = 2.0*L_interv/k_split
	
	! Getting time intervals width
	delta_t = 1.0*T_tot/time_split
	
	! Getting vector space dimension
	dim_spac = k_split+1
	
	! Getting vector time dimension (counting also the starting instant)
	dim_time = time_split+1
	
	! Evaluating the max of frequencies
	max_p = 2.*pi_/h_spac
	
	
	! Printing all variables to debug
	!WRITE(*,*)"L_interv", L_interv
	!WRITE(*,*)"k_split", k_split
	!WRITE(*,*)"T_tot", T_tot
	!WRITE(*,*)"time_split", time_split
	!WRITE(*,*)"h_spac", h_spac
	!WRITE(*,*)"delta_t", delta_t
	!WRITE(*,*)"dim_spac", dim_spac
	!WRITE(*,*)"dim_time", dim_time
	
	! Checking dimensions greater than 0
	IF (Debug_flag.AND.(dim_spac.LE.0)) THEN
		CALL Print_error(Debug_flag, "Dim<=0")
	END IF
	IF (Debug_flag.AND.(dim_time.LE.0)) THEN
		CALL Print_error(Debug_flag, "Dim<=0")
	END IF
	
	! Defining physical constants
	h_bar = 1.0
	mass = 1.0
	omega = 1.0
	
! ---------- ALLOCATING VARIABLES ----------

	ALLOCATE(x_pos(dim_spac))
	
	ALLOCATE(psi0_eval(dim_spac,dim_time))
	ALLOCATE(psi0_eval_t(dim_spac,dim_time))
	
	ALLOCATE(psi0_squ(dim_spac,dim_time))
	ALLOCATE(psi0_squ_t(dim_spac,dim_time))
	
	ALLOCATE(p_moments(dim_spac))
	
	ALLOCATE(app1(dim_spac), app2(dim_spac), app3(dim_spac))
	ALLOCATE(fft1(dim_spac), fft2(dim_spac))
	
! ---------- INITIALIZING VARIABLES ----------

	! Initializing potential maximum position
	max_Vx = 1.0

	! Initializing x_pos vector
	DO ii=1,dim_spac
		x_pos(ii) = -L_interv + (h_spac)*(INT(ii)-1)
	END DO
	
	! Initializing vector of frequencies (p moments in this case!)
	DO ii=1,(dim_spac-1)/2
		p_moments(ii) = max_p/dim_spac * (ii-1)
	END DO
	
	DO ii=(dim_spac+1)/2, dim_spac
		p_moments(ii) = max_p/dim_spac * (ii-1) - max_p
	END DO

	
	! Needed "0" as integer*8 and not integer*4
	ii = 1
	eigf_const = (mass*omega/(h_bar))**(0.5)
	
	! Evaluating the first eigenfunction on the grid
	fact_const = FLOAT(Factorial(ii-1))
	psi0_eval(:,ii) = CMPLX((eigf_const/SQRT(pi_))**0.5/SQRT(fact_const*(2**(ii-1))), 0.)
	psi0_eval(:,ii) = psi0_eval(:,ii)*evalHerm_Poly((eigf_const*x_pos(:)), (ii-1))*EXP(-0.5*(eigf_const*x_pos(:))**2)
	
	! Evaluating also the "theoretical" counterpart (evaluating using whole t_i = ii*delta_t)
	psi0_eval_t(:,ii) = psi0_eval(:,ii)

	! STARTING TRANSFORMING
	DO ii=2,dim_time
	
		! ---------- TRANSFORMING USING "LOCAL" delta_t ----------
		
		! Applying exp(-i*V(x)*delta_t/2)
		!app1(:) = (CMPLX(1,0) - CMPLX(0., -delta_t*0.5*(x_pos(:)-ii*delta_t/T_tot)**2)) * psi0_eval(:,ii-1)
		app1(:) = EXP(CMPLX(0., -0.5*delta_t*0.5*(x_pos(:)-ii*max_Vx*delta_t/T_tot)**2)) * psi0_eval(:,ii-1)
		
		! Fourier transforming from x to p domain
		CALL dfftw_plan_dft_1d(plan, dim_spac, app1, fft1, FFTW_FORWARD, FFTW_ESTIMATE)
		CALL dfftw_execute_dft(plan, app1, fft1)
		CALL dfftw_destroy_plan(plan)
		
		! Applying exp(-i*T(p)*delta_t)
		!app2(:) = (CMPLX(1,0) - CMPLX(0., -delta_t*(p_moments(:))*0.5)) * fft1(:)
		app2(:) = EXP(CMPLX(0., -delta_t*(p_moments(:)**2)*0.5)) * fft1(:)

		! Transforming back
		!  Also, re-normalizing using dim_spac
		!  (see http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029)
		!  The evaluation is not normalized and using FFTW_FORWARD and BACKWARD makes the function be multiplied
		!   the working dimension.
		CALL dfftw_plan_dft_1d(plan, dim_spac, app2, fft2, FFTW_BACKWARD, FFTW_ESTIMATE)
		CALL dfftw_execute_dft(plan, app2, fft2)
		CALL dfftw_destroy_plan(plan)
		! Re-normalizing
		fft2(:) = fft2(:) / dim_spac
		
		! Applying again exp(-i*V(x)*delta_t/2)
		!psi0_eval(:,ii) = (CMPLX(1,0) - CMPLX(0., -delta_t*0.5*(x_pos(:)-ii*delta_t/T_tot)**2)) * fft2(:)
		psi0_eval(:,ii) = EXP(CMPLX(0., -0.5*delta_t*0.5*(x_pos(:)-ii*max_Vx*delta_t/T_tot)**2)) * fft2(:)


		! ---------- TRANSFORMING USING t_i = ii * delta_t ----------
		
		! Applying exp(-i*V(x)*delta_t/2)
		!app1(:) = (CMPLX(1,0) - CMPLX(0., -delta_t*0.5*(x_pos(:)-ii*delta_t/T_tot)**2)) * psi0_eval(:,ii-1)
		app1(:) = EXP(CMPLX(0., -0.5*delta_t*ii*0.5*(x_pos(:)-delta_t*ii*max_Vx/T_tot)**2)) * psi0_eval_t(:,1)
		
		! Fourier transforming from x to p domain
		CALL dfftw_plan_dft_1d(plan, dim_spac, app1, fft1, FFTW_FORWARD, FFTW_ESTIMATE)
		CALL dfftw_execute_dft(plan, app1, fft1)
		CALL dfftw_destroy_plan(plan)
		
		! Applying exp(-i*T(p)*delta_t)
		!app2(:) = (CMPLX(1,0) - CMPLX(0., -delta_t*(p_moments(:))*0.5)) * fft1(:)
		app2(:) = EXP(CMPLX(0., -delta_t*ii*(p_moments(:)**2)*0.5)) * fft1(:)

		! Transforming back
		!  Also, re-normalizing using dim_spac
		!  (see http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029)
		!  The evaluation is not normalized and using FFTW_FORWARD and BACKWARD makes the function be multiplied
		!   the working dimension.
		CALL dfftw_plan_dft_1d(plan, dim_spac, app2, fft2, FFTW_BACKWARD, FFTW_ESTIMATE)
		CALL dfftw_execute_dft(plan, app2, fft2)
		CALL dfftw_destroy_plan(plan)
		! Re-normalizing
		fft2(:) = fft2(:) / dim_spac
		
		! Applying again exp(-i*V(x)*delta_t/2)
		!psi0_eval(:,ii) = (CMPLX(1,0) - CMPLX(0., -delta_t*0.5*(x_pos(:)-ii*delta_t/T_tot)**2)) * fft2(:)
		psi0_eval_t(:,ii) = EXP(CMPLX(0., -0.5*delta_t*ii*0.5*(x_pos(:)-delta_t*ii*max_Vx/T_tot)**2)) * fft2(:)


	END DO

	! Evaluating the square of the eigenfunctions
	DO ii=1,dim_spac
		psi0_squ(ii,:) = ABS(psi0_eval(ii,:))**2
		psi0_squ_t(ii,:) = ABS(psi0_eval_t(ii,:))**2
		
	END DO

! ---------- WRITING RESULTS ONTO FILE ----------

	! Writing psi0 evolution (real and complex) (in columns, so are easier to plot using gnuplot)
	OPEN(12, file="temp_evolut.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dim_spac
		! Writing also spatial position as first column
		WRITE(12,"(F13.8)",advance="No") x_pos(ii)

		DO jj=1,dim_time-1
			WRITE(12,"(F13.8,F13.8)",advance="No") psi0_eval(ii, jj)
		END DO
		WRITE(12,"(F13.8,F13.8)",advance="Yes") psi0_eval(ii,dim_time)
	END DO
	
	CLOSE(12)

	! Writing psi0 evolution (module of wave function) (in columns, so are easier to plot using gnuplot)
	OPEN(12, file="temp_evolut_mod.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dim_spac
		! Writing also spatial position as first column
		WRITE(12,"(F13.8)",advance="No") x_pos(ii)

		DO jj=1,dim_time-1
			WRITE(12,"(F13.8)",advance="No") psi0_squ(ii, jj)
		END DO
		WRITE(12,"(F13.8)",advance="Yes") psi0_squ(ii,dim_time)
	END DO
	
	CLOSE(12)
	

	! Writing psi0_theor evolution (module of wave function) (in columns, so are easier to plot using gnuplot)
	OPEN(12, file="temp_evolut_mod_theo.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dim_spac
		! Writing also spatial position as first column
		WRITE(12,"(F13.8)",advance="No") x_pos(ii)

		DO jj=1,dim_time-1
			WRITE(12,"(F13.8)",advance="No") psi0_squ_t(ii, jj)
		END DO
		WRITE(12,"(F13.8)",advance="Yes") psi0_squ_t(ii,dim_time)
	END DO
	
	CLOSE(12)
	
	! Saving the potential in funtciont of time
	OPEN(12, file="pot_evolut.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dim_spac
		! Writing also spatial position as first column
		WRITE(12,"(F13.8)",advance="No") x_pos(ii)
		DO jj = 1,dim_time-1
			WRITE(12,"(F13.8)",advance="No") (x_pos(ii)-delta_t*jj*max_Vx/T_tot)**2
		END DO
		WRITE(12,"(F13.8)",advance="Yes") (x_pos(ii)-delta_t*dim_time*max_Vx/T_tot)**2
	END DO
	
	CLOSE(12)


! ---------- DEALLOCATING ----------
	
	! Deallocating rises errors
	DEALLOCATE(x_pos)
	DEALLOCATE(psi0_eval)
	DEALLOCATE(psi0_eval_t)
	DEALLOCATE(psi0_squ)
	DEALLOCATE(psi0_squ_t)
	DEALLOCATE(app1)
	DEALLOCATE(app2)
	DEALLOCATE(fft1)
	DEALLOCATE(fft2)
	
	STOP
	
	END PROGRAM Ex7

