	PROGRAM Ex6
	! Program that help in studying the unidimensional quantum harmonic oscillator
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

	IMPLICIT NONE

	REAL*8, DIMENSION(:,:), ALLOCATABLE :: Hamilt
	REAL*8 :: L_interv, h_spac
	INTEGER*8 :: counter, k_split, dims, info_eig
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ev_exp, ev_theor
	REAL*8, DIMENSION(:), ALLOCATABLE :: x_pos
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: eig_fun_eval
	INTEGER*8 :: ii, jj
	INTEGER*8, DIMENSION(2) :: dims_
	LOGICAL :: Debug_flag
	REAL*8 :: h_bar, mass, omega, eigf_const, pi_, fact_const
	
	! Declearing variables to be used in DYSEV subroutine
	REAL*8, DIMENSION(:), ALLOCATABLE ::  WORK_, wk_opt
	INTEGER*8 :: LWORK_

	COMMON counter
	counter = 0

	! Debugging
	Debug_flag = .FALSE.

	! Opening file to read interval dimensions
	OPEN(11, file="L.txt", status='old', action='READ')
	! Reading interval length from file
	READ(11,*) L_interv
	! Closing file
	CLOSE(11)

	! Opening file to read interval split
	OPEN(11, file="k_split.txt", status='old', action='READ')
	! Reading number of intervals from file
	READ(11,*) k_split
	! Closing file
	CLOSE(11)
	
	
	
	! Getting intervals width
	h_spac = 2.0*L_interv/k_split
	
	! Getting matrix dimension
	dims = k_split+1
	
	! Defining dims_
	dims_(1) = dims
	dims_(2) = dims

	
	! Checking dimensions greater than 0
	IF (Debug_flag.AND.(dims.LE.0)) THEN
		CALL Print_error(Debug_flag, "Dim<=0")
	END IF
	
	! Defining physical constants
	h_bar = 1.0
	mass = 1.0
	omega = 1.0

! ---------- ALLOCATING VARIABLES ----------

	! Allocating Hamiltonian matrix using file
	ALLOCATE(Hamilt(dims,dims))

	! Allocating approximated eigenvalues array
	ALLOCATE(ev_exp(dims))

	! Allocating theoretical eigenvalues array
	ALLOCATE(ev_theor(dims))
	
	! Allocating eigenfunction evaluations
	!ALLOCATE(eig_fun_eval(dims,dims))	! In principle one can evaluate all he wants
	!						!  but in practice there are recursive functions acting...
	ALLOCATE(eig_fun_eval(dims,10))

	! Allocating x_pos vector
	ALLOCATE(x_pos((dims)))

! ---------- INITIALIZING VARIABLES ----------
	
! Theoretical results:

	! Theortical eigenvalues array
	DO ii=1,dims
		ev_theor(ii) = h_bar*omega*(0.5 + (ii-1))
	END DO
	
	! Eigenfunctions (only first 10 are taken into account for now)
	
	
	
	! x_pos vector
	CALL Checkpoint(Debug_flag, 'Read', "x_pos", x_pos, dims)
	DO ii=1,dims
		x_pos(ii) = -L_interv + (h_spac)*(INT(ii)-1)
	END DO
	CALL Checkpoint(Debug_flag, 'Read', "x_pos", x_pos, dims)

	
	! Hamilt
	! Starting checkpointing
	CALL Checkpoint(Debug_flag, 'Write', "Hamilt", Hamilt, dims_)

	DO ii=1,dims
		DO jj=1,dims
			! Filling matrix diagonal
			IF (ii.EQ.jj) THEN
				Hamilt(ii,jj) = 2.0*h_bar*h_bar/(2.0*mass*h_spac*h_spac) + 0.5*mass*omega*omega*x_pos(ii)*x_pos(ii)
			ELSEIF (ii.EQ.(jj+1)) THEN
				Hamilt(ii,jj) = -1.0*h_bar*h_bar/(2.*mass*h_spac*h_spac)
			ELSEIF (ii.EQ.(jj-1)) THEN
				Hamilt(ii,jj) = -1.0*h_bar*h_bar/(2.*mass*h_spac*h_spac)
			ELSE
				Hamilt(ii,jj) = 0d0
			
			END IF
		END DO
	END DO

	CALL Checkpoint(Debug_flag, 'Write', "Hamilt", Hamilt, dims_)

	
! ---------- STARTING OPERATIONS ----------

	! Allocating 1 value for wk_opt
	ALLOCATE(wk_opt(1))
	
	! Giving LWORK=-1 the subroutine only evaluates the optimal dimension for WORK_
	LWORK_=-1
	
	! See documentation online
	!  (http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html)
	! Nowm with 
	CALL DSYEV('V', 'U', dims, Hamilt, dims, ev_exp, wk_opt, LWORK_, info_eig)
	
	LWORK_ = INT(wk_opt(1))
	
	! Allocating the optimal dimension for WORK
	ALLOCATE(WORK_(LWORK_))


	! Recalling the subroutine to do the proper task	
	CALL DSYEV('V', 'U', dims, Hamilt, dims, ev_exp, WORK_, LWORK_, info_eig)

	
	! Evaluating first 10 eigenfunctions
	! Pi = 2*arcsin(1)
	pi_ = 2*ASIN(1d0)
	eigf_const = (mass*omega/(h_bar))**(0.5)
		
	DO ii=1,10
		fact_const = FLOAT(Factorial(ii-1))
		eig_fun_eval(:,ii) = (eigf_const/SQRT(pi_))**0.5/SQRT(fact_const*(2**(ii-1)))
		eig_fun_eval(:,ii) = eig_fun_eval(:,ii)*evalHerm_Poly((eigf_const*x_pos), (ii-1))*EXP(-0.5*(eigf_const*x_pos)**2)
	END DO
	
	! Making the last sample from the half of the eigenvectors (not working...)
	!ii = INT((dims-1)/2)
	!fact_const = FLOAT(Factorial(ii-1))
	!eig_fun_eval(:,10) = (eigf_const/SQRT(pi_))**0.5/SQRT(fact_const*(2**(ii-1)))
	!eig_fun_eval(:,10) = eig_fun_eval(:,10)*evalHerm_Poly((eigf_const*x_pos), (ii-1))*EXP(-0.5*(eigf_const*x_pos)**2)
	
	! Scaling eigenvectors so that eigenfunctions are normalized
	Hamilt = Hamilt/SQRT(h_spac)


! ---------- SAVING RESULTS INTO A FILE ----------

	! Opening and REWRITING existing file
	
	! Writing theoretical eigenvalues
	OPEN(11, file="eval_theo.txt", status='REPLACE', action='WRITE')
	
	DO ii=1, dims
		WRITE(11,"(F15.5)") ev_theor(ii)
	END DO
	
	CLOSE(11)
	
	! Writing approximated eigenvalues
	OPEN(12, file="eval_exp.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dims
		WRITE(12,"(F15.5)") ev_exp(ii)
	END DO
	
	CLOSE(12)
	
	
	! Writing approximated and "NORMALIZED" (scaled so that they are normalized) eigenvectors
	!	(which are now in the Hamiltonian matrix because lapack's subroutine overwrite the matrix)
	! (writing them in columns, so are easier to plot using gnuplot)
	OPEN(12, file="evect_exp.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dims
		DO jj=1,dims-1
			WRITE(12,"(F13.8)",advance="No") Hamilt(ii,jj)
		END DO
		WRITE(12,"(F13.8)",advance="Yes") Hamilt(ii,dims)
	END DO
	
	CLOSE(12)


	! Writing theoretical and normalized eigenvectors (in columns, so are easier to plot using gnuplot)
	OPEN(12, file="evect_theo.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,dims
		DO jj=1,10-1
			WRITE(12,"(F13.8)",advance="No") eig_fun_eval(ii,jj)
		END DO
		WRITE(12,"(F13.8)",advance="Yes") eig_fun_eval(ii,jj)
	END DO
	
	CLOSE(12)

! ---------- DEALLOCATING ----------
	
	! Deallocating rises errors
	DEALLOCATE(Hamilt)
	DEALLOCATE(ev_theor)
	DEALLOCATE(ev_exp)
	DEALLOCATE(WORK_)
	DEALLOCATE(wk_opt)
	STOP
	
	END PROGRAM Ex6
