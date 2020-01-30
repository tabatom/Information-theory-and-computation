	PROGRAM Ex9
	! Program that helps in studying a system made of N spin-1/2 particles in a 1-D lattice
	! Variables are:
	!	num_par		: it is the number of particles. Parameters read from file "num_par.txt"
	!	dim_			: it is the dimension of every subsystem. Parameters read from file "dim_.txt"
	!	ii,jj,hh,kk		: simply some numbers to do loops
	!	step			: it is the step used to build the diagonal non-interacting terms
	!	h_diag		: it is a vector used then to initialize matrix diagonal
	!	hamilt		: matrix representing the hamiltonian of the system
	!	hamilt_magn		: matrix representing the non-interaction part of the hamiltonian
	!	hamilt_int		: matrix representing the interaction part of the hamiltonian
	!	app			: matrix used as support to build the hamiltonian
	!	eig_val		: it is the eigenvalues vector used to collect eigenvalues from the llapack subroutine
	!	eig_val_results	: it is a matrix used to store the first four eigenvectors of every iteration
	!	lambda		: variable to change non-interaction contribution to the hamiltonian
	!	Debug_flag		: logical value to enable runtime error warnings and checkpointing
	!	counter		: integer variable to be used to count checkpoints steps

	USE DEBUG
	USE MY_MATH

	IMPLICIT NONE
	
	INTEGER*8 :: num_par, dim_, ii, jj, kk, hh, step
	! Using integer*4 vector to evaluate hamiltonian diagonal faster
	INTEGER*4, DIMENSION(:), ALLOCATABLE :: h_diag
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt_magn
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt_int
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: app
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eig_val
	REAL*8, DIMENSION(21,5) :: eig_val_results
	REAL*8 :: lambda
	! Defining 2-d identity and Pauli's matrices
	COMPLEX*16, DIMENSION(2,2), PARAMETER ::&
		ident_2 = RESHAPE((/ &
		CMPLX(1d0,0d0), CMPLX(0d0,0d0),&
		CMPLX(0d0,0d0), CMPLX(1d0,0d0)/),&
		SHAPE(ident_2)),&
		sigma_x = (RESHAPE((/ &
		CMPLX(0d0,0d0), CMPLX(1d0,0d0),&
		CMPLX(1d0,0d0), CMPLX(0d0,0d0)/),&
		SHAPE(sigma_x))),&
		sigma_y = (RESHAPE((/ &
		CMPLX(0d0,0d0), CMPLX(0d0,1d0),&
		CMPLX(0d0,-1d0), CMPLX(0d0,0d0)/),&
		SHAPE(sigma_y))),&
		sigma_z = (RESHAPE((/ &
		CMPLX(1d0,0d0), CMPLX(0d0,0d0),&
		CMPLX(0d0,0d0), CMPLX(-1d0,0d0) /),&
		SHAPE(sigma_z)))
	

	LOGICAL :: Debug_flag
	INTEGER*8 :: counter
	
	! Declearing variables to be used in DYSEV subroutine
	COMPLEX*16, DIMENSION(:), ALLOCATABLE ::  WORK_, wk_opt
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK_
	INTEGER*8 :: LWORK_, info_eig


	COMMON counter
	counter = 0

	! Debugging
	Debug_flag = .FALSE.


	! Opening file to read number particles
	OPEN(11, file="num_par.txt", status='old', action='READ')
	! Reading number of space intervals from file
	READ(11,*) num_par
	! Closing file
	CLOSE(11)
	
	! Opening file to read dimension of subspaces
	OPEN(11, file="dim_.txt", status='old', action='READ')
	! Reading dim from file
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

	ALLOCATE(hamilt(dim_**num_par,dim_**num_par))
	ALLOCATE(hamilt_magn(dim_**num_par,dim_**num_par))
	ALLOCATE(hamilt_int(dim_**num_par,dim_**num_par))
	
	ALLOCATE(app(dim_**num_par,dim_**num_par))
	
	ALLOCATE(eig_val(dim_**num_par))
	
	ALLOCATE(h_diag(dim_**num_par))
	
	ALLOCATE(RWORK_(3*dim_**num_par-2))


! ---------- DOING TASK ----------

! ---------- Evaluating diagonal terms ----------

	! Initializing diagonal to 0s
	h_diag = 0

	! Initializing hamiltonian to 0s
	hamilt(:,:) = 0
	hamilt_magn(:,:) = 0
	hamilt_int(:,:) = 0


	! Initializing app matrix to 0s
	app = 0

	! Building diagonal
	DO ii=1,num_par
		step = dim_**(num_par+1-ii)
		DO jj=0,dim_**num_par-1,step
			h_diag((jj+1):(jj+step/2)) = h_diag((jj+1):(jj+step/2))+1
			h_diag((jj+step/2+1):(jj+step)) = h_diag((jj+step/2+1):(jj+step))-1
		END DO
	END DO

	! Initializing hamiltonian of sigma_z contribute using the diagonal
	!	representing the interaction
	DO ii=1,dim_**num_par
		hamilt_magn(ii,ii) = h_diag(ii)
	END DO

	! Checking the diagonal
	DO ii=1,dim_**num_par
		WRITE(*,*) h_diag(ii)
	END DO

	!DO ii=1,dim_**num_par
	!	WRITE(*,"(F5.2,F5.2)", advance="no") hamilt(ii,:)
	!END DO


! ---------- Evaluating interaction terms ----------

	! Looping over the interaction couples (which are N-1)
	DO ii=1,num_par-1
		! Restoring the app variable
		app = 0
		! We have to evaluate only (N-1) tensor products (in case of dim^2, we have only 1 tensor product for example)
		!	The trick is to evaluate N tensor products, having the first one trivially evaluated
		DO jj=1,num_par
			! If first step, then do "trivial evaluation"
			IF (jj.EQ.1) THEN
				! The interactions are between particle ii and ii+1
				IF (((num_par+1-jj).EQ.ii).OR.((num_par+1-jj).EQ.(ii+1))) THEN
					app(1:dim_**(jj),1:dim_**(jj)) = sigma_x(:,:)
					! Debugging
					WRITE(*,*) "sigma_x"
				ELSE
					app(1:dim_**(jj),1:dim_**(jj)) = ident_2(:,:)
					! Debugging
					WRITE(*,*) "ident_2"
				END IF
			ELSE
				! The interactions are between particle ii and ii+1
				IF (((num_par+1-jj).EQ.ii).OR.((num_par+1-jj).EQ.(ii+1))) THEN
					! Storing increasing matrices to evaluate the temporary results
					app(1:dim_**(jj),1:dim_**(jj)) = Tensor_product(sigma_x, app( 1:dim_**(jj-1),1:dim_**(jj-1) ))
					! Debugging
					WRITE(*,*) "sigma_x"
				ELSE
					app(1:dim_**(jj),1:dim_**(jj)) = Tensor_product(ident_2, app( 1:dim_**(jj-1),1:dim_**(jj-1) ))
					! Debugging
					WRITE(*,*) "ident_2"
				END IF
			END IF
		END DO
		
		WRITE(*,*) "--------------------------"
		
		! After having evaluated the temporary result, add it to the hamiltonian
		hamilt_int(:,:) = hamilt_int(:,:) + app(:,:)
		
	END DO
	

! Checking HAMILTONIAN results writing to file
	!OPEN(12, file="hamilt.txt", status='REPLACE', action='WRITE')
	!
	!DO ii=1,dim_**(num_par)
	!	WRITE(12,"(F5.2,F7.2)",advance="No") REAL(hamilt_magn(ii,1)+hamilt_int(ii,1))
	!
	!	DO jj=2,(dim_**(num_par)-1)
	!		WRITE(12,"(F7.2,F7.2)",advance="No") REAL(hamilt_magn(ii,jj)+hamilt_int(ii,jj))
	!	END DO
	!	
	!	WRITE(12,"(F7.2,F7.2)",advance="Yes") REAL(hamilt_magn(ii,dim_**num_par)+hamilt_int(ii,dim_**num_par))
	!END DO
	!
	!CLOSE(12)

	
! ---------- LOOPING OVER LAMBDAS ----------

	DO ii=0,20

		lambda = 3.0/(20.)*ii

		WRITE(*,*) "Evaluating lambda = ", lambda

		! Building proper hamiltonian
		hamilt = 0
		hamilt(:,:) =  hamilt_int(:,:) + lambda*hamilt_magn(:,:)

		! DIAGONALIZING
		! Allocating 1 value for wk_opt
		ALLOCATE(wk_opt(1))
		
		! Giving LWORK=-1 the subroutine only evaluates the optimal dimension for WORK_
		LWORK_=-1
		
		! See documentation online
		!  (http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen.html) 
		CALL ZHEEV('N', 'U', dim_**num_par, hamilt, dim_**num_par, eig_val, wk_opt, LWORK_, RWORK_, info_eig)
		
		LWORK_ = INT(wk_opt(1))
		
		! Allocating the optimal dimension for WORK
		ALLOCATE(WORK_(LWORK_))

		! Recalling the subroutine to do the diagonalization	
		CALL ZHEEV('N', 'U', dim_**num_par, hamilt, dim_**num_par, eig_val, WORK_, LWORK_, RWORK_, info_eig)
		
		DEALLOCATE(wk_opt)
		DEALLOCATE(WORK_)
	
		! Storing the first four eigenvalues for every lambda
		eig_val_results(ii+1,1) = lambda
		eig_val_results(ii+1,2:5) = eig_val(1:4)

	END DO

! ---------- STORING RESULTS ----------

	! Saving eigenvalues to file
	OPEN(12, file="eig_val_results.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,21
		WRITE(12,"(F12.8)",advance="No") eig_val_results(ii,1)

		DO jj=2,5-1
			WRITE(12,"(F15.8)",advance="No") eig_val_results(ii,jj)
		END DO
		
		WRITE(12,"(F15.8)",advance="Yes") eig_val_results(ii,5)
	END DO
	
	CLOSE(12)
	

! ---------- DEALLOCATING ----------

	DEALLOCATE(hamilt)
	DEALLOCATE(hamilt_magn)
	DEALLOCATE(hamilt_int)
	
	DEALLOCATE(app)
	
	DEALLOCATE(eig_val)
	
	DEALLOCATE(h_diag)
	
	DEALLOCATE(RWORK_)

	END PROGRAM Ex9
