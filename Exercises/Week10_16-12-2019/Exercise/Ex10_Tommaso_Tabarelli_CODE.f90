	PROGRAM Ex10
	! Program that computes the ground state of a Hamiltonian of N particles in a 1-D lattice using real-space RG
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
	
	INTEGER*8 :: num_par, dim_, ii, jj, kk, hh, step, size_n
	! Using integer*4 vector to evaluate hamiltonian diagonal faster
	INTEGER*4, DIMENSION(:), ALLOCATABLE :: h_diag
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt_magn
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt_int
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: app
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt_l, hamilt_r
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: hamilt_2N, eig_vec
	COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: identity_n
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

	! Initializing size_n
	size_n = dim_**num_par
	
! ---------- ALLOCATING VARIABLES ----------

	ALLOCATE(hamilt(dim_**num_par,dim_**num_par))
	ALLOCATE(hamilt_magn(dim_**num_par,dim_**num_par))
	ALLOCATE(hamilt_int(dim_**num_par,dim_**num_par))
	
	ALLOCATE(app(dim_**num_par,dim_**num_par))

	ALLOCATE(hamilt_l(dim_**(num_par),dim_**(num_par)))
	ALLOCATE(hamilt_r(dim_**(num_par),dim_**(num_par)))
	
	ALLOCATE(hamilt_2N(dim_**(2*num_par),dim_**(2*num_par)))
	ALLOCATE(eig_vec(dim_**(2*num_par),dim_**(2*num_par)))
	!ALLOCATE(hamilt_lr())
	ALLOCATE(identity_n(dim_**(num_par),dim_**(num_par)))
	
	ALLOCATE(eig_val(dim_**(2*num_par)))
	
	ALLOCATE(h_diag(dim_**num_par))
	
	ALLOCATE(RWORK_(3*dim_**(2*num_par)-2))


! ---------- EVALUATING HAMILTONIAN FOR THE "SMALL" SYSTEM ----------

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

	!DO ii=1,dim_**num_par
	!	WRITE(*,"(F5.2,F5.2)", advance="no") hamilt_magn(ii,:)
	!END DO

	! Checking the diagonal
	DO ii=1,dim_**num_par
		WRITE(*,*) h_diag(ii)
	END DO


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

	! Checking interaction hamiltonian
	!WRITE(*,*) "Hamilt_int"	
	!DO ii=1,dim_**num_par
	!	WRITE(*,"(F5.2,F5.2)", advance="no") hamilt_int(ii,1)
	!	DO jj=2,dim_**num_par-1
	!		WRITE(*,"(F5.2,F5.2)", advance="no") hamilt_int(ii,jj)
	!	END DO
	!	WRITE(*,"(F5.2,F5.2)") hamilt_int(ii,dim_**num_par)
	!END DO

! Initializing identity N

	DO ii=1,dim_**num_par
		identity_n(ii,ii) = 1
	END DO


! ---------- LOOPING OVER LAMBDAS ----------

	DO ii=0,10

		lambda = 3.0/(10.)*ii

		WRITE(*,*) "Evaluating lambda = ", lambda

		! Initializing system hamiltonian (for every lambda)
		hamilt = hamilt_int + lambda*hamilt_magn

		! Initializing hamilt_l and hamilt_r (for every lambda)
		hamilt_r = 0
		hamilt_l = 0
		
		DO jj=0,num_par-1
		
			! First step, initialize stuff
			IF (jj.EQ.0) THEN
				hamilt_r(1:dim_**(jj+1),1:dim_**(jj+1)) = sigma_x(:,:)
				hamilt_l(1:dim_**(jj+1),1:dim_**(jj+1)) = ident_2(:,:)
			ELSE
				! Conditions for hamilt_l
				IF (jj.EQ.num_par-1) THEN
					hamilt_r(1:dim_**(jj+1),1:dim_**(jj+1)) = Tensor_product(hamilt_r(1:dim_**jj,1:dim_**jj), ident_2(:,:))
					hamilt_l(1:dim_**(jj+1),1:dim_**(jj+1)) = Tensor_product(hamilt_l(1:dim_**jj,1:dim_**jj), sigma_x(:,:))
				ELSE
					! Enlarging using identities
					hamilt_r(1:dim_**(jj+1),1:dim_**(jj+1)) = Tensor_product(hamilt_r(1:dim_**jj,1:dim_**jj), ident_2(:,:))
					hamilt_l(1:dim_**(jj+1),1:dim_**(jj+1)) = Tensor_product(hamilt_l(1:dim_**jj,1:dim_**jj), ident_2(:,:))
				END IF
			END IF
		END DO
		
		! ---------- RG PROCEDURE ----------
		! Iterating 50 times (conventionally)
		DO hh=1,50
		
			hamilt_2N = Tensor_product(hamilt,identity_n) + Tensor_product(identity_n,hamilt) + Tensor_product(hamilt_l,hamilt_r)
			
			! Storing another matrix in which there will be eigenvectors
			eig_vec = hamilt_2N

			! DIAGONALIZING
			! Allocating 1 value for wk_opt
			ALLOCATE(wk_opt(1))
			
			! Giving LWORK=-1 the subroutine only evaluates the optimal dimension for WORK_
			LWORK_=-1
			
			! See documentation online
			!  (http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen.html) 
			CALL ZHEEV('V', 'U', dim_**(2*num_par), eig_vec, dim_**(2*num_par), eig_val, wk_opt, LWORK_, RWORK_, info_eig)
			
			LWORK_ = INT(wk_opt(1))
			
			! Allocating the optimal dimension for WORK
			ALLOCATE(WORK_(LWORK_))

			! Recalling the subroutine to do the diagonalization	
			CALL ZHEEV('V', 'U', dim_**(2*num_par), eig_vec, dim_**(2*num_par), eig_val, WORK_, LWORK_, RWORK_, info_eig)
			
			DEALLOCATE(wk_opt)
			DEALLOCATE(WORK_)
			
			! Updating hamilt using truncation (first dim_**num_par eigenvectors)
			!	(hamilt_2N is [dim_**(2*num_par), dim_**(2*num_par)] )
			hamilt = MATMUL( TRANSPOSE(CONJG(eig_vec(:,1:size_n))), MATMUL(hamilt_2N, eig_vec(:,1:size_n)) )

			! Updating hamilt_l and hamilt_r (first N eigenvectors)
			hamilt_l = MATMUL( TRANSPOSE(CONJG(eig_vec(:,1:size_n))), MATMUL(Tensor_product(hamilt_l,identity_n),eig_vec(:,1:size_n)) )
			hamilt_r = MATMUL( TRANSPOSE(CONJG(eig_vec(:,1:size_n))), MATMUL(Tensor_product(identity_n,hamilt_r),eig_vec(:,1:size_n)) )
		
		END DO
	
	
		! Storing the first four eigenvalues for every lambda
		eig_val_results(ii+1,1) = lambda
		eig_val_results(ii+1,2) = eig_val(1)/(num_par*2**50)

		WRITE(*,*) "Ground state energy:", eig_val(1)/(num_par*2**50)
		WRITE(*,*) ""
	END DO


! ---------- STORING RESULTS ----------

	! Saving eigenvalues to file
	OPEN(12, file="eig_val_results.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,11
		WRITE(12,"(F12.8)",advance="No") eig_val_results(ii,1)
		WRITE(12,"(F15.8)",advance="Yes") eig_val_results(ii,2)
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

	END PROGRAM Ex10
	
	
	! Checking hamiltonian writing it
	!WRITE(*,*) "Hamilt"	
	!DO hh=1,dim_**num_par
	!	WRITE(*,"(F5.2,F5.2)", advance="no") hamilt(hh,1)
	!	DO jj=2,dim_**num_par-1
	!		WRITE(*,"(F5.2,F5.2)", advance="no") hamilt(hh,jj)
	!	END DO
	!	WRITE(*,"(F5.2,F5.2)") hamilt(hh,dim_**num_par)
	!END DO
