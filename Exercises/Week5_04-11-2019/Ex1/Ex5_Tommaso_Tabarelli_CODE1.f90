	PROGRAM Ex5
	! Program that study random hermitian matrices
	! Variables are:
	! r_h_m	:	random hermitian matrix
	! dims	:	matrix dimension (matrix is squared)
	! info_eig	:	is an output variable to be used in CHEEV subroutine; if info_eig==0,
	!				then eigenvalues are stored in eig_val in ascending order
	! eig_val	:	vector of eigenvalues of HERMITIAN MATRIX (they are real, not complex)
	! ii,jj	:	indeces/numbers to be used when needed
	! dims_	:	helpful variable to store both dimension matrix in a single variable
	!				and to use debug subroutine
	! Debug_flag:	logical value to enable runtime error warnings and checkpoints
	! real_	:	temporary variable used to initialize real parts of r_h_m elements
	! imag_	:	same as real, but for imaginary part
	! delta_eig	:	vector of differences between consecutives arrays [ Δλ_i = λ_(i+1) - λ_i ]
	! s_i		:	the normalized spacing between eigenvalues [ s_i = Δλ_i / AVG(Δλ_i) ]
	! s_i_a	:	a matrix built using different local spacings
	! spacings	:	list of possible spacings to be used to evaluate local AVG(Δλ_i)
	! s_i_n	:	For a generic n, average local spacing using n eig_val to evaluate
	! counter	:	integer variable to be used to count checkpoints steps

	USE DEBUG
	USE RNG
	USE SORT

	IMPLICIT NONE

	COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: r_h_m
	INTEGER*8 :: counter, dims, info_eig
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eig_val
	INTEGER*2 :: ii,jj
	INTEGER*8, DIMENSION(2) :: dims_
	LOGICAL :: Debug_flag
	REAL*8 :: real_, imag_
	REAL*8, DIMENSION(:), ALLOCATABLE :: delta_eig, s_i
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: s_i_a
	
	! Declearing dummy variables to be used in CHEEV
	COMPLEX*16, DIMENSION(:), ALLOCATABLE ::  WORK_
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK_
	INTEGER*8, DIMENSION(3) :: spacings
	REAL*8, DIMENSION(:), ALLOCATABLE :: spac200, spac100, spac50, spac25, spac10, spac5, spac2
	INTEGER*8 :: temp_spac

	COMMON counter
	counter = 0

	! Debugging
	Debug_flag = .FALSE.

	! Defining spacings (it is intended to be fixed)
	!spacings = [200, 100, 50, 25, 10, 5, 2]
	spacings = [10, 5, 2]

	! Opening file to read dimensions
	OPEN(10, file="dim_file.txt", status='old', action='READ')
	
	! Reading dimensions from file
	READ(10,*) dims
	
	! Closing file
	CLOSE(10)
	
	! Defining dims_
	dims_(1) = dims
	dims_(2) = dims
	
	! Checking dimensions greater than 0
	IF (Debug_flag.AND.(dims.LE.0)) THEN
		CALL Print_error(Debug_flag, "Dim<=0")
	END IF


! ---------- ALLOCATING VARIABLES ----------

	! Allocating SQUARE matrix using file
	ALLOCATE(r_h_m(dims,dims))
	
	! Allocating eigenvalues array
	ALLOCATE(eig_val(dims))
	
	! Allocating eig differeces and spacings array
	ALLOCATE(delta_eig(dims-1), s_i(dims-1))
	
	! Allocating matrix for alternatives spacings
	ALLOCATE(s_i_a(SIZE(spacings), dims-1))

	! Allocating correspondig vectors for spacing
	ALLOCATE(spac200(dims/200), spac100(dims/100), spac50(dims/50))
	ALLOCATE(spac25(dims/25), spac10(dims/10), spac5(dims/5), spac2(dims/2))
	

! ---------- INITIALIZING MATRICES (with random numbers) ----------

	! r_h_m

	! Starting checkpointing (and eventually initializing to 0 all elements)
	CALL Checkpoint(Debug_flag, 'Read', "r_h_m", r_h_m, dims_)
!	CALL RANDOM_NUMBER(r_h_m)

	! Making the matrix hermitian
	DO ii=1,dims
		DO jj=1,ii
			CALL Normal_rng( real_, imag_ )
			
			IF (ii.EQ.jj) THEN
				r_h_m(ii,jj) = CMPLX(real_, 0e0)
			ELSE
				r_h_m(ii,jj) = CMPLX(real_, imag_)
				r_h_m(jj,ii) = CMPLX(real_, -imag_)
			END IF
		END DO
	END DO

	CALL Checkpoint(Debug_flag, 'Write', "r_h_m", r_h_m, dims_)

	
! ---------- STARTING OPERATIONS ----------

	! Allocating variables to make subroutine works with dummy variables
	ALLOCATE(WORK_(2*dims-1), RWORK_(3*dims-2))
	
	! CHEEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)
	! Arguments are:
	!  JOBZ	: JOBZ is CHARACTER*1. If = 'N': Compute eigenvalues only; if = 'V': Compute eigenvalues and eigenvectors.
	!  UPLO	: UPLO is CHARACTER*1. If = 'U':  Upper triangle of A is stored; if = 'L': Lower triangle of A is stored.
	!  N		: N is INTEGER. The order of the matrix A.  N >= 0.
	CALL ZHEEV('N', 'L', dims, r_h_m, dims, eig_val, WORK_, 2*dims-1, RWORK_, info_eig)
	
	
	! Evaluating differences between consecutive eigenvalues
	delta_eig(:) = eig_val(2:) - eig_val(:dims-1)
	
	s_i = delta_eig/(SUM(delta_eig)/(dims-1))
	
	! Evaluating different spacings
	DO ii=1,SIZE(spacings)
		temp_spac = dims/spacings(ii)
		DO jj=1, temp_spac
			
			IF (temp_spac.EQ.(dims/200)) THEN
				spac200(jj) = SUM(eig_val(2:(spacings(ii)*jj+1)) - eig_val((spacings(ii)*(jj-1)+1):(spacings(ii)*jj)))/spacings(ii)
				! Saving new results
				s_i_a(ii,(spacings(ii)*(jj-1)+1):(spacings(ii)*jj)) = delta_eig((spacings(ii)*(jj-1)+1):(spacings(ii)*jj))/spac200(jj)
				
			ELSEIF (temp_spac.EQ.(dims/100)) THEN
				spac100(jj) = SUM(eig_val(2:spacings(ii)*jj+1) - eig_val((spacings(ii)*(jj-1)+1):(spacings(ii)*jj)))/spacings(ii)
				! Saving new results
				s_i_a(ii,(spacings(ii)*(jj-1)+1):(spacings(ii)*jj)) = delta_eig((spacings(ii)*(jj-1)+1):(spacings(ii)*jj))/spac100(jj)

			ELSEIF (temp_spac.EQ.(dims/50)) THEN
				spac50(jj) = SUM(eig_val(2:spacings(ii)*jj+1) - eig_val((spacings(ii)*(jj-1)+1):(spacings(ii)*jj)))/spacings(ii)
				! Saving new results
				s_i_a(ii,(spacings(ii)*(jj-1)+1):(spacings(ii)*jj)) = delta_eig((spacings(ii)*(jj-1)+1):(spacings(ii)*jj))/spac50(jj)

			ELSEIF (temp_spac.EQ.(dims/25)) THEN
				spac25(jj) = SUM(eig_val(2:spacings(ii)*jj+1) - eig_val((spacings(ii)*(jj-1)+1):(spacings(ii)*jj)))/spacings(ii)
				! Saving new results
				s_i_a(ii,(spacings(ii)*(jj-1)+1):(spacings(ii)*jj)) = delta_eig((spacings(ii)*(jj-1)+1):(spacings(ii)*jj))/spac25(jj)

			ELSEIF (temp_spac.EQ.(dims/10)) THEN
				spac10(jj) = SUM(eig_val(2:spacings(ii)*jj+1) - eig_val((spacings(ii)*(jj-1)+1):(spacings(ii)*jj)))/spacings(ii)
				! Saving new results
				s_i_a(ii,(spacings(ii)*(jj-1)+1):(spacings(ii)*jj)) = delta_eig((spacings(ii)*(jj-1)+1):(spacings(ii)*jj))/spac10(jj)

			ELSEIF (temp_spac.EQ.(dims/5)) THEN
				spac5(jj) = SUM(eig_val(2:spacings(ii)*jj+1) - eig_val((spacings(ii)*(jj-1)+1):(spacings(ii)*jj)))/spacings(ii)
				! Saving new results
				s_i_a(ii,(spacings(ii)*(jj-1)+1):(spacings(ii)*jj)) = delta_eig((spacings(ii)*(jj-1)+1):(spacings(ii)*jj))/spac5(jj)

			ELSEIF (temp_spac.EQ.(dims/2)) THEN
				spac2(jj) = SUM(eig_val(2:spacings(ii)*jj+1) - eig_val((spacings(ii)*(jj-1)+1):(spacings(ii)*jj)))/spacings(ii)
				! Saving new results
				s_i_a(ii,(spacings(ii)*(jj-1)+1):(spacings(ii)*jj)) = delta_eig((spacings(ii)*(jj-1)+1):(spacings(ii)*jj))/spac2(jj)
			
			END IF
		END DO
	END DO
	
	!IF (Debug_flag == .TRUE.) THEN
	!WRITE(*,*) "Info: ", info_eig
	!WRITE(*,*) eig_val
	!WRITE(*,*) delta_eig
	!WRITE(*,*) s_i_a
	!END IF
	

! ---------- SAVING RESULTS INTO A FILE ----------
	
	! Opening and REWRITING existing file
	OPEN(11, file="s_i_distr_comp.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,(dims-1)
		WRITE(11,"(F13.7)")s_i(ii)
	END DO
	
	CLOSE(11)
	
	STOP
	
! ---------- DEALLOCATING ----------
	
	! Deallocating rises errors
	DEALLOCATE(r_h_m)
	DEALLOCATE(eig_val, delta_eig, s_i, s_i_a)
	DEALLOCATE(spac200, spac100, spac50, spac25, spac10, spac5, spac2)
	STOP
	
	END PROGRAM Ex5
