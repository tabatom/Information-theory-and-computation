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
	INTEGER*8 :: ii,jj
	INTEGER*8, DIMENSION(2) :: dims_
	LOGICAL :: Debug_flag
	REAL*8 :: real_, imag_
	REAL*8, DIMENSION(:), ALLOCATABLE :: delta_eig, s_i
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: s_i_a
	
	! Declearing dummy variables to be used in CHEEV
	COMPLEX*16, DIMENSION(:), ALLOCATABLE ::  WORK_
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK_
	INTEGER*8, DIMENSION(7) :: spacings
	REAL*8 :: temp_res
	! ASSUMING NAMES FOR SPACING FILES ARE NOT WITH MORE THAN 4 DIGITS
	CHARACTER(4) :: file_name

	COMMON counter
	counter = 0

	! Debugging
	Debug_flag = .FALSE.

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


	! Defining spacings (it is intended to be fixed)
	spacings = [dims/200, dims/100, dims/50, dims/25, dims/10, dims/5, dims/2]

	! Writing "spacings" (intervals) to file
	OPEN(11, file="spacings.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,SIZE(spacings)
		WRITE(11,"(I6.5)")spacings(ii)
	END DO
	
	CLOSE(11)

! ---------- ALLOCATING VARIABLES ----------

	! Allocating SQUARE matrix using file
	ALLOCATE(r_h_m(dims,dims))
	
	! Allocating eigenvalues array
	ALLOCATE(eig_val(dims))
	
	! Allocating eig differeces and spacings array
	ALLOCATE(delta_eig(dims-1), s_i(dims-1))
	
	! Allocating matrix for alternatives spacings
	ALLOCATE(s_i_a(SIZE(spacings), dims-1))


! ---------- INITIALIZING MATRICES (with random numbers) ----------

	! r_h_m

	! Starting checkpointing (and eventually initializing to 0 all elements)
	CALL Checkpoint(Debug_flag, 'Read', "r_h_m", r_h_m, dims_)

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
	
	! ZHEEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)
	! See documentation online
	CALL ZHEEV('N', 'U', dims, r_h_m, dims, eig_val, WORK_, 2*dims-1, RWORK_, info_eig)
	
	
	! Evaluating differences between consecutive eigenvalues
	delta_eig(:) = eig_val(2:) - eig_val(:dims-1)
	
	s_i = delta_eig/(SUM(delta_eig)/(dims-1))
	

	! Evaluating different spacings
	DO ii=1,SIZE(spacings)

		! Looping over all possible differences [ Δλ_i = λ_(i+1) - λ_i ]
		DO jj=1, (dims-1)
			
			CALL Help_spacings(spacings, ii, eig_val, jj, dims, temp_res)
			s_i_a(ii,jj) = temp_res
			
		END DO
	END DO


! ---------- SAVING RESULTS INTO A FILE ----------
	
	! Opening and REWRITING existing file
	
	! Writing "standard" results
	OPEN(11, file="s_i_distr_comp.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,(dims-1)
		WRITE(11,"(F13.7)")s_i(ii)
	END DO
	
	CLOSE(11)
	
	
	! Writing "local spacing" results
	DO ii=1,SIZE(spacings)
		! Converting spacings(ii) to string
		WRITE(file_name,"(I4.4)") spacings(ii)
		
		OPEN(11, file="s_i_distr_comp_local"//TRIM("_"//file_name)//".txt", status='REPLACE', action='WRITE')
		
		DO jj=1,dims-1
			WRITE(11,"(F13.7)")s_i_a(ii,jj)
		END DO
		
		CLOSE(11)
	END DO
	
	!STOP

! ---------- DEALLOCATING ----------
	
	! Deallocating rises errors
	DEALLOCATE(r_h_m)
	DEALLOCATE(eig_val)
	DEALLOCATE(delta_eig)
	DEALLOCATE(s_i)
	DEALLOCATE(s_i_a)
	STOP
	
	END PROGRAM Ex5


	
! ---------------------------------------------------------------
! ---------------------------------------------------------------

! ---------- SUBROUTINE TO HELP IN EVALUATING SPACINGS ----------


	! Arguments:
	! spac		: INPUT, the vector of different spacings
	! i_spac		: INPUT, the index to move into spacings vector
	! eig_values	: INPUT, the vector storing eigenvalues
	! eig_index		: INPUT, the index to move in eig_values vector
	!		N.B.	: 1 <= eig_index <= eig_num-1, where
	! eig_num 		: INPUT, is the number of eigenvalues
	! sub_result	: OUTPUT, the local spacing (already evaluated when returned by subroutine)
	
	SUBROUTINE Help_spacings(spac, i_spac, eig_values, eig_index, eig_num, sub_result)
	
		INTEGER*8, DIMENSION(*), INTENT(IN) :: spac
		INTEGER*8, INTENT(IN) :: i_spac, eig_index, eig_num
		DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: eig_values
		INTEGER*8 :: temp_spac, min_index, max_index
		REAL*8 :: temp_sum, temp_delta, temp_size, sub_result
		
		temp_spac = spac(i_spac)
		
		min_index = MAX(1,eig_index-(temp_spac/2))
		max_index = MIN(eig_index+(temp_spac/2), eig_num-1)
		
		temp_delta = eig_values(eig_index+1) - eig_values(eig_index)
		temp_sum = SUM(eig_values((min_index+1):(max_index+1)) - eig_values(min_index:max_index))
		temp_size = SIZE(eig_values((min_index+1):(max_index+1)) - eig_values(min_index:max_index))
		
		sub_result = temp_delta / (temp_sum/temp_size)
		
	END SUBROUTINE Help_spacings
