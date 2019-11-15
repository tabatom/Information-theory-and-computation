	PROGRAM Ex5_real
	! Program that generates random real diagonal matrices (it generates ONLY A VECTOR)
	! Variables are:
	! diag	:	random vector (representing a diagonal matrix)
	! dims	:	vector dimension (matrix dimension)
	! eig_val	:	vector of eigenvalues of the matrix (since it is diagonal, it is simply "diag" itself)
	! ii,jj	:	indeces/numbers to be used when needed
	! Debug_flag:	logical value to enable runtime error warnings and checkpoints
	! delta_eig	:	vector of differences between consecutives arrays [ Δλ_i = λ_(i+1) - λ_i ]
	! s_i		:	the normalized spacing between eigenvalues [ s_i = Δλ_i / AVG(Δλ_i) ]
	! spacings	:	list of possible spacings to be used to evaluate local AVG(Δλ_i)
	! s_i_n	:	For a generic n, average local spacing using n eigenvalues to evaluate
	! counter	:	integer variable to be used to count checkpoints steps

	USE DEBUG
	USE RNG
	USE SORT

	IMPLICIT NONE

	REAL*8, DIMENSION(:), ALLOCATABLE :: diag, delta_eig, s_i
	REAL*8 :: num1, num2
	INTEGER*8 :: counter, dim_, ii, jj
	LOGICAL :: Debug_flag
	INTEGER*8, DIMENSION(7) :: spacings
	! Auxiliary variables to evaluate spacings (Fortran complains about long lines)
	REAL*8 :: temp_res
	INTEGER*8 :: temp_spac
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: s_i_a
	CHARACTER(4) :: file_name
	
	COMMON counter
	counter = 0

	! Debugging
	Debug_flag = .FALSE.
	
	! Opening file to read dimensions
	OPEN(10, file="dim_file.txt", status='old', action='READ')
	
	! Reading dimensions from file
	READ(10,*) dim_
	
	! Closing file
	CLOSE(10)
	
	! Checking dimensions greater than 0
	IF (Debug_flag.AND.(dim_.LE.0)) THEN
		CALL Print_error(Debug_flag, "Dim<=0")
	END IF
	
	
	! Defining spacings (it is intended to be fixed)
	spacings = [dim_/200, dim_/100, dim_/50, dim_/25, dim_/10, dim_/5, dim_/2]

	! Writing "spacings" to file so that python can read it
	OPEN(11, file="spacings.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,SIZE(spacings)
		WRITE(11,"(I6.5)")spacings(ii)
	END DO
	
	CLOSE(11)

! ---------- ALLOCATING VARIABLES ----------

	! Allocating diag vector using file
	ALLOCATE(diag(dim_))
	
	! Allocating differeces vector using file
	ALLOCATE(delta_eig(dim_))
	
	! Allocating matrix for alternatives spacings
	ALLOCATE(s_i_a(SIZE(spacings), dim_-1))

! ---------- INITIALIZING MATRIX (with gaussian random numbers) ----------

	! Case dim_ is odd
	IF (MOD(dim_,2).EQ.1) THEN
		DO ii=1,((dim_-1)/2)
			CALL Normal_rng(num1,num2)
			diag(2*ii-1) = num1
			diag(2*ii) = num2
		END DO
		CALL Normal_rng(num1,num2)
		diag(dim_) = num1
	ELSE
	! Case dim_ is even
		DO ii=1,(dim_/2)
			CALL Normal_rng(num1,num2)
			diag(2*ii-1) = num1
			diag(2*ii) = num2
		END DO
	END IF


! ---------- SORTING ----------

	CALL BubbleSort(diag)



! ---------- COMPUTING DIFFERENCES ----------

	! Evaluating differences between consecutive eigenvalues
	delta_eig(:) = diag(2:) - diag(:dim_-1)
	
	s_i = delta_eig/(SUM(delta_eig)/(dim_-1))


! ---------- EVALUATING RESULTS FOR DIFFERENT SPACINGS ----------

	! Evaluating different spacings
	DO ii=1,SIZE(spacings)
	
		temp_spac = spacings(ii)
		
		! Looping over all possible differences [ Δλ_i = λ_(i+1) - λ_i ]
		DO jj=1, (dim_-1)
			
			CALL Help_spacings(spacings, ii, diag, jj, dim_, temp_res)
			s_i_a(ii,jj) = temp_res

		END DO
	END DO

! ---------- PRINTING TO FILE ----------

	! Opening and REWRITING existing file
	OPEN(11, file="s_i_distr_real.txt", status='REPLACE', action='WRITE')
	
	DO ii=1,(dim_-1)
		WRITE(11,"(F13.7)") s_i(ii)
	END DO
	
	CLOSE(11)
	
	! Writing "local spacing" results
	DO ii=1,SIZE(spacings)
		! Converting spacings(ii) to string
		WRITE(file_name,"(I4.4)") spacings(ii)
		
		OPEN(11, file="s_i_distr_real_local"//TRIM("_"//file_name)//".txt", status='REPLACE', action='WRITE')
		
		DO jj=1,dim_-1
			WRITE(11,"(F13.7)")s_i_a(ii,jj)
		END DO
		
		CLOSE(11)
	END DO
	
	CLOSE(11)
	
! ---------- DEALLOCATING ----------
	
	DEALLOCATE(diag)
	DEALLOCATE(s_i)
	DEALLOCATE(s_i_a)
	DEALLOCATE(delta_eig)
	
	END PROGRAM Ex5_real

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
