	PROGRAM Ex3_debug
	! Program that evaluate manually the matrix multiplication
	! Variables are:
	! matrix_L	:	leftmost matrix in matrix moltiplication
	! matrix_R	:	rightmost matrix "    "       "
	! result_m	:	the result matrix using standard intuitive algorithm
	! result_F	:	result matrix using "MATMUL" function
	! transp_L	:	transpose of left matrix, to try to use a different algorithm (fortran allocates matrixes with columns elements near)
	! result_T	:	result matrix using transp_L instead of matrix_L
	! ii,jj,kk	:	indeces/numbers to be used when needed
	! Debug	:	logical value to enable runtime error warnings
	! test_dim	:	logical value to make test (fortran compels about tests with more conditions (I do not know why) )

	IMPLICIT NONE

      REAL*8, DIMENSION(:,:),ALLOCATABLE :: matrix_L,matrix_R,result_m
      REAL*8, DIMENSION(:,:),ALLOCATABLE :: result_F,transp_L,result_T
      INTEGER*2 :: ii,jj,kk
      INTEGER*8, DIMENSION(2) :: dim_L,dim_R,dim_res,counter
      LOGICAL :: Debug, test_dim

! Variables to measure the time spent
      REAL*8 :: start_m,finish_m,start_F,finish_F,start_T,finish_T,app

      COMMON counter
      counter = 0

	! Debugging
	Debug = .TRUE.

! Asking the user to put dimensions of rows and columns
! N.B.: here the dimensions should be acceptable to do a matrix multiplication (i.e. col_L = row_R)
      !PRINT*,"Insert row dimension of left matrix"
	!READ*,row_L
      !PRINT*,"Insert column dimension of left matrix"
	!READ*,col_L
      !PRINT*,"Insert row dimension of right matrix"
	!READ*,row_R
      !PRINT*,"Insert column dimension of right matrix"
	!READ*,col_R

	! To avoid input requests, setting dimensions to 500 (arbitrary)
	dim_L(1)=5
	dim_L(2)=5
	dim_R(1)=5
	dim_R(2)=5
	
	dim_res(1)=dim_L(1)
	dim_res(2)=dim_R(2)

	! Checking dimensions greater than 0 (doing it on more lines
	!	because the condition is too long for fortran...)
	IF (Debug) THEN
		test_dim = ((dim_L(1).LE.0).OR.(dim_L(2).LE.0))
		test_dim = test_dim.OR.((dim_R(1).LE.0).OR.(dim_R(2).LE.0))
		IF (test_dim) THEN
			CALL Print_error(Debug, "Dim<=0")
		END IF
	END IF

	! Checking corresponding dimensions for matrices multiplication
	IF ((Debug) .AND. (dim_L(2).NE.dim_R(1))) THEN
		CALL Print_error(Debug, "No match dim")
	END IF

      ! Allocating matrix memories
      ALLOCATE(matrix_L(dim_L(1),dim_L(2)),matrix_R(dim_R(1),dim_R(2)))
      ALLOCATE(result_m(dim_L(1),dim_R(2)),result_F(dim_L(1),dim_R(2)))
	ALLOCATE(transp_L(dim_L(2),dim_L(1)),result_T(dim_L(1),dim_R(2)))
      

! ---------- INITIALIZING MATRICES ----------

! matrix_L
	
	! Starting checkpointing (and eventually initializing to 0 all elements)
	CALL Checkpoint(Debug, 'Read', "result_L", matrix_L, dim_L)

	! Initializing with incremental values (arbitrary choice)
	kk = 0
	DO ii=1,dim_L(1)
		DO jj=1,dim_L(2)
			CALL Checkpoint(Debug, 'Read', 'matrix_L', matrix_L, dim_L)
			matrix_L(ii,jj)=kk
			CALL Checkpoint(Debug, 'Write', 'matrix_L', matrix_L, dim_L)
			kk = kk + 1
		END DO
	END DO
	kk = 0

! matrix_R

	! Starting checkpointing (and eventually initializing to 0 all elements)
	CALL Checkpoint(Debug, 'Read', "result_R", matrix_R, dim_R)

	kk = 0
	DO ii=1,dim_R(1)
		DO jj=1,dim_R(2)
			CALL Checkpoint(Debug, 'Read', 'matrix_R', matrix_R, dim_R)
			matrix_R(ii,jj)=kk
	         	CALL Checkpoint(Debug,'Write','matrix_R',matrix_R,dim_R)
			kk = kk + 1
	   END DO
	END DO
	kk = 0


! result_m

	! Starting checkpointing (and eventually initializing to 0 all elements)
	CALL Checkpoint(Debug, 'Read', "result_m", result_m, dim_res)
	CALL Checkpoint(Debug, 'Write', "result_m", result_m, dim_res)


! ---------- STARTING OPERATIONS ----------

! MANUAL OPERATION: measuring time
      CALL CPU_TIME(start_m)
	! This operations are considered starting from the "column" of right matrix (indented operations)

      DO ii=1,dim_R(2)
         DO jj=1,dim_L(1)
		! Resetting temporary result to 0
		app=0
            DO kk=1,dim_L(2)       ! col_L and row_R should be the same
		   app=app+matrix_L(jj,kk)*matrix_R(kk,ii)
            END DO
            
      	CALL Checkpoint(Debug, 'Read', "result_m", result_m, dim_res)
		result_m(jj,ii)=app
		CALL Checkpoint(Debug, 'Write', "result_m", result_m, dim_res)
         END DO
      END DO
      CALL CPU_TIME(finish_m)


! FUNCTION OPERATION: measuring time
	CALL Checkpoint(Debug, 'Read', "result_F", result_F, dim_res)
      CALL CPU_TIME(start_F)
      result_F=MATMUL(matrix_L,matrix_R)
      CALL CPU_TIME(finish_F)
	CALL Checkpoint(Debug, 'Write', "result_F", result_F, dim_res)
      
	
! TRANPOSING, THEN OPERATING
	! Trasposing left matrix and repeating calculous
	CALL Checkpoint(Debug, 'Read', "result_T", result_T, dim_res)

	CALL CPU_TIME(start_T)

	! Transposing left matrix
	DO ii=1,dim_L(2)
		DO jj=1,dim_L(1)
			! Reading from colums of matrix_L and writing to rows of transp_L
			transp_L(ii,jj)=matrix_L(jj,ii)
		END DO
	END DO

	DO ii=1,dim_R(2)
		DO jj=1,dim_L(1)
			! Resetting temporary result to 0
			app=0
			DO kk=1,dim_L(2)       ! col_L and row_R should be the same
				app=app+transp_L(kk,jj)*matrix_R(kk,ii)
			END DO
			result_T(jj,ii)=app
		END DO
      END DO
	CALL CPU_TIME(finish_T)
	CALL Checkpoint(Debug, 'Write', "result_T", result_T, dim_res)

!	Printing results
!      PRINT*, "Left matrix:"
!      PRINT*, matrix_L
!      PRINT*, "Right matrix:"
!      PRINT*, matrix_R
!      PRINT*, "Manual result matrix:"
!      PRINT*, result_m
!      PRINT*, "Transpose result matrix:"
!      PRINT*, result_T
!      PRINT*, "Function result matrix:"
!      PRINT*, result_F

	! Check if times are meaningful
	IF ((Debug).AND.(finish_m.LE.start_m)) THEN
		CALL Print_error(Debug, "Timem<0")
	END IF
	IF ((Debug).AND.(finish_T.LE.start_T)) THEN
		CALL Print_error(Debug, "TimeT<0")
	END IF
	IF ((Debug).AND.(finish_F.LE.start_F)) THEN
		CALL Print_error(Debug, "TimeF<0")
	END IF

      PRINT*,"Manual process took:"
      PRINT*,finish_m-start_m

      PRINT*,"Transposed manual process took:"
      PRINT*,finish_T-start_T	

      PRINT*,"Function process took:"
      PRINT*,finish_F-start_F

      ! Here it reises ERRORS
      DEALLOCATE(matrix_L,matrix_R,result_m,result_F,transp_L,result_T)
	STOP
      END PROGRAM Ex3_debug


	INCLUDE "Debug.f"
