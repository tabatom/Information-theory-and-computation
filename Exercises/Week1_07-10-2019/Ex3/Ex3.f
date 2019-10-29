      PROGRAM Ex3
!     Program that evaluate manually the matrix multiplication
!
!	Variables are:
!	matrix_L	:	leftmost matrix in matrix moltiplication
!	matrix_R	:	rightmost matrix "    "       "
!	result_m	:	the result matrix using standard intuitive algorithm
!	result_F	:	result matrix using "MATMUL" function
!	transp_L	:	transpose of left matrix, to try to use a different algorithm (fortran allocates matrixes with columns elements near)
!	result_T	:	result matrix using transp_L instead of matrix_L
!	ii,jj,kk	:	indeces/numbers to be used when needed

      REAL*8, DIMENSION(:,:),ALLOCATABLE :: matrix_L,matrix_R,result_m
      REAL*8, DIMENSION(:,:),ALLOCATABLE :: result_F,transp_L,result_T
      INTEGER*2 :: ii,jj,kk,row_L,col_L,row_R,col_R

! Variables to measure the time spent
      REAL*8 :: start_m,finish_m,start_F,finish_F,start_T,finish_T,app

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
	row_L=500
	col_L=500
	row_R=500
	col_R=500

      ! Allocating matrix memories
      ALLOCATE(matrix_L(row_L,col_L),matrix_R(row_R,col_R))
      ALLOCATE(result_m(row_L,col_R),result_F(row_L,col_R))
	ALLOCATE(transp_L(col_L,row_L),result_T(row_L,col_R))
      
! Initializing matrices with incremental values (arbitrary choice)
! matrix_L
      kk = 0
      DO ii=1,row_L
         DO jj=1,col_L
            matrix_L(ii,jj)=kk
            kk = kk + 1
         END DO
      END DO
      kk = 0

! matrix_R
      kk = 0
      DO ii=1,row_R
         DO jj=1,col_R
            matrix_R(ii,jj)=kk
            kk = kk + 1
         END DO
      END DO
      kk = 0

! results
      DO ii=1,row_L
         DO jj=1,col_R
            ! Setting the result matrix empty, so I can just sum terms during the operation
            result_m(ii,jj)=0
         END DO
      END DO

      
! Checking dimensions of matrices
      IF (col_L.NE.row_R) THEN
         PRINT*,"Error: matrices dimensions not compatible!"
	   STOP
      END IF

! MANUAL OPERATION: measuring time
      CALL CPU_TIME(start_m)
!     This operations are considered starting from the "column" of right matrix (indented operations)

      DO ii=1,col_R
         DO jj=1,row_L
		! Resetting temporary result to 0
		app=0
            DO kk=1,col_L       ! col_L and row_R should be the same
		   app=app+matrix_L(jj,kk)*matrix_R(kk,ii)
!               result_m(jj,ii)=result_m(jj,ii)+
!     +              matrix_L(jj,kk)*matrix_R(kk,ii)
            END DO
		result_m(jj,ii)=app
         END DO
      END DO
      CALL CPU_TIME(finish_m)

      ! Measuring time for Function operation
      CALL CPU_TIME(start_F)
      result_F=MATMUL(matrix_L,matrix_R)
      CALL CPU_TIME(finish_F)
      
	
	! Trasposing left matrix and repeating calculous
	CALL CPU_TIME(start_T)

	! Transposing left matrix
	DO ii=1,col_L
		DO jj=1,row_L
			! Reading from colums of matrix_L and writing to rows of transp_L
			transp_L(ii,jj)=matrix_L(jj,ii)
		END DO
	END DO

	DO ii=1,col_R
         DO jj=1,row_L
		! Resetting temporary result to 0
		app=0
            DO kk=1,col_L       ! col_L and row_R should be the same
		    app=app+transp_L(kk,jj)*matrix_R(kk,ii)
!               result_T(jj,ii)=result_T(jj,ii)+
!     +              transp_L(kk,jj)*matrix_R(kk,ii)
            END DO
		result_T(jj,ii)=app
         END DO
      END DO
	CALL CPU_TIME(finish_T)

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

      PRINT*,"Manual process took:"
      PRINT*,finish_m-start_m

      PRINT*,"Transposed manual process took:"
      PRINT*,finish_T-start_T	

      PRINT*,"Function process took:"
      PRINT*,finish_F-start_F

      ! Here it reises ERRORS
      DEALLOCATE(matrix_L,matrix_R,result_m,result_F,transp_L,result_T)
	STOP
      END PROGRAM Ex3
