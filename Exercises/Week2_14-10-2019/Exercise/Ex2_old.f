! Program	: Ex2 (Fortran90)
!			The program should prepare a structure to deal with
!			matrices and their main characteristics
! Author	: Tommaso Tabarelli
! Date	: 15/10/2019


! Module and type to be used
	MODULE MATRICES
		IMPLICIT NONE

		TYPE DMATRIX
			COMPLEX*8, DIMENSION(:,:), ALLOCATABLE :: Elem
			INTEGER, DIMENSION(2) :: Dims
			COMPLEX*8 :: Trace, Det
		END TYPE DMATRIX

	! Defining operator
		INTERFACE OPERATOR(.Trc.)
			MODULE PROCEDURE Mat_trc
		END INTERFACE

		INTERFACE OPERATOR(.Adj.)
			MODULE PROCEDURE Mat_adj
		END INTERFACE

! TRACE FUNCTION
		CONTAINS
		FUNCTION Mat_trc(M_Values, M_Dim)
			! M_Values	: matrix elements
			! M_Dim	: matrix dimensions

			IMPLICIT NONE

			COMPLEX*8 :: Mat_trc
			INTEGER, DIMENSION(2), INTENT(IN) :: M_Dim
			COMPLEX*8, DIMENSION(M_Dim(1),M_Dim(2)), INTENT(IN) :: M_Values
			INTEGER :: N, ii

		! N.B.: SUPPOSING THE MATRIX IS SQUARE (so collecting only 1 dimension)
			N = M_Dim(1)
			Mat_trc = 0d0

			! Initialize matrix (or suppose it is already initialized)

			! Evaluating trace
			DO ii=1,N
				Mat_trc = Mat_trc + M_Values(ii,ii)
			END DO

			RETURN
		END FUNCTION Mat_trc


! ADJOINT FUNCTION
		FUNCTION Mat_adj(Matrix_)
			! Matrix_ is of TYPE DMATRIX	: it is the input object

			IMPLICIT NONE

		! N.B.: SUPPOSING THE MATRIX IS SQUARE (so collecting only 1 dimension)
			INTEGER :: ii
			TYPE(DMATRIX), INTENT(IN) :: Matrix_
			TYPE(DMATRIX) :: Mat_adj

			! Initialize matrix (or suppose it is already initialized)

			Mat_adj%Dims=Matrix_%Dims

			Mat_adj%Det = Matrix_%Det

			Mat_adj%Elem = CONJG(TRANSPOSE(Matrix_%Elem))

			RETURN
		END FUNCTION Mat_adj

! INITIALIZING FUNCTION
		FUNCTION Mat_init(Dims, Value_)
			! Dims	: matrix dimensions
			! Value_	: value to initialize ALL matrix elements
			
			INTEGER, DIMENSION(2) :: Dims
			COMPLEX*8 :: Value_
			TYPE(DMATRIX) Mat_init
			INTEGER :: ii,jj
			
			! Here it should be better to use a WHILE (I cant use it for now)
			IF ( (Dims(1).LE.0).OR.(Dims(2).LE.0) ) THEN
				WRITE(*,*) "ERROR: matrix dimensions CAN NOT BE negative or 0"
				STOP
			END IF
			
			ALLOCATE( Mat_init%Elem(Dims(1),Dims(2)) )
			
			DO ii=0,Dims(1)
				DO jj=0,Dims(1)
					Mat_init%Elem(ii,jj)=Value_
				END DO
			END DO
			
			RETURN
			
		END FUNCTION Mat_init

	END MODULE MATRICES


! Starting program
	PROGRAM Ex2_MATRICES
		
		STOP
	END PROGRAM Ex2_MATRICES
