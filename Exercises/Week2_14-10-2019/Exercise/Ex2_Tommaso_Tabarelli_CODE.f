! Program	: Ex2 (Fortran90)
!			The program should prepare a structure to deal with
!			matrices and their main characteristics
! Author	: Tommaso Tabarelli
! Date	: 15/10/2019


	! Module and type to be used
	MODULE MATRICES
		IMPLICIT NONE

		TYPE DMATRIX
			! Elem	: matrix elements
			! Dims	: matrix dimensions (Elem effective dimension and Dims should be equal)
			! Trace	: matrix trace
			! Det		: matrix determinant
			COMPLEX*8, DIMENSION(:,:), ALLOCATABLE :: Elem
			INTEGER, DIMENSION(2) :: Dims
			COMPLEX*8 :: Trace, Det
		END TYPE DMATRIX

		! Defining interface to use to initialize DMATRIX objects
		INTERFACE init
			MODULE PROCEDURE Init_matr
		END INTERFACE

		! Defining operator
		INTERFACE OPERATOR(.Trc.)
			MODULE PROCEDURE Mat_trc,Trace_num
		END INTERFACE

		INTERFACE OPERATOR(.Adj.)
			MODULE PROCEDURE Mat_adj
		END INTERFACE

		CONTAINS
		
		! TRACE FUNCTION returning type: DMATRIX
		FUNCTION Mat_trc(Matrix_)
			! Matrix_	: TYPE DMATRIX

			IMPLICIT NONE

			TYPE(DMATRIX), INTENT(IN) :: Matrix_
			TYPE(DMATRIX) :: Mat_trc
			COMPLEX*8 :: Trc
			INTEGER :: ii

			! Verifiyng the matrix is square and returning a "NULL" result if not
			IF (Matrix_%Dims(1).NE.(Matrix_%Dims(2))) THEN
				WRITE(*,*) "The matrix is not square. Returning a NULL result."
				Mat_trc%Dims = (0,0)
				ALLOCATE( Mat_trc%Elem(0,0) )
				Mat_trc%Elem = (0d0,0d0)
				Mat_trc%Trace = (0d0,0d0)
				Mat_trc%Det = (0d0,0d0)
				RETURN
			END IF
			
			Mat_trc%Dims(1) = Matrix_%Dims(1)
			Mat_trc%Dims(2) = Matrix_%Dims(2)
			
			ALLOCATE( Mat_trc%Elem(Mat_trc%Dims(1),Mat_trc%Dims(2)) )


			! Initialize matrix (or suppose it is already initialized)
			Mat_trc%Elem = Matrix_%Elem
			
			Trc = 0d0
			
			Mat_trc%Trace = 0d0

			! Evaluating trace (supposing the matrix is square)
			DO ii=1,Mat_trc%Dims(1)
				Trc = Trc + Matrix_%Elem(ii,ii)
			END DO

			Mat_trc%Trace = Trc

			RETURN
		END FUNCTION Mat_trc


		! TRACE FUNCTION returning type: COMPLEX*8
		FUNCTION Trace_num(Matrix_elem)
			! Matrix_elem	: TYPE COMPLEX*8 bidimensional array

			IMPLICIT NONE

			!TYPE(DMATRIX), INTENT(IN) :: Matrix_
			COMPLEX*8, DIMENSION(:,:), INTENT(IN) :: Matrix_elem
			COMPLEX*8 :: Trace_num
			INTEGER :: ii
			INTEGER*8, DIMENSION(2) :: dims

			Trace_num = 0d0
			
			! Getting proper dimensions
			dims = SHAPE(Matrix_elem)

			! Verifiyng the matrix is square and returning a "NULL" result if not
			IF (dims(1).NE.(dims(2))) THEN
				WRITE(*,*) "The matrix is not square. Returning a NULL result."
				Trace_num = (0d0,0d0)
				RETURN
			END IF

			! Evaluating trace (supposing the matrix is square, so
			!	looping only in 1 dimension)
			DO ii=1,dims(1)
				Trace_num = Trace_num + Matrix_elem(ii,ii)
			END DO

			RETURN
		END FUNCTION Trace_num


		! ADJOINT FUNCTION
		FUNCTION Mat_adj(Matrix_)
			! Matrix_ is of TYPE DMATRIX	: it is the input object

			IMPLICIT NONE

			INTEGER :: ii
			TYPE(DMATRIX), INTENT(IN) :: Matrix_
			TYPE(DMATRIX) :: Mat_adj

			Mat_adj%Dims(1)=Matrix_%Dims(2)
			Mat_adj%Dims(2)=Matrix_%Dims(1)
			
			Mat_adj%Trace = CONJG(Matrix_%Trace)

			Mat_adj%Det = CONJG(Matrix_%Det)

			Mat_adj%Elem = CONJG(TRANSPOSE(Matrix_%Elem))

			RETURN
		END FUNCTION Mat_adj

		! INITIALIZING FUNCTION
		FUNCTION Init_matr(Dims, Value_)
			! Dims	: matrix dimensions
			! Value_	: value to initialize ALL matrix elements
			
			INTEGER, DIMENSION(2) :: Dims
			COMPLEX*8 :: Value_
			TYPE(DMATRIX) Init_matr
			INTEGER :: ii,jj

			Init_matr%Dims(1) = Dims(1)
			Init_matr%Dims(2) = Dims(2)
			
			! Here it should be better to use a WHILE (I cant use it for now)
			! 	to eventually ask again for dimensions, or eventually
			!	return a "NULL" result to be collected by a loop the function is in
			IF ( (Dims(1).LE.0).OR.(Dims(2).LE.0) ) THEN
				WRITE(*,*) "ERROR: matrix dimensions CAN NOT BE negative or 0"
				STOP
			END IF
			
			ALLOCATE( Init_matr%Elem(Dims(1),Dims(2)) )
			
			DO ii=0,Dims(1)
				DO jj=0,Dims(2)
					Init_matr%Elem(ii,jj)=Value_
				END DO
			END DO
			
			! Initializing trace and determinant to 0
			Init_matr%Trace = (0d0, 0d0)
			Init_matr%Det = (0d0, 0d0)
			
			RETURN
			
		END FUNCTION Init_matr

		! Subroutine to print matrix
		SUBROUTINE Dump_Matr(Matrix_, File_name)
			! Matrix_	: TYPE DMATRIX

			IMPLICIT NONE

			TYPE(DMATRIX), INTENT(IN) :: Matrix_
			CHARACTER(LEN=*), INTENT(IN) :: File_name
			INTEGER*8 :: ii

			! Opening file (ID=51)
			OPEN(UNIT=51,FILE=File_name,STATUS="unknown")
			
			WRITE(51,*) "Matrix dimensions are: "
			WRITE(51,*) Matrix_%Dims
			WRITE(51,*) "Matrix trace is: "
			WRITE(51,*) Matrix_%Trace
			WRITE(51,*) "Matrix determinant is: "
			WRITE(51,*) Matrix_%Det
			
			WRITE(51,*) "Matrix elements are: "
			
			DO ii=1,Matrix_%Dims(1)
				WRITE(51,*) Matrix_%Elem(ii,:)
			END DO
			
			CLOSE(51)
			
		END SUBROUTINE Dump_Matr

	END MODULE MATRICES



	! Starting program
	PROGRAM Ex2_MATRICES
		USE MATRICES
		IMPLICIT NONE
		
		TYPE(DMATRIX) My_try, Trace_
		INTEGER, DIMENSION(2) :: Dims
		COMPLEX*8 :: comp_Val
		INTEGER*8, DIMENSION(2) :: dim1
		
		
		WRITE(*,*) "Insert the 2 dimensions"
		READ(*,*) Dims(1), Dims(2)
		
		WRITE(*,*) "Insert the value to initialize the matrix"
		READ(*,*) comp_Val
		
		My_try = init(Dims, comp_Val)
		
		! Want to overwrite the DMATRIX onto itself
		!My_try = .Trc.My_try
		
		Trace_ = .Trc.My_try
		My_try%Trace = .Trc.My_try%Elem

		WRITE(*,*) My_try%Elem
		WRITE(*,*) Trace_%Elem
		WRITE(*,*) My_try%Trace
		WRITE(*,*) Trace_%Trace
		
		CALL Dump_Matr(My_try, "Matrix.txt")
		
		CALL Dump_Matr(.Adj.My_try, "Adjoint_Matrix.txt")
		
		STOP
	END PROGRAM Ex2_MATRICES
