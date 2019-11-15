	MODULE DEBUG
	
	INTERFACE Checkpoint
		MODULE PROCEDURE CP_matrix_real8, CP_matrix_comp8, CP_matrix_comp16, CP_vec_real8
	END INTERFACE
	
	CONTAINS
	
	SUBROUTINE Print_error(Debug, Error)
		! Subroutine used to DEBUG
		! Arguments:
		!	Debug		: logical, flag which confirms the call for debugging
		!	Error		: string, to be used to identify the error and print
		!					proper message
		
		IMPLICIT NONE
		
		LOGICAL :: Debug
		CHARACTER(LEN=*) :: Error
		
		IF (Debug) THEN

			SELECT CASE (Error)

			CASE("Dim<=0")
				WRITE(*,*) "ERROR: dimensions are negative or 0!"
				STOP

			CASE("No match dim")
				WRITE(*,*) "ERROR: matrices dimensions mismatch: can not proceed with matrix product."
				STOP

			CASE("Timem<0")
				WRITE(*,*) "WARNING: time for manual operation is can be meaningless!"
			
			CASE("TimeT<0")
				WRITE(*,*) "WARNING: time for manual operation using tranposed left matrix can be meaningless!"
     			
     			CASE("TimeF<0")
				WRITE(*,*) "WARNING: time for automatic operation using matmul can be meaningless!"
     			
			CASE DEFAULT
				! Nothing to do
				WRITE(*,*) "WRONG ERROR IDENTIFIER: PLEASE, CHECK THE ERROR STRING"
				STOP

			END SELECT
		END IF
		
		RETURN
		
	END SUBROUTINE Print_error
	
	
	! Writing CHECKPOINT subroutine for REAL MATRICES
	SUBROUTINE CP_matrix_real8(Debug, IO, label, real_matrix, dims)
		! Subroutine used to save/retrieve stuff from file
		! Arguments are:
		!	Debug		: logical, flag which confirms the call for debugging
		!	IO		: string, should be "Read" or "Write", tells if read or write on log file
		!	label		: string. File names are: ".Status_" + label + ".log"
		!	real_matrix	: real*8 2D array (in this specific case); is a variable to be stored
		!	dims		: integer*8 array with 2 elements. Used to eventually print a 0s matrix in the file.
		!				dims should correspond to actual real_matrix dimensions.
		
		
		IMPLICIT NONE
			
		LOGICAL, INTENT(IN) :: Debug
		CHARACTER(LEN=*), INTENT(IN) :: IO, label
		CHARACTER(50) :: file_name
		INTEGER*8, DIMENSION(2), INTENT(IN) :: dims
		REAL*8, DIMENSION(dims(1),dims(2)), INTENT(INOUT) :: real_matrix
		INTEGER*8 :: stat, ii, jj, counter
		COMMON counter
		
		file_name = ".Status_" // TRIM(label) // ".log"
		file_name = TRIM(file_name)
	
	
		IF (Debug) THEN	
			! Retrieving data
			IF (IO.EQ."Read") THEN
			
				! Opening file to read previous state
				!  (iostat is a collecting error variable)
				OPEN(10, file=file_name, status='old', action='READ', iostat=stat)
				
				IF (stat.EQ.0) THEN
					READ(10,*) real_matrix
				ELSE
					! If no data to retrieve, then start from scratch
					DO ii=1,dims(1)
						DO jj=1,dims(2)
							real_matrix(ii,jj)=0d0
						END DO
					END DO
				END IF
				
				! Closing file
				CLOSE(10)

			END IF


			! Writing/saving data
			IF (IO.EQ."Write") THEN
				
				! Opening file to write current state
				OPEN(10, file=file_name, status="REPLACE", action="WRITE", iostat=stat)
				
				! Writing data to file
				WRITE(10,*) real_matrix
				
				! Closing file
				CLOSE(10)


				! Opening another file in which save whole history
				
				! Creating file if not existing
				IF (counter.EQ.0) THEN
					OPEN(10, file=".debug_history.log", status="REPLACE")
					CLOSE(10)
				END IF
				
				OPEN(10, file=".debug_history.log", status="old", action="WRITE", iostat=stat, position="append")
				WRITE(10,*) ""
				WRITE(10,*) "Call n: ", counter, label
				counter = counter + 1
				WRITE(10,*) ""
				WRITE(10,*) real_matrix
				
				CLOSE(10)
				
				
			END IF
		END IF
	
	END SUBROUTINE CP_matrix_real8
	
	! Writing CHECKPOINT subroutine for COMPLEX*8 MATRICES
	SUBROUTINE CP_matrix_comp8(Debug, IO, label, complex_matrix, dims)
		! Subroutine used to save/retrieve stuff from file
		! Arguments are:
		!	Debug			: logical, flag which confirms the call for debugging
		!	IO			: string, should be "Read" or "Write", tells if read or write on log file
		!	label			: string. File names are: ".Status_" + label + ".log"
		!	complex_matrix	: complex*8 2D array (in this specific case); is a variable to be stored
		!	dims			: integer*8 array with 2 elements. Used to eventually print a 0s matrix in the file.
		!					dims should correspond to actual real_matrix dimensions.
		
		
		IMPLICIT NONE
			
		LOGICAL, INTENT(IN) :: Debug
		CHARACTER(LEN=*), INTENT(IN) :: IO, label
		CHARACTER(50) :: file_name
		INTEGER*8, DIMENSION(2), INTENT(IN) :: dims
		COMPLEX*8, DIMENSION(dims(1),dims(2)), INTENT(INOUT) :: complex_matrix
		INTEGER*8 :: stat, ii, jj, counter
		COMMON counter
		
		file_name = ".Status_" // TRIM(label) // ".log"
		file_name = TRIM(file_name)
	
	
		IF (Debug) THEN	
			! Retrieving data
			IF (IO.EQ."Read") THEN
			
				! Opening file to read previous state
				!  (iostat is a collecting error variable)
				OPEN(10, file=file_name, status='old', action='READ', iostat=stat)
				
				IF (stat.EQ.0) THEN
					READ(10,*) complex_matrix
				ELSE
					! If no data to retrieve, then start from scratch
					DO ii=1,dims(1)
						DO jj=1,dims(2)
							complex_matrix(ii,jj)=(0d0,0d0)
						END DO
					END DO
				END IF
				
				! Closing file
				CLOSE(10)

			END IF


			! Writing/saving data
			IF (IO.EQ."Write") THEN
				
				! Opening file to write current state
				OPEN(10, file=file_name, status="REPLACE", action="WRITE", iostat=stat)
				
				! Writing data to file
				WRITE(10,*) complex_matrix
				
				! Closing file
				CLOSE(10)


				! Opening another file in which save whole history
				
				! Creating file if not existing
				IF (counter.EQ.0) THEN
					OPEN(10, file=".debug_history.log", status="REPLACE")
					CLOSE(10)
				END IF
				
				OPEN(10, file=".debug_history.log", status="old", action="WRITE", iostat=stat, position="append")
				WRITE(10,*) ""
				WRITE(10,*) "Call n: ", counter, label
				counter = counter + 1
				WRITE(10,*) ""
				
				! Writing matrix row-wise
				DO ii=1,dims(1)
					WRITE(10,*) complex_matrix(ii,:)
				END DO
				
				CLOSE(10)
				
				
			END IF
		END IF

	END SUBROUTINE CP_matrix_comp8
	
	! Writing CHECKPOINT subroutine for COMPLEX*8 MATRICES
	SUBROUTINE CP_matrix_comp16(Debug, IO, label, complex_matrix, dims)
		! Subroutine used to save/retrieve stuff from file
		! Arguments are:
		!	Debug			: logical, flag which confirms the call for debugging
		!	IO			: string, should be "Read" or "Write", tells if read or write on log file
		!	label			: string. File names are: ".Status_" + label + ".log"
		!	complex_matrix	: complex*16 2D array (in this specific case); is a variable to be stored
		!	dims			: integer*8 array with 2 elements. Used to eventually print a 0s matrix in the file.
		!					dims should correspond to actual real_matrix dimensions.
		
		
		IMPLICIT NONE
			
		LOGICAL, INTENT(IN) :: Debug
		CHARACTER(LEN=*), INTENT(IN) :: IO, label
		CHARACTER(50) :: file_name
		INTEGER*8, DIMENSION(2), INTENT(IN) :: dims
		COMPLEX*16, DIMENSION(dims(1),dims(2)), INTENT(INOUT) :: complex_matrix
		INTEGER*8 :: stat, ii, jj, counter
		COMMON counter
		
		file_name = ".Status_" // TRIM(label) // ".log"
		file_name = TRIM(file_name)
	
	
		IF (Debug) THEN	
			! Retrieving data
			IF (IO.EQ."Read") THEN
			
				! Opening file to read previous state
				!  (iostat is a collecting error variable)
				OPEN(10, file=file_name, status='old', action='READ', iostat=stat)
				
				IF (stat.EQ.0) THEN
					READ(10,*) complex_matrix
				ELSE
					! If no data to retrieve, then start from scratch
					DO ii=1,dims(1)
						DO jj=1,dims(2)
							complex_matrix(ii,jj)=(0d0,0d0)
						END DO
					END DO
				END IF
				
				! Closing file
				CLOSE(10)

			END IF


			! Writing/saving data
			IF (IO.EQ."Write") THEN
				
				! Opening file to write current state
				OPEN(10, file=file_name, status="REPLACE", action="WRITE", iostat=stat)
				
				! Writing data to file
				WRITE(10,*) complex_matrix
				
				! Closing file
				CLOSE(10)


				! Opening another file in which save whole history
				
				! Creating file if not existing
				IF (counter.EQ.0) THEN
					OPEN(10, file=".debug_history.log", status="REPLACE")
					CLOSE(10)
				END IF
				
				OPEN(10, file=".debug_history.log", status="old", action="WRITE", iostat=stat, position="append")
				WRITE(10,*) ""
				WRITE(10,*) "Call n: ", counter, label
				counter = counter + 1
				WRITE(10,*) ""
				
				! Writing matrix row-wise
				DO ii=1,dims(1)
					WRITE(10,*) complex_matrix(ii,:)
				END DO
				
				CLOSE(10)
				
				
			END IF
		END IF

	END SUBROUTINE CP_matrix_comp16
	
	! Writing CHECKPOINT subroutine for REAL*8 VECTORS
	SUBROUTINE CP_vec_real8(Debug, IO, label, real_vec, dim_)
		! Subroutine used to save/retrieve stuff from file
		! Arguments are:
		!	Debug		: logical, flag which confirms the call for debugging
		!	IO		: string, should be "Read" or "Write", tells if read or write on log file
		!	label		: string. File names are: ".Status_" + label + ".log"
		!	real_vec	: real*8 1D array (in this specific case); is a variable to be stored
		!	dim_		: integer*8. Used to eventually print a 0s vector in the file.
		!				dims should correspond to actual real_vec dimension.
		
		
		IMPLICIT NONE
			
		LOGICAL, INTENT(IN) :: Debug
		CHARACTER(LEN=*), INTENT(IN) :: IO, label
		CHARACTER(50) :: file_name
		INTEGER*8, INTENT(IN) :: dim_
		REAL*8, DIMENSION(dim_), INTENT(INOUT) :: real_vec
		INTEGER*8 :: stat, ii, counter
		COMMON counter
		
		file_name = ".Status_" // TRIM(label) // ".log"
		file_name = TRIM(file_name)
	
	
		IF (Debug) THEN	
			! Retrieving data
			IF (IO.EQ."Read") THEN
			
				! Opening file to read previous state
				!  (iostat is a collecting error variable)
				OPEN(11, file=file_name, status='old', action='READ', iostat=stat)
				
				IF (stat.EQ.0) THEN
					READ(11,*) real_vec
				ELSE
					! If no data to retrieve, then start from scratch
					DO ii=1,dim_
						real_vec(ii)=0d0
					END DO
				END IF
				
				! Closing file
				CLOSE(11)

			END IF


			! Writing/saving data
			IF (IO.EQ."Write") THEN
				
				! Opening file to write current state
				OPEN(11, file=file_name, status="REPLACE", action="WRITE", iostat=stat)
				
				! Writing data to file
				WRITE(11,*) real_vec
				
				! Closing file
				CLOSE(11)


				! Opening another file in which save whole history
				
				! Creating file if not existing
				IF (counter.EQ.0) THEN
					OPEN(11, file=".debug_history.log", status="REPLACE")
					CLOSE(11)
				END IF
				
				OPEN(11, file=".debug_history.log", status="old", action="WRITE", iostat=stat, position="append")
				WRITE(11,*) ""
				WRITE(11,*) "Call n: ", counter, label
				counter = counter + 1
				WRITE(11,*) ""
				WRITE(11,*) real_vec
				
				CLOSE(11)
				
				
			END IF
		END IF
	
	END SUBROUTINE CP_vec_real8
	
	END MODULE DEBUG

