
SUBROUTINE Checkpoint(Debug, Error)
	! Subroutine used to DEBUG
	
	IMPLICIT NONE
	
	LOGICAL :: Debug
	CHARACTER(LEN=*) :: Error
	
	IF (Debug) THEN
		SELECT CASE (Error)

		CASE("Dim=0")
			WRITE(*,*) "ERROR: dimensions are set to 0!"
			STOP

		CASE("dim<0")
			WRITE(*,*) "ERROR: dimensions are negative!"
			STOP

		CASE("Another error")
			WRITE(*,*) "ERROR: generic error."
			STOP

		CASE DEFAULT
			! Nothing to do

		END SELECT
	END IF
	
	RETURN
	
END SUBROUTINE Checkpoint
