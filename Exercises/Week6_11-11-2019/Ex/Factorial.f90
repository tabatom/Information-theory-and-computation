	MODULE MY_MATH
	
	CONTAINS
	
	FUNCTION Factorial(n)
		
		IMPLICIT NONE
		
		INTEGER*8, INTENT(IN) :: n
		INTEGER*8 :: f_count
		INTEGER*8 :: Factorial
		
		IF (n.LT.0) THEN
			WRITE(*,*)"Factorial is not defined for n<0."
			Factorial = 0
			RETURN

		ELSEIF ((n.EQ.0).OR.(n.EQ.1)) THEN
			Factorial = 1
			RETURN
		ELSE
			Factorial = 1
			DO f_count=n,1,-1
				Factorial = Factorial*f_count
			END DO
			RETURN
		END IF
	END FUNCTION Factorial
	
	END MODULE MY_MATH
