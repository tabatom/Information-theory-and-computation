	MODULE MY_MATH
	
	CONTAINS

	FUNCTION Tensor_product(matrix_1, matrix_2)
	
		IMPLICIT NONE
		
		COMPLEX*16, DIMENSION(:,:) :: matrix_1, matrix_2
		COMPLEX*16, DIMENSION(SIZE(matrix_1,1)*SIZE(matrix_2,1),SIZE(matrix_1,2)*SIZE(matrix_2,2)) :: Tensor_product
		INTEGER*8 :: ii, jj, l_r, l_c, r_r, r_c
		
		l_r = SIZE(matrix_1,1)
		l_c = SIZE(matrix_1,2)
		
		r_r = SIZE(matrix_2,1)
		r_c = SIZE(matrix_2,2)
		
		DO ii=0, l_r-1
			DO jj=0, l_c-1
				Tensor_product(ii*r_r+1:(ii+1)*r_r,jj*r_c+1:(jj+1)*r_c) = matrix_1(ii+1,jj+1)*matrix_2(:,:)
			END DO
		END DO
		
		RETURN
		
	END FUNCTION Tensor_product



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
