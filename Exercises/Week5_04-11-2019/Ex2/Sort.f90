	MODULE SORT

	INTERFACE BubbleSort
		MODULE PROCEDURE BS_real8
	END INTERFACE BubbleSort
	
	CONTAINS
	
	SUBROUTINE BS_real8(array)
		! Simple implementation of bubblesort for 1D array
		REAL*8, DIMENSION(:), INTENT(INOUT) :: array
		INTEGER*8 :: i1, i2
		REAL*8 :: temp_num
		
		DO i2=1,SIZE(array)
			DO i1=1,(SIZE(array)-1)
				IF (array(i1).GE.array(i1+1)) THEN
					temp_num = array(i1)
					array(i1) = array(i1+1)
					array(i1+1) = temp_num
				END IF
			END DO
		END DO
		RETURN
	END SUBROUTINE BS_real8
	
	END MODULE SORT
