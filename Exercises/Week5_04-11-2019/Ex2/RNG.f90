	MODULE RNG
		! Module to generate Normal random number
		IMPLICIT NONE
		
		INTERFACE Normal_rng
			MODULE PROCEDURE Normal_rng_comp8
		END INTERFACE Normal_rng
		
		CONTAINS
			SUBROUTINE Normal_rng_comp8(x, y)
				REAL*8 :: u1, u2, r, theta
				REAL*8, INTENT(OUT) :: x, y
				
				CALL RANDOM_NUMBER(u1)
				CALL RANDOM_NUMBER(u2)
				
				r = SQRT(2*(-LOG(u1)))
				! π = 2*arcsin(1)
				theta = 2*(2*ASIN(1.)) * u2
				
				x = r*COS(theta)
				y = r*SIN(theta)

			END SUBROUTINE Normal_rng_comp8

	END MODULE RNG
