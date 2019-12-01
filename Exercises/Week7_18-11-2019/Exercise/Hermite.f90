! REFERENCES:
! https://sukhbinder.wordpress.com/hermite-polynomials-fortran-module/
! http://mathworld.wolfram.com/HermitePolynomial.html
! https://en.wikipedia.org/wiki/Hermite_polynomials


	MODULE HERMITE
	
	CONTAINS

	RECURSIVE FUNCTION HermitePoly(n) RESULT(hp_coeff)
	
	INTEGER*8 :: n
	REAL*8, DIMENSION(n+1) ::  hp_coeff, coef1, coef2

	IF(n .EQ. 0) THEN
		hp_coeff(1)=1.0
		RETURN
	END IF

	IF(n .EQ. 1) THEN
		hp_coeff(1)=2.0
		hp_coeff(2)=0.0
	ELSE

		coef1(1:n+1)=0.0
		coef1(1:n)=2.0*HermitePoly(n-1)

		coef2(1:n+1)=0.0
		coef2(3:)=2.0*(n-1)*HermitePoly(n-2)

		hp_coeff = coef1 - coef2

	END IF
	
	RETURN

	END FUNCTION
	
	FUNCTION evalHerm_Poly(xi,n) RESULT(y)
	
	INTEGER*8 :: n, pow_i, ii, jj
	REAL*8 :: xi(:), y(size(xi)),h_coeff(n+1)

	k=size(xi)

	h_coeff=HermitePoly(n)

	y(1:k)=h_coeff(n+1)
	
	pow_i=1
	
	DO ii=n,1,-1
		DO jj=1,k
			! Adding to the point in j-th position the 
			!  contribution given by the pow_i-th power
			!  of x_i multiplied by the proper coefficient
			y(jj) = y(jj)+h_coeff(ii)*xi(jj)**pow_i
		END DO
		! For simplicity, updating the power here
		pow_i = pow_i + 1
	END DO
	
	END FUNCTION

	
!	SUBROUTINE Herm0(x, y)
!		REAL*8, INTENT(IN) :: x
!		REAL*8, INTENT(OUT) :: y
!		y = 1
!	END SUBROUTINE Herm0
!	
!	SUBROUTINE Herm1(x, y)
!		REAL*8, INTENT(IN) :: x
!		REAL*8, INTENT(OUT) :: y
!		y = x
!	END SUBROUTINE Herm1
!	
!		SUBROUTINE Herm2(x, y)
!		REAL*8, INTENT(IN) :: x
!		REAL*8, INTENT(OUT) :: y
!		y = x**2 - 1
!	END SUBROUTINE Herm2
!	
!	SUBROUTINE Herm3(x, y)
!		REAL*8, INTENT(IN) :: x
!		REAL*8, INTENT(OUT) :: y
!		y = x**3 - 3*x
!	END SUBROUTINE
!	
!	SUBROUTINE Herm4(x, y)
!		REAL*8, INTENT(IN) :: x
!		REAL*8, INTENT(OUT) :: y
!		y = x**4 - 6*x**2 + 3
!	END SUBROUTINE
!	
!	SUBROUTINE Herm5(x, y)
!		REAL*8, INTENT(IN) :: x
!		REAL*8, INTENT(OUT) :: y
!		y = x**5 - 10*x**3 + 
!	END SUBROUTINE

	
	END MODULE HERMITE

