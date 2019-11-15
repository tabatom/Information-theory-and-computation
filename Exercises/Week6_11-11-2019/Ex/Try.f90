	MODULE hermite

	CONTAINS
	RECURSIVE FUNCTION HermitePoly(n) RESULT(hp)

!
! More info
! http://mathworld.wolfram.com/HermitePolynomial.html
!

	REAL hp(n+1)
	REAL hp1(n+1),hp2(n+1)

	IF(n .EQ. 0) THEN
		hp(1)=1.0
		RETURN
	END IF

	IF(n .EQ. 1) THEN
		hp(1)=2.0
		hp(2)=0.0
	ELSE

		hp1(1:n+1)=0.0
		hp1(1:n)=2.0*HermitePoly(n-1)

		hp2(1:n+1)=0.0
		hp2(3:)=2.0*(n-1)*HermitePoly(n-2)

		hp=hp1-hp2

	END IF

	END FUNCTION

	FUNCTION evalHermitePoly(ix,n) RESULT(y)
	INTEGER n,ip
	REAL ix(:),y(size(ix)),h(n+1)

	k=size(ix)

	h=HermitePoly(n)

	y(1:k)=h(n+1)
	ip=1
	DO i=n,1,-1
	  DO j=1,k
	   y(j)=y(j)+h(i)*ix(j)**ip
	  END DO

	  ip=ip+1
	END DO
	END FUNCTION
	END MODULE

	PROGRAM Hermite_
	USE IFqwin
	USE hermite
	INTEGER,parameter :: inum=100
	INTEGER :: pgopen,pgcurs
	INTEGER(4) :: RESULT

	TYPE(qwinfo) :: winfo
	REAL x(inum),y(inum)
	CHARACTER*1 ch
	winfo%TYPE=qwin$MAX

	RESULT=SETWSIZEQQ(QWIN$FRAMEWINDOW,winfo)
	RESULT=ABOUTBOXQQ("Hello\rVersion 1.0"c)

	IF(pgopen('?') .LE. 0) STOP

	CALL pgenv(-5.0,5.0,-200.0,200.0,0,1)
	x(1)=-5.0
	ainc=0.1
	! print *, ainc
	DO i=2,inum
	x(i)=x(i-1)+ainc
	END DO

	   DO i=1,6

	   CALL pgsci(i+1)
	   y=evalHermitePoly(x,i)

	   CALL pgline(inum,x,y)

	   END DO

	CALL pgclos

	END PROGRAM
