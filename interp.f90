!c     this file contains the quasi-vanilla and slightly modified
!c     versions of FORTRAN Numerical Recipe routines

!c     SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
!c     SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
!c     SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)

!c     two dimentional interpolation of a tabulated function, YA=F(X1A,X2A)
!c     SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
!c     SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)

!c     MODIFICATIONS from the vanilla Numerical Recipes: all of the
!c     following routines have been trivially modified in two ways from
!c     the original vanilla; (1) the declarations are now IMPLICIT NONE
!c     or IMPLICIT DOUBLE PRECISION; (2) communication of errors which
!c     may occur in the form of WRITE statements with hard STOPs; some of
!     the routines have been further modified (see comments in each
!c     routine for details)


!c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)

!c     perform polynomial interpolation; given arrays XA and YA of length
!c     N and a value X, this routine returns a value Y and an error
!c     estimate DY. if P(X) is the polynomial of degree N-1 such that
!c     YA=P(XA), then the returned value is Y=P(X)

!c..............................................................................

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
!c      check division by zero error; finite difference fail
          IF(DEN.EQ.0.) THEN
           WRITE(6,*) 'ERROR(polint): division by zero error due to'
           WRITE(6,*) 'ill-posed finite differencing'
           WRITE(6,*) '*** terminating program  ***'
           STOP
          END IF
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


!c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)

!c     computes 2D second derivatives of a grid for interpolation; given
!c     an M by N tabulated function YA, and tabulated independent
!c     variables X1A (an array of length M) and X2A (an array of length
!c     N), this routine constructs one dimensional natural cubic splines
!c     of the rows of YA and returns the second derivatives in the table
!c     Y2A (of dimensions M by N); the Y2A table is employed by the
!c     routine 'SPLIN2' (below), which preforms the 2D cubic spline

!c     this routine is only called once and must be called in order for
!c     routine 'SPLIN2' to work

!c..............................................................................

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=600)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)
      DO 13 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL SPLINE(X2A,YTMP,N,1.D30,1.D30,Y2TMP)
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
13    CONTINUE
      RETURN
      END


!c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)

!c     perform 2D interpolation of a table; given the arrays of the
!c     independent variables X1A (of length M) and X2A (of length N), and
!c     given the function YA (dimension M by N) that tablulates a
!c     monotonic function of X1A and X2A, and given the array Y2A (of
!c     dimension M by N), which contains the second derivatives of YA and
!c     is the output of routine 'SPLINE' (above), and given a single
!c     value of the independent variable X1 and a single value of the
!c     independent variable X2, this routine returns the value Y using
!c     two dimensional cubic spline interpolation of the table Y2A; Note-
!c     the array X1A and X2A arrays must be in ascending order or the
!c     indices of X1A that bracket X1 and of X2A that bracket X2 will not
!c     be correctly determined

!c..............................................................................

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=600)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN), &
               Y2TMP(NN),YYTMP(NN)
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
11      CONTINUE
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
12    CONTINUE
      CALL SPLINE(X1A,YYTMP,M,1.D30,1.D30,Y2TMP)
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
      RETURN
      END


!c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)

!c     compute the 1D second derivatives of the table Y which is a
!c     function of the array X to set up interpolation by routine
!c     'SPLINE' (below); given arrays X and Y of length N containing a
!c     tabulated function i.e y=f(x), with X monotonically increasing,
!c     and given values for YP1 and YPN, which are the first derivatives
!c     at the points N=1 and N, respectively, this routine returns the
!c     array Y2 of length N which contains the second derivatives of the
!c     Y at the tabulated points X;

!c     regarding the choice of values for YP1 and YPN; if YP1 and/or YPN
!c     are set larger than 1.0E30, the routine automatically sets the
!c     boundary condition for a "natural spline" (one with nulled second
!c     derivatives at the boundary)

!c     this routine is only called once and must be called in order for
!c     routine 'SPLINE' to work

!c..............................................................................

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=600)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99D30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
           /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99D30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END


!c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)

!c     perform 1D interpolation of a table; given the arrays XA and YA of
!c     length N, where XA is the independent variable and YA tablulates a
!c     monotonic function, and given the array Y2A, which contains the
!     second derivatives of YA and is the output of routine 'SPLINE'
!c     (above), and given a single value of the independent variable X,
!c     this routine returns the value Y using cubic spline interpolation
!c     of the table YA; Note- the array XA must be in ascending order or
!c     the indices of XA that bracket X will not be correctly determined

!c..............................................................................

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
!c     check for bad XA input
      IF (H.EQ.0.) THEN
       WRITE(6,*) 'ERROR(splint): degenerate indices for bracketing'
       WRITE(6,*) 'dependent variable on table for 1D interpolation'
       WRITE(6,*) 'lower index = ',KLO,' XA(KLO) = ', XA(KLO)
       WRITE(6,*) 'upper index = ',KHI,' XA(KHI) = ', XA(KHI)
       WRITE(6,*) '*** terminating program ***'
       STOP
      END IF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+ &
           ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
