!     this file contains the quasi-vanilla and slightly modified
!     versions of FORTRAN Numerical Recipe routines

!c     DOUBLE PRECISION FUNCTION ZBRENT(FUNC,X1,X2,TOL)

!c     MODIFICATIONS from the vanilla Numerical Recipes: all of the
!c     following routines have been trivially modified in two ways from
!c     the original vanilla; (1) the declarations are now IMPLICIT NONE
!c     or IMPLICIT DOUBLE PRECISION; (2) communication of errors which
!c     may occur in the form of WRITE statements with hard STOPs; some of
!c     the routines have been further modified (see comments in each
!c     routine for details)

!c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      FUNCTION ZBRENT(FUNC,X1,X2,TOL)

!c     return the root of the function FUNC using Brent's method where
!c     the root of the function is assumed lie between X1 and X2; the
!c     root is refined until its accuracy is TOL and is returned as
!c     ZBRENT

!c..............................................................................

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ITMAX=100,EPS=3.E-16)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
!c     error: Root not bracketed 
      IF(FB*FA.GT.0.) THEN
       WRITE(6,*) 'ERROR(zbrent): root of FUNC not bracketed.'
       WRITE(6,*) 'lower bracket A = ',A,' FUNC(A) = ',FUNC(A)
       WRITE(6,*) 'upper bracket B = ',B,' FUNC(B) = ',FUNC(B)
       WRITE(6,*) 'the sign of FUNC(A) and FUNC(B) must be'
       WRITE(6,*) 'opposite in order to bracket the root.'
       WRITE(6,*) 'this might be remedied by changing A and/or B'
       WRITE(6,*) 'to broaden your search interval. OR, it may'
       WRITE(6,*) 'indicate an error in your function.'
       WRITE(6,*) '*** terminating program ***'
       STOP
      END IF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
!c     error message for too many iterations
      WRITE(6,*) 'ERROR(zbrent): number of iterations required to'
      WRITE(6,*) 'converge on the root of FUNC exceeded maximum'
      WRITE(6,*) 'allowed.  This can be remedied by increasing'
      WRITE(6,*) 'parameter ITMAX in the file or by relaxing'
      WRITE(6,*) 'the convrgence tolerance TOL'
      WRITE(6,*) '*** terminating program  ***'
      STOP
      ZBRENT=B
      RETURN
      END
