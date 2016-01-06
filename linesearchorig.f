

      SUBROUTINE LINESEARCH(N,XCUR,FCUR,FVEC,GRADIENT,P,SCALE,NOSCALE,
     &                 Sf,MAXSTEP,STEPTOL,RETCODE,XNEW,FNEW,FVp,
     &                 MAX_TAKEN)
C----------------------------------------------------------------------
C---FIND THE APPROXIMATE MINIMUM ALONG A DESCENT LINE FROM XCUR.
C---THAT IS,GIVEN A VECTOR P SUCH THAT GRADIENTTRANSPOSE*P<0
C---FIND A LAMBDA SUCH THAT XNEW=XCUR+LAMBDA*P YIELDS
C---FNEW<FCUR+ALPHA*LAMBDA*GTRANSPOSE*P,WHERE FNEW AND FCUR ARE THE
C---OBJECTIVE FUNCTION EVALUATED AT XNEW AND XCUR RESPECTIVELY.
C---INPUT
C  N                  NUMBER OF ELEMENTS IN XCUR,P,SCALE,XNEW
C  XCUR(I)            I=1,2...N CURRENT SOLUTION POINT
C  FCUR               VALUE OF OBJECTIVE FUNCTION AT XCUR
C  FVEC                 SUBROUTINE NAME WHICH CALCULATES FNEW
C  GRADIENT(I)        GRADIENT OF OBJECTIVE FUNCTION AT XCUR
C  P(I)               (ANY) DESCENT DIRECTION FROM XCUR
C  SCALE(I)           DIAGONAL SCALING MATRIX
C  MAXSTEP            MAX ALLOWABLE STEP SIZE
C  STEPTOL            STEP TOLERANCE FOR DISTINGUISHING XNEW,XCUR
C  NOSCALE            IF .TRUE. DONT USE SCALE,
C                     IF .FALSE. SCALE CALCULATIONS USING SCALE
C  Sf(n)              Typ size of F(n)
C
C---OUTPUT:
C  RETCODE              INTEGER VARIABLE WITH MEANING AS FOLLOWS
C                       RETCODE=0   VALID XNEW FOUND
C                               1   FAILED TO FIND A SATISFACTORY XNEW,
C                                   SUFFICIENTLY DISTINC FROM XCUR
C  XNEW(I)              I=1,2..N THE NEW SOLUTION POINT
C  FNEW                 VALUE OF OBJECTIVE FUNCTION AT XNEW
C  MAX_TAKEN            LOGICAL VARIABLE =.TRUE. IF MAX STEP
C                       WAS TAKEN
C
C
C------------------------------------------------------------------HSJ
      
      USE copy
      IMPLICIT NONE
      INTEGER N,RETCODE,J
      REAL*8  SCALE(n),
     &     MAXSTEP,STEPTOL,ALPHA,STEPLENGTH,
     &     RELLENGTH,INITSLOPE,MINLAMBDA,LAMBDA,ABS,
     &     DUMY1,DUMY2,DUMY3,LAMBDATEMP,LAMBDAPREV,
     &     DISC,A,B,LAMBDASQINV,LAMBDAPREVSQINV,
     &     XCUR(n),Sf(n),FVp(n),
     &         GRADIENT(n),P(n),XNEW(n),DDDD,
     &         FCUR,FNEW,FNEWPREV,SDOT
      LOGICAL MAX_TAKEN,NOSCALE
      EXTERNAL FVEC
C
C

      MAX_TAKEN=.FALSE.
      RETCODE=2
      ALPHA=1.E-4
C---GET LENGTH OF DESCENT VECTOR:
      CALL EUCLIDNORM(P,SCALE,N,NOSCALE,STEPLENGTH)
      IF(STEPLENGTH .GT. MAXSTEP)THEN !DECREASE LENGTH IF TOO LARGE
          DDDD= MAXSTEP/STEPLENGTH
          DO J=1,N
              P(J)=P(J)*DDDD
          ENDDO
          STEPLENGTH=MAXSTEP
      ENDIF
      INITSLOPE=SDOT(GRADIENT,P,N)
      RELLENGTH=0.0
      DO J=1,N
          DUMY1=1.0
          IF( .NOT. NOSCALE)DUMY1=1./SCALE(J)
          DUMY1=MAX(ABS(XCUR(J)),DUMY1)
          RELLENGTH=MAX(RELLENGTH,ABS(P(J))/DUMY1)
      ENDDO
      MINLAMBDA=STEPTOL/RELLENGTH
      LAMBDA=1.d0
C----------------------------------------------------------------------
C---CHECK IF XNEW=XCUR+LAMBDA*P IS SATISFACTORY. IF NOT GENERATE NEW LAMBDA
C---UNTIL IT IS O.K.:
C----------------------------------------------------------------------
      DO WHILE (RETCODE .GE. 2)
          DO J=1,N
              XNEW(J)=XCUR(J)+LAMBDA*P(J)
          ENDDO
           CALL FVEC(N,XNEW,SF,FVp,FNEW)
           print *,'FNEW,fcur ,lambda,initslope',
     .            FNEW,fcur,lambda,initslope

          IF(FNEW .LT. FCUR +ALPHA*LAMBDA*INITSLOPE)THEN
              RETCODE=0   !NEW STEP IS OK
              IF(LAMBDA .EQ. 1.0 .AND. (STEPLENGTH .GT. 0.99*MAXSTEP))
     &                                                  MAX_TAKEN=.TRUE.
          ELSE IF(LAMBDA .LT. MINLAMBDA)THEN !XNEW,XCUR NOT SUFFICIENTLY
                  RETCODE=1                  !DISTINCT
                  CALL DCOPY(N,XCUR,1,XNEW,1)
          ELSE  !REDUCE LAMBDA
                  IF(LAMBDA .EQ. 1.0)THEN !QUADRATIC FIT ON FIRST TRY
                      LAMBDATEMP=-INITSLOPE/(2.*(FNEW-FCUR-INITSLOPE))
                  ELSE  !CUBIC FIT THEREAFTER
                      DUMY1=FNEW-FCUR-LAMBDA*INITSLOPE
                      DUMY2=FNEWPREV-FCUR-LAMBDAPREV*INITSLOPE
                      DUMY3=1./(LAMBDA-LAMBDAPREV)
                      LAMBDASQINV=1./(LAMBDA*LAMBDA)
                      LAMBDAPREVSQINV=1./(LAMBDAPREV*LAMBDAPREV)
                      A=DUMY3*(LAMBDASQINV*DUMY1-LAMBDAPREVSQINV*DUMY2)
                      B=DUMY3*(-LAMBDAPREV*LAMBDASQINV*DUMY1
     &                                +LAMBDA*LAMBDAPREVSQINV*DUMY2)
                      DISC=B**2-3.*A*INITSLOPE
                      IF(A .EQ. 0.0)THEN
                          LAMBDATEMP=-INITSLOPE/(2.*B)
                      ELSE
                          LAMBDATEMP=(-B+SQRT(DISC))/(3.*A)
                      ENDIF
                      IF(LAMBDATEMP .GT. 0.5*LAMBDA)LAMBDATEMP
     &                                                       =0.5*LAMBDA
                  ENDIF
                  LAMBDAPREV=LAMBDA
                  FNEWPREV=FNEW
                  IF(LAMBDATEMP .LE. 0.1*LAMBDA)THEN
                      LAMBDA=0.1*LAMBDA
                  ELSE
                      LAMBDA=LAMBDATEMP
                  ENDIF
          ENDIF
      ENDDO
C
C
      RETURN
      END
