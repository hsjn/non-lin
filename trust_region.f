C
C
C
      SUBROUTINE TRUST_REGION(N,XCUR,FCUR,FVEC,GRADIENT,L,S,SCALE,
     &         NOSCALE,Sf,
     &         NEWTON_STEP,STEPLENGTH_MAX,STEP_TOL,STEP_TYPE,HESSIAN,
     &         DELTA,RETCODE,XPREV,FPREV,XNEW,FNEW,MAX_TAKEN)
C-----------------------------------------------------------------------
C   GIVEN THE STEP S(I),I=1,2..N,PRODUCED BY THE HOOKSTEP OR DOGLEG
C   METHODS,DECIDE IF XNEW=XCUR+S IS ACCEPTABLE AS THE NEXT ITERATE.
C   IF SO THEN ADJUST THE INITIAL TRUST REGION FOR THE NEXT ITERATION
C   APPROPRIATELY. IF NOT THEN DECREASE OR INCREASE THE TRUST REGION
C   RADIUS FOR THE CURRENT ITERATION.
C---INPUT
C  N                  NUMBER OF ELEMENTS IN XCUR,GRADIENT,S,SCALE,XNEW
C  XCUR(I)            I=1,2...N CURRENT SOLUTION POINT
C  FCUR               VALUE OF OBJECTIVE FUNCTION AT XCUR
C  FVEC                 SUBROUTINE NAME WHICH CALCULATES FNEW
C  GRADIENT(I)        GRADIENT OF OBJECTIVE FUNCTION AT XCUR
C  S(I)               I=1,2..N,THE STEP TO BE TESTED.
C  L(N,N))            CHOLESKY DECOMPOSITION OF HESSIAN
C  SCALE(I)           DIAGONAL SCALING MATRIX
C  NOSCALE            IF .TRUE. DONT USE SCALE,
C                     IF .FALSE. SCALE CALCULATIONS USING SCALE
C  NEWTON_STEP        BOOLEAN,=.TRUE. ONLY IF FULL NEWTON STEP WAS TAKEN
C  STEPLENGTH_MAX     MAX ALLOWABLE STEP SIZE
C  STEP_TOL           STEP TOLERANCE LEVEL
C  STEP_TYPE          =1  FOR HOOKSTEP,=2 FOR DOUBLE DOGLEG
C  HESSIAN           in symmetric storage mode

C---INPUT AND OUTPUT:
C  DELTA                THE TRUST REGION RADIUS
C  RETCODE              INTEGER VARIABLE WITH MEANING AS FOLLOWS
C                       RETCODE=0   VALID XNEW FOUND,DELTA SET TO NEW VALUE
C                               1   FAILED TO FIND A SATISFACTORY XNEW,
C                                   SUFFICIENTLY DISTINC FROM XCUR
C                                   BECAUSE LENGTH OF XNEW-XCUR IS LESS
C                                   THAN STEP TOL. THIS MEANS THAT NO
C                                   FURTHER PROGRESS TOWARD THE SOLUTION
C                                   IS POSSIBLE. IN ALMOST ALL CASES THIS
C                                   MEANS THAT XCUR IS AN ACCEPTABLE SOLUTION.
C                               2   FNEW IS TOO LARGE. THIS MEANS THAT
C                                   THE TRUST REGION RADIUS WAS SO LARGE
C                                   THAT THE DESCENT STEP TOTALLY OVERSHOT
C                                   THE SOLUTION POINT. WE HAVE TO DECREASE
C                                   THE TRUST REGION RADIUS DELTA AND
C                                   TRY AGAIN
C                               3   FNEW IS SUFFICIENTLY SMALL BUT
C                                   THE PROBABILITY OF TAKING A LONGER
C                                   SUCCESSFUL STEP IS LARGE. HENCE
C                                   CONTINUE THE CURRENT ITERATION WITH
C                                   A LARGER TRUST REGION RADIUS DELTA.
C  XPREV
C  FPREV
C---OUTPUT:
C  XNEW(I)              I=1,2..N THE NEW SOLUTION POINT
C  FNEW                 VALUE OF OBJECTIVE FUNCTION AT XNEW
C  MAX_TAKEN            LOGICAL VARIABLE =.TRUE. IF MAX STEP
C                       WAS TAKEN
C
C
C------------------------------------------------------------------HSJ
C
C
      USE copy
      IMPLICIT NONE
      INTEGER N,RETCODE,I,J,STEP_TYPE,
     &        N2,IHIJ,IHII,ILJI
      REAL*8 SCALE(*),
     &     STEPLENGTH_MAX,STEP_TOL,ALPHA,STEPLENGTH,
     &     RELLENGTH,INITSLOPE,ABS,
     &     DELTAF,DELTAF_PREDICTED,
     &     DUMY,DELTA,SUM
      REAL *8 L(N,N),XCUR(*),GRADIENT(*),SDOT,
     &      S(*),HESSIAN(*),Sf(N),FVp(n),
     &      XPREV(*),XNEW(*),FNEW,FCUR,FPREV
      LOGICAL MAX_TAKEN,NOSCALE,NEWTON_STEP
      EXTERNAL FVEC
C
C
      N2=2*N
      MAX_TAKEN=.FALSE.
      ALPHA=1.E-04
      CALL EUCLIDNORM(S,SCALE,N,NOSCALE,STEPLENGTH)
      DO J=1,N
        XNEW(J)=XCUR(J)+S(J)
      ENDDO
C---GET THE FUNCTION VALUE AT XNEW,DERIVATIVES ARE NOT USED:
      CALL FVEC(N,XNEW,Sf,FVp,FNEW)
      DELTAF=FNEW-FCUR
      INITSLOPE=SDOT(GRADIENT,S,N)
      IF(RETCODE .NE. 3)FPREV=0.0
      IF((RETCODE .EQ. 3) .AND. ((FNEW .GE. FPREV) .OR. (DELTAF .GT.
     &                                         ALPHA*INITSLOPE)))THEN
          RETCODE=0
          CALL DCOPY(N,XPREV,1,XNEW,1)
          FNEW=FPREV
          DELTA=DELTA*0.5
C     ELSE IF(DELTAF .GE. ALPHA*INITSLOPE .AND. NEWTON_STEP .EQ.
C    &                       .FALSE. )THEN   !FNEW IS TOO LARGE
      ELSE IF(DELTAF .GE. ALPHA*INITSLOPE)THEN
      RELLENGTH=0.0    !ABOVE IS MODIFIED VERSION
          DO J=1,N
              DUMY=1.0
              IF( .NOT. NOSCALE)DUMY=1./SCALE(J)
              DUMY=MAX(ABS(XNEW(J)),DUMY)
              RELLENGTH=MAX(RELLENGTH,ABS(S(J))/DUMY)
          ENDDO
          IF(RELLENGTH .LT. STEP_TOL)THEN !XNEW-XCUR TOO SMALL
              RETCODE=1
              CALL DCOPY(N,XCUR,1,XNEW,1)
          ELSE
C---REDUCE DELTA,CONTINUE GLOBAL STEP
              RETCODE=2
              DUMY=-(INITSLOPE*STEPLENGTH)/(2.0*(DELTAF-INITSLOPE))
              IF(DUMY .LT. 0.1*DELTA)THEN
                  DELTA=0.1*DELTA
              ELSE IF(DUMY .GT. 0.5*DELTA)THEN
                  DELTA=0.5*DELTA
              ELSE
                  DELTA=DUMY
              ENDIF
          ENDIF
      ELSE !FNEW IS SUFFICIENTLY SMALL
          DELTAF_PREDICTED=INITSLOPE
          IF(STEP_TYPE .EQ. 1)THEN
C---HOOKSTEP CALCULATE ST*HESSIAN*S:
              DO I=1,N
                  SUM=0.0
                  DO J=I+1,N
                      IHIJ=((I-1)*(N2-I)+2*J)/2
                      SUM=SUM+HESSIAN(IHIJ)*S(J)
                  ENDDO
                  IHII=(I*(N2-I+3))/2-N
                  DUMY=S(I)*(0.5*HESSIAN(IHII)*S(I)+SUM)
                  DELTAF_PREDICTED=DELTAF_PREDICTED+DUMY
              ENDDO
          ELSE
C---DOGLEG CALCULATE ST*L*LT*S:
              DO I=1,N
                  SUM=0.0
                  DO J=I,N
                      SUM=SUM+L(j,i)*S(J)
                  ENDDO
                  DELTAF_PREDICTED=DELTAF_PREDICTED+0.5*SUM*SUM
              ENDDO
          ENDIF
          DUMY=ABS(DELTAF_PREDICTED-DELTAF)
          IF(RETCODE .NE. 2 .AND. ((DUMY .LE. 0.1*ABS(DELTAF))
     &        .OR. (DELTAF .LE. INITSLOPE)) .AND.( NEWTON_STEP .EQV.
     &              .FALSE.) .AND. (DELTA .LE. 0.99*STEPLENGTH_MAX))THEN
C---DOUBLE DELTA AND CONTINUE THE ITERATION
              RETCODE=3
              CALL DCOPY(N,XNEW,1,XPREV,1)
              FPREV=FNEW
              DELTA=MIN(2.*DELTA,STEPLENGTH_MAX)
          ELSE
C---ACCEPT XNEW AS NEW ITERATE,CHOOSE NEW TRUST REGION RADIUS DELTA:
              RETCODE=0
              IF(STEPLENGTH .GT. .99*STEPLENGTH_MAX)MAX_TAKEN=.TRUE.
              IF(DELTAF .GE. 0.1*DELTAF_PREDICTED)THEN
                  DELTA=0.5*DELTA
              ELSE IF(DELTAF .LE. 0.75*DELTAF_PREDICTED)THEN
                  DELTA=MIN(2.*DELTA,STEPLENGTH_MAX)
              ENDIF
          ENDIF
      ENDIF
C
C
      RETURN
      END
