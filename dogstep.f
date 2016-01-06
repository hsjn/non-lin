C
          SUBROUTINE DOGSTEP(N,GRADIENT,L,STEP,SCALE,NOSCALE,STEPLENGTH,
     &                       STEPLENGTH_MAX,DELTA,FIRSTDOG,CAUCHYLENGTH,
     &                       ETA,SSD,V,S,NEWTON_STEP,IOUNIT)
C-----------------------------------------------------------------------
C---HERE WE FIND AN APPROXIMATE SOLUTION TO
C    MINMIZE ( GT*S+0.5*ST*L*LT*S) SUBJECT TO THE CONSTRAINT THAT
C    THE (SCALED) STEP LENGTH IS LESS THAN DELTA.
C    THE OUTPUT STEP S(I),I=1,2..N WILL HAVE LENGTH DELTA IF THE
C    FULL NEWTON STEP COULD NOT BE TAKEN (BECAUSE IT IS OUTSIDE
C    THE TRUST REGION). IF THE FULL NEWTON STEP WAS TAKEN,INDICATED
C    BY NEWTON_STEP =TRUE,THEN THE OUTPUT STEP S(I) IS IDENTICAL
C    TO THE INPUT STEP(I).
C---INPUT
C  N                  NUMBER OF ELEMENTS IN STEP,SCALE,GRADIENT
C  GRADIENT(I)        GRADIENT OF OBJECTIVE FUNCTION AT XCUR
C  L(N,N)               CHOLESKY FACTOR OF HESSIAN
C  STEP(I)            I=1,2..N  THE CURRENT STEP VECTOR
C  SCALE(I)           DIAGONAL SCALING MATRIX
C  NOSCALE            IF .TRUE. DONT USE SCALE,
C                     IF .FALSE. SCALE CALCULATIONS USING SCALE
C  STEPLENGTH         MAGNITUDE OF STEP
C  STEPLENGTH_MAX     MAX ALLOWABLE STEP SIZE
C
C
C---INPUT AND OUTPUT:
C  FIRSTDOG           IS SET TO FALSE TO INDICATE THAT THE NEXT
C                     CALL TO THIS SUBROUTINE DOES NOT HAVE TO
C                     CALCULATE THE DOUBLE DOGLEG CURVE
C  DELTA              THE CURRENT TRUST REGION RADIUS
C  CAUCHYLENGTH       LENGTH OF THE CAUCHY STEP
C  ETA
C  SSD(I)             I=1,2..N,THE STEEPEST DESCENT STEP
C  V(I)
C
C
C---OUTPUT:
C  NEWTON_STEP        SET TO TRUE IF FULL NEWTON STEP WAS TAKEN
C                     OTHERWISE SET TO FALSE.
C  S(I)               I=1,2..N  THE NEW STEP VECTOR
C---------------------------------------------------------------------HSJ
C}}}
      USE copy
      IMPLICIT NONE
      INTEGER N,I,J,K,IOUNIT
      REAL*8
     &     STEPLENGTH,STEPLENGTH_MAX,DELTA,CAUCHYLENGTH,ETA,
     &     DUMY,DUMY1,DUMY2,LAMBDA,ALPHA,BETA,SCALE(*),
     &      L(N,N),GRADIENT(*),STEP(*),SSD(*),V(*),S(*),SDOT
      LOGICAL 
     &     NOSCALE,NEWTON_STEP,FIRSTDOG
C
      IF(STEPLENGTH .LE. DELTA)THEN
          NEWTON_STEP=.TRUE.
          CALL DCOPY(N,STEP,1,S,1)
          DELTA=STEPLENGTH
          IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'NEWTON STEP TAKEN'
      ELSE
C---NEWTON STEP TOO LONG,GET S ON DOUBLE DOGLEG CURVE:
          NEWTON_STEP=.FALSE.
          IF(FIRSTDOG)THEN
C---CALCULATE THE DOUBLE DOGLEG CURVE ON FIRST TRY
              FIRSTDOG=.FALSE.  !SET FOR NEXT CALL TO THIS ROUTINE
              IF(.NOT. NOSCALE)THEN
                  DO J=1,N
                      S(J)=GRADIENT(J)/SCALE(J)
                  ENDDO
              ELSE
                  CALL DCOPY(N,GRADIENT,1,S,1)
              ENDIF
              CALL EUCLIDNORM(S,SCALE,N,.TRUE.,ALPHA)
              ALPHA=ALPHA*ALPHA
              BETA=0.0
              DO I=1,N
                  DUMY2=0.0
                  DO J=I,N
                      DUMY=1.0
                      IF(.NOT. NOSCALE)DUMY=SCALE(J)
                      DUMY2=DUMY2+L(j,i)*GRADIENT(J)/(DUMY*DUMY)
                  ENDDO
                  BETA=BETA+DUMY2*DUMY2
              ENDDO
              DUMY=-ALPHA/BETA
              DO J=1,N
                  SSD(J)=DUMY*S(J)   ! SSD is Cauchy step 
              ENDDO
              CAUCHYLENGTH=ALPHA*SQRT(ALPHA)/BETA
              DUMY=SDOT(GRADIENT,STEP,N)
              ETA=0.2+0.8*ALPHA*ALPHA/(BETA*ABS(DUMY))
              DO J=1,N
                  DUMY=ETA
                  IF( .NOT. NOSCALE)DUMY=SCALE(J)*DUMY
                  V(J)=DUMY*STEP(J)-SSD(J)
              ENDDO
              IF(DELTA .EQ. -1.d0)
     .                   DELTA=MIN(CAUCHYLENGTH,STEPLENGTH_MAX)
          ENDIF
          IF(ETA*STEPLENGTH .LE. DELTA)THEN
C---PARTIAL LENGTH STEP IN NEWTON DIRECTION:
              DUMY=DELTA/STEPLENGTH
              DO J=1,N
                  S(J)=DUMY*STEP(J)
              ENDDO
               IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'PARTL NEWTON STEP TAKEN'
          ELSE IF(CAUCHYLENGTH .GE. DELTA)THEN
C---TAKE STEP IN STEEPEST DESCENT DIRECTION:
              DUMY=DELTA/CAUCHYLENGTH
              DO J=1,N
                  DUMY1=1.0
                  IF(.NOT. NOSCALE)DUMY1=1./SCALE(J)
                  S(J)=DUMY*SSD(J)*DUMY1
              ENDDO
               IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'STEPST DESCT STEP TAKEN'
          ELSE
C---GET THE CONVEX COMBINATION THAT HAS SCALED LENGTH DELTA:
              DUMY=SDOT(V,SSD,N)
              DUMY1=SDOT(V,V,N)
              LAMBDA=(-DUMY+SQRT(DUMY**2-DUMY1*
     &                             (CAUCHYLENGTH**2-DELTA**2)))/DUMY1
              DO J=1,N
                  DUMY=1.
                  IF(.NOT. NOSCALE)DUMY=1./SCALE(J)
                  S(J)=DUMY*(SSD(J)+LAMBDA*V(J))
              ENDDO
                  IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'DOGLEG STEP TAKEN'
          ENDIF
      ENDIF
C
C
      RETURN
      END