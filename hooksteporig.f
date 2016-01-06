
          SUBROUTINE HOOKSTEP(N,GRADIENT,L,HESSIAN,STEP,SCALE,NOSCALE,
     &                EPSILOND,DELTA,MU,DELTAPREV,PHI,PHIP,FIRSTHOOK,
     &                PHIPINIT,S,NEWTON_STEP,STEPLENGTH,IOUNIT)
C---------------------------------------------------------------------
C
C
      USE copy
      IMPLICIT none
      INTEGER N,I,J,K,N2,IOUNIT
      REAL *8 TEMPVEC(n),S(*),L(N,N),HESSIAN(*),STEP(*),GRADIENT(*),
     &        EPSILOND
      REAL *8 HI,LOW,STEPLENGTH,MU,DELTAPREV,PHI,PHIP,SCALE(*),
     &      MUUP,DELTA,PHIPINIT,DUMY,MULOW,ADDMAX,
     &      STEPLS,POSDEF
      LOGICAL NEWTON_STEP,FIRSTHOOK,NOSCALE,DONE
C
C
      N2=2*N
      HI=1.5
      LOW=0.75
      IF(STEPLENGTH .LE. HI*DELTA)THEN !SET S TO NEWTON STEP
          NEWTON_STEP=.TRUE.
          CALL DCOPY(N,STEP,1,S,1)
          MU=0.0
          DELTA=MIN(DELTA,STEPLENGTH) 
         IF(IOUNIT .NE. 0)THEN
             WRITE(IOUNIT,*)'NEWTON STEP TAKEN in HOOKSTEP'
             WRITE(IOUNIT,2)STEPLENGTH,DELTA
    2         FORMAT(2X,'STEPLENGTH,DELTA =',2(1PE14.5))
         ENDIF
      ELSE
          NEWTON_STEP=.FALSE.
          IF(MU .GT. 0.0)
     &        MU=MU-(PHI+DELTAPREV)*((DELTAPREV-DELTA)+PHI)/(DELTA*PHIP)
          PHI=STEPLENGTH-DELTA
          IF(FIRSTHOOK)THEN
              FIRSTHOOK=.FALSE.
              DO J=1,N
                  DUMY=1.0
                  IF(.NOT. NOSCALE)DUMY=SCALE(J)
                  TEMPVEC(J)=STEP(J)*DUMY*DUMY
              ENDDO

              CALL LSOLVE(N,TEMPVEC,L,TEMPVEC)
              CALL EUCLIDNORM(TEMPVEC,SCALE,N,.TRUE.,PHIPINIT)
              PHIPINIT=-PHIPINIT*PHIPINIT/STEPLENGTH
          ENDIF
          MULOW=-PHI/PHIPINIT
          DO J=1,N
              DUMY=1.0
              IF(.NOT. NOSCALE)DUMY=1./SCALE(J)
              TEMPVEC(J)=GRADIENT(J)*DUMY
          ENDDO
          CALL EUCLIDNORM(TEMPVEC,SCALE,N,.TRUE.,MUUP)
          MUUP=MUUP/DELTA
          DONE=.FALSE.
          DO WHILE (.NOT. DONE)
              IF((MU .LT. MULOW) .OR. (MU .GT. MUUP))
     &               MU=MAX(SQRT(MULOW*MUUP),1.E-3*MUUP)
              DO I=1,N
                  K=(I*(N2-I+3))/2-N
                  DUMY=1.0
                  IF(.NOT. NOSCALE)DUMY=SCALE(I)
                  HESSIAN(K)=HESSIAN(K)+MU*DUMY*DUMY
              ENDDO
              POSDEF=0.0
              CALL CHOLDECOMP(N,HESSIAN,POSDEF,EPSILOND,L,ADDMAX)
              CALL CHOLSOLVE(N,GRADIENT,L,S)
C             RESTORE HESSIAN :
              DO I=1,N
                  K=(I*(N2-I+3))/2-N
                  DUMY=1.0
                  IF(.NOT. NOSCALE)DUMY=SCALE(I)
                  HESSIAN(K)=HESSIAN(K)-MU*DUMY*DUMY
              ENDDO
              CALL EUCLIDNORM(S,SCALE,N,NOSCALE,STEPLS)
              PHI=STEPLS-DELTA
              DO J=1,N
                  DUMY=1.0
                  IF(.NOT. NOSCALE)DUMY=SCALE(J)
                  TEMPVEC(J)=S(J)*DUMY*DUMY
              ENDDO
              CALL LSOLVE(N,TEMPVEC,L,TEMPVEC)
              CALL EUCLIDNORM(TEMPVEC,SCALE,N,.TRUE.,PHIP)
              PHIP=-PHIP*PHIP/STEPLS
              IF(PHIP .EQ. 0.0)DONE=.TRUE. !HSJ 5/2/94 CANT DO ANY BETTER
              IF( .NOT. DONE)THEN          !SO ACCEPT THIS AND RETURN
                  IF((STEPLS .GE. LOW*DELTA) .AND. (STEPLS .LE.
     &               HI*DELTA) .OR. (MUUP-MULOW .LE. 0.0))THEN
                      DONE=.TRUE.    !ACCEPT S AS THE STEP
                  ELSE        !S NOT ACCEPTABLE,CALCULATE NEW MU,MULOW,MUUP:
                      MULOW=MAX(MULOW,MU-PHI/PHIP)
C                     IF(PHI .GT. 0.0)MULOW=MU
                      IF(PHI .LT. 0.0)MUUP=MU
                      MU=MU-(STEPLS*PHI)/(DELTA*PHIP)
                  ENDIF
              ENDIF
          ENDDO
C         IF(IOUNIT .NE. 0)THEN
C             WRITE(IOUNIT,1)MU
C   1         FORMAT(2X,'HOOKSTEP TAKEN,MU=',2X,1PE12.4)
C             WRITE(IOUNIT,2)STEPLENGTH,DELTA
C         ENDIF
      ENDIF
C
C
      RETURN
      END
