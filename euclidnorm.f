
C
      SUBROUTINE EUCLIDNORM(VEC,SCALE,N,NOSCALE,LENGTH)
C---------------------------------------------------------------------
C---CALCULATE THE EUCLIDEAN NORM OF VECTOR VEC.
C---IF NOSCALE= .FALSE. THEN SCALE VEC BY SCALE BEFORE
C---DOING THE CALCULATION. THIS LENGTH DOES NOT HAVE TO
C---BE VERY PRECISE SO WE DONT USE SDOT HERE.
C------------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      LOGICAL NOSCALE
      INTEGER N,J
      REAL *8 VEC(*),SUM
      REAL*8 LENGTH,SCALE(*)
C
C
      SUM=0.0D0
      IF(.NOT. NOSCALE)THEN
          DO J=1,N
              SUM=SUM+(SCALE(J)*VEC(J))**2
          ENDDO
      ELSE
          DO J=1,N
              SUM=SUM+VEC(J)*VEC(J)
          ENDDO
      ENDIF
      LENGTH=SQRT(SUM)
C
C
      RETURN
      END
