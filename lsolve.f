C
C
      SUBROUTINE LSOLVE(N,B,L,Y)
C---------------------------------------------------------------------
C---SOLVE L*Y=B
C---WHERE L IS LOWER TRIANGULAR
C--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I,J,K
      REAL *8 L(N,N),Y(*),B(*),SUMD
C
C
      Y(1)=B(1)/L(1,1)
      K=1
      DO I=2,N
        SUMD=0.0D0
        DO J=1,I-1
          K=K+1
          SUMD=SUMD+L(i,j)*Y(J)
        ENDDO
        K=K+1
        Y(I)=(B(I)-SUMD)/L(i,i)
      ENDDO
C
C
      RETURN
      END
