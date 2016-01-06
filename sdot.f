C
C
C
C
      REAL*8  FUNCTION SDOT(A,B,N)
C----------------------------------------------------------------------
C---CALCULATE THE DOT PRODUCT OF VECTORS A AND B IN DOUBLE PRECISION
C---HERE WE UNROLL THE DO LOOP FOR BETTER PERFORMANCE
C--------------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      REAL*8  SUMD,A(*),B(*)
      INTEGER M,N,MP1,I
C
C
      SUMD=0.0D0
      M=MOD(N,5)
      IF(M .EQ. 0)GO TO 40
      DO I=1,M
         SUMD=SUMD+A(I)*B(I)
      ENDDO
      IF(N .LT. 5)GO TO 50
   40 MP1=M+1
      DO I=MP1,N,5
         SUMD=SUMD+A(I)*B(I)+A(I+1)*B(I+1)
     &          +A(I+2)*B(I+2)+A(I+3)*B(I+3)
     &          +A(I+4)*B(I+4)
      ENDDO
   50 SDOT=SUMD
C
C
      RETURN
      END
