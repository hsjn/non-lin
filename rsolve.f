      SUBROUTINE rsolve(n,M,M2,B)
c--------------------------------------------------------------------------
c  Solve Rx = b 
c  R is stored in upper triangle of M and diagonal of R is in M2
c on output B contains solution x
c------------------------------------------------------------------HSJ
      IMPLICIT none
      INTEGER 
     .       n
      REAL *8
     .       M(n,n)
c LOCAL:
       INTEGER j,i
       REAL*8
     .       sumd

!OUTPUT:
      REAL*8
     .     M2(n),B(n)


      B(n) =B(n)/M2(n)
      DO i =n-1,1,-1
         sumd =0.0D0
         DO j = i+1,n
            sumd =sumd+M(i,j)*B(j)
         ENDDO
         B(i) = (B(i) - sumd)/M2(i)
      ENDDO

      RETURN
      END

