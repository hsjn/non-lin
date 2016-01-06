      SUBROUTINE  qrupdate(n,u,v,method,Z,LM)
c------------------------------------------------------------------------
c--- This routine is part of the scheme for updateing  A to A+.
c--- Given QR factorization of matrix A, get the factorization of
c--- Q+*R+ of A+ = Q(R+u*vt) 
C--- LM contains R in upper triangular
c--- LM contains R+ on output

c------------------------------------------------------------------------



      IMPLICIT none
c  INPUT:
      INTEGER n,method
      REAL *8 u(n),v(n)
C  LOCAL:
      INTEGER
     .        i,j,k
      REAL*8 
     .       uip1,mip
c INPUT/OUTPUT:
      REAL*8 LM(N,N),Z(n,n)


      k = n
      DO WHILE (u(k) .eq. 0.0d0 .and. k .gt. 1)
         k = k -1
      ENDDO

      DO i = k-1,1,-1
         uip1 = - u(i+1)
         CALL jackrotate(n,i,u(i),uip1,method,Z,LM)
         IF(u(i) .eq. 0.0D0)THEN
            u(i) = ABS(u(i+1))
         ELSE
            u(i) = SQRT(u(i)**2 + u(i+1)**2)
         ENDIF
      ENDDO

      DO j =1,n
       LM(1,j) = LM(1,j)+U(1)*v(j)
      ENDDO
      
      DO i =1,k-1
        mip = -LM(i+1,i)
         CALL jackrotate(n,i,LM(i,i),mip,method,Z,LM)
      ENDDO
      

      RETURN
      END

