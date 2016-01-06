      SUBROUTINE  jackrotate(n,i,alpha,beta,method,Z,LM)
c------------------------------------------------------------------------
c--- premultiply M,Z by Jacobi rotation matrix J(i,i+1,alpha,beta)
c--- this uses upper triangular LM plus the subdiagonal(of M).
c--- for diagoanl elment (i,i) the subdiagoanl is (i+1,i). this subdiagonal
c--  element is found in postion (n*(n+1))/2 + (i*(i+1))/2
c--  in vector LM (it waas put there in routine loadRT)

c----------------------------------------------------------HSJ-----------



      IMPLICIT none
c  INPUT:
      INTEGER n,method,i
      REAL *8 alpha,beta
C  LOCAL:
      INTEGER
     .        j
      REAL*8 
     .      den,c,s ,y,w
c INPUT/OUTPUT:
      REAL*8 LM(N,N),Z(n,n)


      IF(alpha .eq. 0.0d0)THEN
         c = 0.0D0
         s = sign(1.d0,beta)
      ELSE
         den =SQRT(alpha**2+beta**2)
         c = alpha/den
         s = beta/den
      ENDIF
      DO j =i,n
         y = LM(i,j)
         w = LM(i+1,j)
         LM(i,j) = c*y -s*w
         LM(i+1,j) = s*y + c*w
      ENDDO
      
      IF(method .eq. 1)THEN
         DO j=1,n
            y=Z(i,j)
            w = Z(i+1,j)
            Z(i,j) = c*y -s*w
            Z(i+1,j) = s*y +c*w
         ENDDO
      ENDIF
      
      RETURN
      END

