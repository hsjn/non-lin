      SUBROUTINE condest(n,M,M2,est)
c--------------------------------------------------------------------------
cEstimate L1 condition number of upper triangular matrix R
cstored in upper triangle of M and diagonal M2
c------------------------------------------------------------------HSJ
      IMPLICIT none
      INTEGER 
     .       n
      REAL *8
     .       M(n,n)
c LOCAL:
       INTEGER j,i
       REAL*8
     .       est,temp,tempm,xnorm,x(n),p(n),pm(n),
     .       xp,xm

!OUTPUT:
      REAL*8
     .     M2(n)


      est =ABS(M2(1))
      DO j =2,n
         temp =0.0d0
         DO i =1,j-1
            temp =temp +ABS(M(i,j))
         ENDDO
         temp =temp +ABS(M2(j))
         est = MAX(temp,est)
      ENDDO
      x(1)=1.D0/M2(1)
      DO i =2,n
         p(i) = M(1,i)*x(1)
      ENDDO
      DO j=2,n
         xp =(1.d0-p(j))/M2(j)
         xm = -(1.d0+p(j))/M2(j)
         temp = ABS(xp)
         tempm = ABS(xm)
         DO i = j+1,n
            pm(i) = p(i) +M(j,i)*xm
            tempm = tempm + ABS(pm(i)/M2(i))
            p(i) =p(i) +M(j,i)*xp
            temp =temp +ABS(p(i)/M2(i))
         ENDDO
         IF(temp .ge. tempm)THEN
            x(j) = xp  ! ej = 1
         ELSE
            x(j) = xm
            DO i = j+1,n
               p(i) = pm(i)
            ENDDO
         ENDIF
      ENDDO
      xnorm =0.0D0
      DO j =1,n
         xnorm =xnorm +ABS(x(j))
      ENDDO
      est =est /xnorm
      CALL rsolve(n,M,M2,x)
      xnorm =0.0D0
      DO j =1,n
         xnorm =xnorm +ABS(x(j))
      ENDDO

      est=est*xnorm

      RETURN
      END




