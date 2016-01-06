      SUBROUTINE  broyfac(n,xc,xp,fc,fp,eta,Sx,Sf,Z,M,M2)
c------------------------------------------------------------------------
c--- Factored form of Broydens update for Jacobian
c--- INPUT:
c--- INPUT/OUTPUT:
c--- z(n,n)  contains Q transpose of A on input
c            Q+, transpose of A+ on output
c--- LM( )   LM contains R (factor of matrix a)
c            on input and R+, factor of matrix A+ on output
c            
c -- M2(n)   contains diagonals of M
c------------------------------------------------------------------------



      IMPLICIT none
c  INPUT:
      INTEGER n
      REAL *8 xc(n),xp(n),fc(n),fp(n),Sx(n),Sf(n),eta
C  LOCAL:
      INTEGER
     .        i,j
      LOGICAL skipupdate
      REAL*8
     .       s(n),t(n),w(n),sumd,denom
c INPUT/OUTPUT:
      REAL*8 M(N,N),Z(n,n),M2(n)



      DO i =1,n
         M(i,i) = M2(i)
          s(i) = xp(i)-xc(i)
      ENDDO
      skipupdate = .true.

      DO i =1,n
         sumd =0.0D0
         DO j = i,n
            sumd =sumd +M(i,j)*s(j)
         ENDDO
         t(i) = sumd   !t = R*s  R in upper triangle of M
      ENDDO
      DO i =1,n
         sumd =0.0D0
         DO j = 1,n
            sumd =sumd +Z(j,i)*t(j)
         ENDDO
         w(i) = Sf(i)*(Fp(i) -Fc(i)) - sumd
         IF(w(i) .gt. eta * Sf(i)*(ABS(Fp(i)+ABS(Fc(i)))))THEN
            skipupdate = .false.
         ELSE
            w(i) = 0.0D0
         ENDIF
      ENDDO
      

      IF(skipupdate .EQV. .false.)THEN
         DO i =1,n
            sumd =0.0D0
            DO j = 1,n
               sumd =sumd +Z(i,j)*w(j)
            ENDDO
            t(i) = sumd
         ENDDO
         sumd =0.0D0
         DO i=1,n
            sumd = sumd + (Sx(i)*s(i))**2
         ENDDO
         denom = sumd
         DO i =1,n
            s(i) = Sx(i)**2*s(i)/denom
         ENDDO

         CALL QRUPDATE(n,t,s,1,Z,M)

         DO i =1,n
           M2(i) = M(i,i)
         ENDDO
      ENDIF


      RETURN
      END

         
