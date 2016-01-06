      SUBROUTINE fdjac(n,xc,fc,FVEC,Sx,Sf,eta,jacf)
c------------------------------------------------------------------------
c--- evaluate finite difference Jacobian


c------------------------------------------------------------------------



      IMPLICIT none
c  INPUT:
      INTEGER n
      REAL *8 xc(n),fc(n),Sx(n),Sf(n),eta
C  LOCAL:
      INTEGER
     .        j,i

      REAL*8
     .       sqrteta,stepsizej,tempj,Fj(n),ssqrc

      LOGICAL test
c  OUTPUT:
      REAL*8 jacf(n,n)

      EXTERNAL FVEC


      
      sqrteta =SQRT(eta)
      DO j = 1,n
         !calculate column j of JACF
         stepsizej = sqrteta*MAX(ABS(xc(j)),1./Sx(j))*SIGN(1.D0,xc(j))
         tempj = xc(j)
         xc(j) = xc(j) + stepsizej
         stepsizej = xc(j)-tempj   !reduces finite precision error slightly 
         CALL FVEC(n,xc,Sf,Fj,ssqrc)
         
         DO i = 1, n
            jacf(i,j) = (Fj(i) -Fc(i))/stepsizej
         ENDDO
         xc(j) = tempj     !set xc back to original value
      ENDDO

      RETURN
      END
