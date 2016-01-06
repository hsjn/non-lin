      SUBROUTINE nestop0(n,F,Sf,fvectol,termcode,consecmax)
c-----------------------------------------------------------------------
c---- determine if we should stop at iteration 0
c---- because the initial guess is good enough.
c---- return termcode =0 if the guess x0 is not sufficient
c---- return termcode =1 if x0 is a sufficient approximation to the
c---- root of F(x0) =0.0
c---- condition is that maximum component is .01*fvectol or less

c--------------------------------------------------------------HSJ------

      IMPLICIT none
c  INPUT:
      INTEGER
     .        n
      REAL*8
     .        f(n),Sf(n),fvectol
c  LOCAL
      INTEGER
     .        i
      REAL*8
     .         mmax,factor
c  OUTPUT:
      INTEGER
     .        termcode,consecmax

      
      factor = 1.D-2
      consecmax = 0
      mmax =0.0d0
      DO i = 1,n
         mmax = MAX(mmax,Sf(i)*ABS(F(i)))
      ENDDO
      IF(mmax .lt. factor*fvectol)THEN
         termcode = 1
      ELSE
         termcode =0
      ENDIF


      RETURN
      END
