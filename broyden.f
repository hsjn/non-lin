      SUBROUTINE f_broyden(n,xc,Sf,F,fp)
c------------------------------------------------------------------------
c---evaluate Jacobian at xc
c---return results in J (implemented for model test  problems only)
c------------------------------------------------------------------------

      IMPLICIT none
c  INPUT:
      INTEGER n,i
      REAL *8 xc(n),F(n),Sf(n),fp
C  LOCAL:


c  OUTPUT:

      F(1) = xc(1) + xc(2) -3.d0
      F(2)=  xc(1)**2 + xc(2)**2 -9.D0

      fp = 0.0d0
      do i =1,n
         fp =fp+(Sf(i)*F(i))**2
      enddo
      fp =fp*0.5d0


      RETURN
      END
      
      
