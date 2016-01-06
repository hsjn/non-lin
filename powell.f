      SUBROUTINE f_powell(n,xc,Sf,F,fp)
c------------------------------------------------------------------------
c---evaluate Jacobian at xc
c---return results in J (implemented for model test  problems only)
c------------------------------------------------------------------------

      IMPLICIT none
c  INPUT:
      INTEGER n,i,n4,i4
      REAL *8 xc(n),F(n),Sf(n),fp
C  LOCAL:


c  OUTPUT:
      F(:) = 0.0D0
      !Extended Powell singular function
      n4 = n/4

      Do i = 1,n4
         i4 = 4*i
         F(i4)   =  sqrt(10.d0)*(xc(i4-3) -xc(i4))**2
         F(i4-1) =  (xc(i4-2) -2.d0*xc(i4-1))**2
         F(i4-2) =  sqrt(5.d0)*(xc(i4-1) -xc(i4))
         F(i4-3) =  xc(i4-3)+10.d0*xc(i4-2)
      enddo
      fp = 0.0d0
      do i =1,n
         fp =fp+(sf(i)*f(i))**2
      enddo
      fp =fp*0.5d0


      RETURN
      END
      
      
