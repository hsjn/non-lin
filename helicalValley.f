      SUBROUTINE f_helical(n,xc,Sf,F,fp)
c------------------------------------------------------------------------
c---evaluate Jacobian at xc
c---return results in J (implemented for model test  problems only)
c------------------------------------------------------------------------

      IMPLICIT none
c  INPUT:
      INTEGER n,i
      REAL *8 xc(n),F(n),Sf(n),fp,theta
C  LOCAL:


c  OUTPUT:
      F(:) = 0.0D0
      F(1) = 10.d0*(xc(3)-10.d0*theta(xc(1),xc(2)))
      F(2)= 10.d0*(SQRT(xc(1)**2+xc(2)**2)-1.d0)
      F(3) = xc(3)
      fp = 0.0d0
      do i =1,n
         fp =fp+(Sf(i)*F(i))**2
      enddo
      fp =fp*0.5d0


      RETURN
      END
      
      
      real *8 function theta(x1,x2)
      
c-------------------------------------------
      IMPLICIT none
      REAL *8 x1,x2,tpi
      tpi =6.283185308D0
      if(x1 .gt. 0.0d0)then
         theta =(1.d0/tpi)*ATAN(x2/x1)
      else
         theta =(1.d0/tpi)*ATAN(x2/x1)+0.5
      endif


      return
      end
