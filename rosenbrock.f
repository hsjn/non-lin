      SUBROUTINE f_rosen(n,xc,Sf,F,fp)
c------------------------------------------------------------------------
c---evaluate Jacobian at xc
c---return results in J (implemented for model test  problems only)
c------------------------------------------------------------------------
      IMPLICIT none
c  INPUT:
      INTEGER n,i,j,k
      REAL *8 xc(n),F(n),Sf(n),fp
C  LOCAL:


c  OUTPUT:

      F(:) = 0.0D0
      !Rosenbrock extended function
         DO i =1,n/2
            j= 2*i
            k = j-1
            F(k) = 10.0d0*(xc(j)-xc(k)**2)
            F(j) = 1.d0-xc(k)
         ENDDO
         fp =0.0D0
         do i =1 ,n
             fp =fp + (Sf(i)*F(i))**2
         enddo
         fp = fp*0.5D0


      RETURN
      END
