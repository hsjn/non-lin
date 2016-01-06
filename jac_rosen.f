      SUBROUTINE jac_rosen(n,xc,JACR)
c------------------------------------------------------------------------
c---evaluate Jacobian at xc
c---return results in J (implemented for model test  problems only)
c------------------------------------------------------------------------
      USE terminate
      IMPLICIT none
c  INPUT:
      INTEGER n,i,j,k
      REAL *8 xc(n)
C  LOCAL:


c  OUTPUT:
      REAL*8 JACR(n,n)

      JACR(:,:) = 0.0D0
      !Rosenbrock extended function
         DO i =1,n/2
            j= 2*i
            k = j-1
            JACR(k,k) = -20.D0*xc(k)
            JACR(k,j) = 10.0D0
            JACR(j,k) = -1.D0
            JACR(j,j) = 0.D0
         ENDDO

      RETURN
      END
