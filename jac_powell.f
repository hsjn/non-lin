      SUBROUTINE jac_powell(n,xc,JACP)
c------------------------------------------------------------------------
c---evaluate Jacobian at xc
c---return results in J (implemented for model test  problems only)
c------------------------------------------------------------------------
      USE terminate
      IMPLICIT none
c  INPUT:
      INTEGER n
      REAL *8 xc(n)
C  LOCAL:
      INTEGER i,n4,i4

c  OUTPUT:
      REAL*8 JACP(n,n)

      JACP(:,:) = 0.0D0
      !powell extended function
      n4 =n/4
         DO i =1,n4
            i4 = 4*i
            JACP(i4-3,i4-3) = 1.D0
            JACP(i4-3,i4-2) = 10.d0
            JACP(i4-2,i4-1) = sqrt(5.d0)
            JACP(i4-2,i4)   = -sqrt(5.d0)
            JACP(i4-1,i4-2) =2.d0*(xc(i4-2)-2.d0*xc(i4-1))
            JACP(i4-1,i4-1) =-4.d0*(xc(i4-2)-2d0*xc(i4-1))
            JACP(i4,i4-3) = 2.d0*sqrt(10.d0)*(xc(i4-3) -xc(i4))
            JACP(i4,i4) = -2.d0*sqrt(10.d0)*(xc(i4-3) -xc(i4))
         ENDDO
      RETURN
      END
