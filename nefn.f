      SUBROUTINE nefn(n,xp,fp,FVEC)
c------------------------------------------------------------------------
c     get 0.5* sum of squares of residuals of nonlinear eautions
c
c INPUT:
c   n
c   xp(n)
c
c
c OUTPUT:
c
c fp         = 0.5*Norm(Sf*fVp)**2


c-------------------------------------------------------------HSJ-----------
      USE COM     !Sf(n),FVp(n) brought in here
      IMPLICIT NONE

      INTEGER
     .        n,i
      REAL*8 xp(n),fp

      EXTERNAL FVEC



c     evaluate the set of equations at xp:
      CALL FVEC(n,xp,Sf,FVp,fp)
      print *,'Fvp(1..N) in nefn'
      do i =1,n
         print *,'i,xp(i),Sf(i),FVp(i)',i,xp(i),Sf(i),FVp(i)
      enddo
      print *,'fp in nefn =',fp
      RETURN
      END
