      SUBROUTINE broyunfac(n,xc,xp,FVc,Fvp,eta,Sx,Jc)
c---------------------------------------------------------------
c --- unfactored Broyden method,not verys useful ??
c -- not implemented at this time 
c-----------------------------------------------------------HSJ
      USE terminate
      IMPLICIT NONE
!INPUT:
      INTEGER 
     .       n
      REAL *8
     .       xc(n),xp(n),Jc(n,n),FVc(n),FVp(n),
     .       Sx(n),eta
      CALL STOP("SUB BROYUNFAC NOT IMPLEMENTED",6)
      RETURN
      END
