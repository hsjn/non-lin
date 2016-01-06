      SUBROUTINE neinck(n,macheps,x0,typx,typf,Fdigits,fvectol,
     .                   steptol,mintol,maxstep,itnlimit,printcode,
     .                   delta, global,analjac,cheapf,factsec,
     .                  Sx,Sf,eta,termcode)
c--------------------------------------------------------------------

c----- check input
c
c-------------------------------------------------------------HSJ---





      IMPLICIT NONE
c INPUT:
      INTEGER
     .        n,Fdigits,itnlimit,printcode,global
      REAL *8
     .        x0(n),typx(n),typf(n),delta,
     .         macheps, fvectol,steptol,mintol,maxstep
      LOGICAL
     .        analjac,cheapf,factsec
c LOCAL:
      INTEGER i
      
      REAL*8
     .       a

c OUTPUT
      INTEGER
     .      termcode
      REAL *8
     .      Sx(n),Sf(n),eta

c-------------------------------------------------------------------


      IF(n .lt. 1)THEN
         termcode = -1
         RETURN
      ENDIF
      DO i=1,n
         Sx(i) = 1./typx(i)
      ENDDO
      DO i=1,n
         Sf(i) = 1./typf(i)
      ENDDO
      IF(Fdigits .eq. -1)THEN
         eta = macheps
      ELSE
         a = 10.d0**(-Fdigits)
         eta = MAX(macheps,a)
      ENDIF
      if(global .eq. 2 .or. global .eq. 3 .and. delta .eq. 0.0D0)
     .         delta = -1.D0

      
      RETURN
      END
