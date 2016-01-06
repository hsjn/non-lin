      subroutine nedriver(n,x0,FVEC,JAC,global,analjac,cheapF,
     .                    bw_low,bw_up,sparse_jac,
     .                    factsec,typx,typF,Fdigits,fvectol,steptol,
     .                    mintol,maxstep,itnlimit,delta,printcode,
     .                    noscale,xf,termcode)
c--------------------------------------------------------------------------
c INPUT
c n          # unknowns
c x0(1...n)  initial guess at unknowns
c FVEC       External subroutine name for evaluating F(1...n)
c            and sum of squares of residulas
C JAC        External subroutine name for evaluating Jacobian
c            (used if analjac = true. Otherwise the Jacobian
c             is calculated from finite difference or secant
c             approximations using only FVEC)
c global  =   1 ! line search 
c global  =   2 ! Hookstep trust region
c global  =   3 Dogleg trust region 
c analjac    logical true  if anlaytic Jacobian is to be used
c            analjac = true results in a call to ???
c cheapf     used only if analjac = .flase.
c              set cheapf = .true.  if f is cheap to evaluate
c             (uses finite differences if cheapf = .true.,
c              uses secant updates if cheapf = .false.)
c            if /  uses secants otherwise
c factsec    used only if analjac = .false. and cheapf = .false.
c            if factsec = .true. use factored secant approximation
c            to update a QR factorization of the Jacobian.
c            factsec = .false. is not implemented !
c typx(1...n)  gives typical magnitude of unknown i
c noscale 
c typx
c typF
c Fdigits
c fvectol
c steptol,
c mintol
c maxstep
c itnlimit    max iterations
c printcode   Fortranunit n0 for error messages
c
c OUTPUT
c xf(i)      i=1,..n solution vector
c termcode   solution indicator
c------------------------------------------------------------HSJ-----------
      USE allocate
      USE COM
      USE copy
      USE terminate
      IMPLICIT NONE
!INPUT:
      INTEGER 
     .       n,global,Fdigits,itnlimit,printcode,
     .       itncount,bw_low,bw_up
      REAL *8
     .       x0(n),typf(n),
     .       fvectol,steptol,mintol,maxstep,delta,
     .       macheps,Sx(n),eta
      LOGICAL
     .       analjac,cheapf,factsec,noscale,sparse_jac

      EXTERNAL FVEC,JAC
c LOCAL:
       REAL *8 fc,fp, FVc(n),gc(n),M2(n),M(n,n),
     .       typx(n),sumd,delta_prev,mu,phi,phip,gradmax
       REAL *8 ,dimension(:,:),allocatable :: Jc
       REAL *8 ,dimension(:)  ,allocatable :: xc,Sn,xp,Hc
       INTEGER i,j,retcode,istat,consecmax,nn
       LOGICAL 
     .       restart,maxtaken
!OUTPUT:
      INTEGER termcode
      REAL*8
     .     xf(n)
 

      

!may be allocated with different size from previous call 
      IF( allocated(Sf) )                      
     .  deallocate (Sf, STAT = istat)
        allocate (Sf(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Sf, nedriver",0,istat)


      IF(allocated(FVp))
     .  deallocate (FVp,STAT = istat)
        allocate (FVp(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("FVp, nedriver",0,istat)



      IF(allocated(Jc))
     .   deallocate (Jc,STAT = istat)
        allocate(Jc(n,n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Jc, nedriver",0,istat)

      IF(allocated(xc))
     .  deallocate (xc,STAT = istat)
        allocate (xc(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("xc, nedriver",0,istat)

      IF( allocated(xp))
     .  deallocate (xp,STAT = istat)
        allocate (xp(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("xp, nedriver",0,istat)

      IF( allocated(Sn))
     .  deallocate (Sn,STAT = istat)
        allocate (Sn(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Sn, nedriver",0,istat)

      nn = (n*(n+1))/2
      IF( allocated(Hc))
     .  deallocate (Hc,STAT = istat)
        allocate (Hc(nn),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Hc, nedriver",0,istat)




      CALL MACHINE_EPSD(macheps)
      CALL NEINCK(n,macheps,x0,typx,typf,Fdigits,fvectol,
     .                   steptol,mintol,maxstep,itnlimit,printcode,
     .                   delta, global,analjac,cheapf,factsec,
     .                   Sx,Sf,eta,termcode)


      IF(termcode .lt. 0)then
         call copya(x0,xf,n)      !copy x0 to xf
         write(printcode,'("error in input to nedriver",/, 
     .   "termcode = ",i5     )')termcode
         CALL STOP('subroutine NEDRIVER: input error', 6)
      ENDIF

      itncount = 0
      CALL nefn(n,x0,fc,FVEC)  !fc = 0.5*Norm(Sf*Fvp)**2,FVp loaded in common
      CALL nestop0(n,FVp,Sf,fvectol,termcode,consecmax)

      IF(termcode .gt. 0)then
         CALL copya(x0,xf,n)          !nestop0 returned termcode = 1
      ELSE  !get inital jacobian      !copy x0 to xf and exit below
         If(analjac) then
            CALL JAC(n,x0,Jc)
         ELSE IF(sparse_jac)then
            CALL fdjac_sparse(n,x0,FVp,FVEC,Sf,Jc,bw_low,bw_up,fvectol)
         ELSE
            CALL fdjac(n,x0,FVp,FVEC,Sx,Sf,eta,Jc)
         ENDIF
         DO i =1,n
            sumd = 0.0d0
            DO j=1,n
               sumd =sumd+Jc(j,i)*FVp(j)*Sf(j)**2
            ENDDO
            gc(i) = sumd
         ENDDO
         CALL copya(FVp,FVc,n)
         CALL copya(x0,xc,n)
      ENDIF

      restart =.true.  !used only for factored broyden method



c     iteration section follows
      DO WHILE (termcode .eq. 0)
         itncount =itncount+1
         IF(analjac .or. cheapf .or. .NOT. factsec)THEN
            CALL nemodel(n,FVc,Jc,gc,Sf,Sx,macheps,global,
     .                   M,Hc,Sn)

         ELSE
            CALL nemodelfac(n,FVc,gc,Sf,Sx,macheps,global,
     .                      restart,M,M2,Jc,Hc,Sn) 

         ENDIF
         IF(global .eq. 1)THEN
            print *,'entering linesearch'
            CALL linesearch(n,xc,fc,FVEC,gc,Sn,Sx,noscale,Sf,maxstep,
     .                      steptol,retcode,xp,fp,FVp,maxtaken)
            print *,'exit linesearch,retcode =',retcode
         ELSEIF(global .eq. 2)THEN
            print *,'entering hookdriver'
            CALL hookdriver(n,xc,fc,FVEC,gc,M,Hc,Sn,Sx,noscale,SF,
     .                    maxstep,steptol,itncount,macheps,delta,
     .                    mu,delta_prev,phi,phip,retcode,xp,fp,
     .                    maxtaken,printcode)
            call FVEC(N,XP,Sf,FVp,fp)     !this update needs to be done in
                                          !drive_dogleg
            print *,'done hookdriver,retcode =',retcode

         ELSE   !global = 3
            print *,'entering drive_dogleg'

            CALL drive_dogleg(n,xc,fc,FVEC,gc,M,Sn,Sx,noscale,Sf,
     .                        maxstep,steptol,delta,retcode,xp,
     .                        fp,maxtaken,printcode)
            call FVEC(N,XP,Sf,FVp,fp)     !this update needs to be done in
                                          !drive_dogleg
            print *,'done drive_dogleg,retcode =',retcode

         ENDIF

         print *,'retcode ,rstart,analjac,cheapftermcode=',
     .   retcode,restart,analjac,cheapf,termcode


         !retcode = 1 means satisfactory new step not found
         !otherwise retcode = 0
 10      IF(retcode .eq. 1 .and. (.not. restart) .and. (.not. analjac)
     .         .and. (.not. cheapf))THEN 
             !this section is for factsec = true cases 
             !( .not. analjac and .not. cheapf ==> factsec = true)
             IF(sparse_jac)then
                CALL fdjac_sparse(n,xc,FVc,FVEC,Sf,Jc,bw_low,
     .                            bw_up,fvectol)
             ELSE
                CALL fdjac(n,xc,FVc,FVEC,Sx,Sf,eta,Jc)
             ENDIF
            DO i =1,n
               sumd = 0.0d0
               DO j=1,n
                  sumd =sumd+Jc(j,i)*FVc(j)*Sf(j)**2
               ENDDO
               gc(i) = sumd
            ENDDO
            IF(global .eq. 2 .or. global .eq. 3)delta = -1
            restart = .true.
            termcode = 0    !HSJ added
         ELSE   !complete the iteration 
            IF(analjac)then
               CALL JAC(n,xp,Jc)
            ELSEIF(cheapf)then
               IF(sparse_jac)then
                   CALL fdjac_sparse(n,xp,FVp,FVEC,Sf,Jc,bw_low,
     .                            bw_up,fvectol)
               ELSE
                   CALL fdjac(n,xp,FVp,FVEC,Sx,Sf,eta,Jc)
               ENDIF
            ELSEIF(factsec)then
               print *,'calling broyfac'
               CALL broyfac(n,xc,xp,FVc,Fvp,eta,Sx,Sf,Jc,
     .                      M,M2)
            ELSE !factsec =.false.
               CALL broyunfac(n,xc,xp,FVc,Fvp,eta,Sx,Jc)
            ENDIF
            IF(factsec)THEN      !get gc ( gradient) from QR factorization
               !gc = Jc*Sf*FVp     !note that Jc contains Q transpose here
               Do i =1,n  
                  sumd =0.0d0
                  DO j=1,n
                     sumd =sumd+Jc(i,j)*FVp(j)*Sf(j)
                  ENDDO
                  gc(i) = sumd
               ENDDO
               Do i =n,1,-1
                  sumd =0.0d0
                  DO j=1,i
                     sumd =sumd+M(j,i)*gc(j)
                  ENDDO
                  gc(i) = sumd
               ENDDO
            ELSE                 !get gc ( gradient) from definition
               DO i =1,n
                  sumd = 0.0d0
               DO j=1,n
                  sumd =sumd+Jc(j,i)*FVp(j)*Sf(j)**2
               ENDDO
                  gc(i) = sumd
               ENDDO
            ENDIF
            CALL check_convergence(n,xc,xp,FVp,fp,gc,Sx,Sf,retcode,
     .                  fvectol,steptol,noscale,itncount,itnlimit,
     .                  maxtaken,consecmax,termcode,printcode,gradmax)
            print *,'termcode afte check_conver',termcode
            IF(termcode .eq. 2 .and. (.not. restart) .and. 
     .                (.not. analjac) .and. (.not. cheapf))THEN
               !restart
               retcode = 1
               print *,'restart,retcode ,termcode =',
     .                  restart,retcode,termcode
               go to 10
            ELSEIF (termcode .gt. 0)THEN
               call copya(xp,xf,n)
            ELSE
               restart = .FALSE.
            ENDIF
            CALL copya(xp,xc,n)
            fc = fp
            FVc(:) = FVp(:)
         ENDIF
      ENDDO  ! termcode = 0 DO WHILE




      RETURN
      END

               
