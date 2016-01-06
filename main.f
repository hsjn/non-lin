      program non_lin_solve
c     example use of non linear equation solver system
c-------------------------------------------------------------HSJ
c
       USE allocate
       USE terminate
       IMPLICIT none
       INTEGER n,nmax,i,j,k,l,istat,global,itnlimit,n2,
     .         printcode,termcode,Fdigits,n4,nmin,intv,
     .         bw_low,bw_up
       REAL *8 factor, fvectol, steptol,mintol,
     .         maxstep, delta
       REAL *8 ,dimension(:)  ,allocatable :: x0,typx,typF,xf
       LOGICAL  noscale,analjac,cheapf,factsec,sparse_jac
       CHaracter *24 FUNCTION
       EXTERNAL f_rosen, jac_rosen,f_powell,jac_powell,
     .          f_helical,jac_helical,f_broyden,jac_broyden
      REAL *8 ,dimension(:)  ,allocatable :: xcRandom

        FUNCTION = "Powell"
        !FUNCTION = "Broyden"
        !FUNCTION = "Rosenbrock"
        !FUNCTION = "Helical"
       printcode = 6  !io for messages
       noscale = .false.
      
      nmax =800
      !nmax=2
       !nmax = 3
       global = 1                   ! line search 
        !global  = 2                  ! Hookstep trust region
       !global  = 3                  !Dogleg trust region 
      analjac = .true.
      analjac = .false.
      cheapf = .true.                 !enables broyfac if false 
      factsec = .false.  
      !factsec = .true.        !used if analjac = .false. and cheapf = .false.
      sparse_jac = .true.
      if( .not. analjac .and. .not. cheapf)factsec = .true.            
      fdigits = 16            !could try -1 here
      fvectol = 1.d-18
      steptol = 1.d-16
      mintol  = 1.d-14
      maxstep = 1.d0
      itnlimit = 5000
      delta = -1.0
      intv = 2
      IF(FUNCTION .eq. 'Helical')then
         intv = 1
         nmin =nmax
         n2 = nmax
      ELSE
         intv = 2
         nmin = nmax
         n2 = nmax/2
      ENDIF

      DO n=nmin,nmax,intv
         IF( allocated(xcRandom) )
     .          deallocate (xcRandom, STAT = istat)
            allocate ( xcRandom(n),STAT = istat)

         IF( allocated(x0) )
     .          deallocate (x0, STAT = istat)
            allocate ( x0(n),STAT = istat)
         IF(istat .ne. 0)
     .          call allocate_error("xc, main",0,istat)
         IF( allocated(xf) )
     .          deallocate (xf, STAT = istat)
            allocate ( xf(n),STAT = istat)
         IF(istat .ne. 0)
     .          call allocate_error("xf, main ",0,istat)
         IF( allocated(typx) )
     .          deallocate (typx, STAT = istat)
            allocate ( typx(n),STAT = istat)
         IF(istat .ne. 0)
     .          call allocate_error("typx, main ",0,istat)
         IF( allocated(typF) )
     .          deallocate (typF, STAT = istat)
           allocate ( typF(n),STAT = istat)
         IF(istat .ne. 0)
     .          call allocate_error("typF, main ",0,istat)

         factor = 1.0D0

         CALL Random_Number(xcRandom)         !uses same seed for repeat cases
         DO l = 1,1



            SELECT CASE(FUNCTION)



            CASE ( "Rosenbrock")





            factor = factor *10.0D0**(l-1)

         if(2*n2 .ne. n) 
     .       CALL STOP('n must be divisible by 2 for Rosenbrock',
     .                  printcode)
         

            DO i = 1,n2
               j=2*i
               k= j-1
               x0(k) = DBLE(xcRandom(k)) *factor
               x0(j) = DBLE(xcRandom(j)) *factor
               typx(k) = 1.0
               typx(j) = 1.0
               typF(k) =1.0
               typF(j) = 1.0
            ENDDO
            bw_low = 1
            bw_up  = 1
             call nedriver(n,x0,f_rosen,jac_rosen,global,analjac,
     .                     cheapF,bw_low,bw_up,sparse_jac,
     .                    factsec,typx,typF,Fdigits,fvectol,steptol,
     .                    mintol,maxstep,itnlimit,delta,printcode,
     .                    noscale,xf,termcode)
            print *,'SOLUTION For Rosenbrock function IS'
            DO i = 1,n
              print *,'xf(i) = ',xf(i)
            ENDDO




         CASE ( "Powell")


         !extenden Powel Singular Function
         !n must be a multiple of 4
         n4=n/4
         if(4*n4 .ne. n) 
     .       CALL STOP('n must be divisible by 4 for Powell',
     .                  printcode)

            factor = factor *10.0D0**(l-1)
            DO i = 1,n2
               j=2*i
               k= j-1
               x0(k) = xcRandom(k) *factor
               x0(j) = -xcRandom(j) *factor
               typx(k) = 1.0
               typx(j) = 1.0
               typF(k) =1.0
               typF(j) = 1.0
            ENDDO

             bw_low = 3
             bw_up = 2

             call nedriver(n,x0,f_powell,jac_powell,global,analjac,
     .                     cheapF,bw_low,bw_up,sparse_jac,
     .                    factsec,typx,typF,Fdigits,fvectol,steptol,
     .                    mintol,maxstep,itnlimit,delta,printcode,
     .                    noscale,xf,termcode)
           print *,'SOLUTION For Powell function IS'
           DO i = 1,n
              print *,'xf(i) = ',xf(i)
           ENDDO


         CASE ( "Helical")


         !Hellical valley  Function
         !n must be 3
         if(n  .ne. 3) 
     .       CALL STOP('n must be 3 for hellical',printcode)

            DO i = 1,n
               x0(i) =  15.0D0 *factor
               typx(i) = 1.0D0
               typF(i) = 1.0D0
            ENDDO
              bw_low = 1
              bw_up = 2
             call nedriver(n,x0,f_helical,jac_helical,global,analjac,
     .                     cheapF,bw_low,bw_up,sparse_jac,
     .                    factsec,typx,typF,Fdigits,fvectol,steptol,
     .                    mintol,maxstep,itnlimit,delta,printcode,
     .                    noscale,xf,termcode)
           print *,'SOLUTION For Helical Valley function IS'
           DO i = 1,n
              print *,'xf(i) = ',xf(i)
           ENDDO

         CASE ( "Broyden")


         !boyden update test case
         !n must be 2
         if(n  .ne. 2) 
     .       CALL STOP('n must be 2 for broyden',printcode)

               x0(1) =  1.0D0
               x0(2) =  5.0D0
               typx(:) = 1.0D0
               typF(:) = 1.0D0
              bw_low = 1
              bw_up = 1
             call nedriver(n,x0,f_broyden,jac_broyden,global,analjac,
     .                     cheapF,bw_low,bw_up,sparse_jac,
     .                    factsec,typx,typF,Fdigits,fvectol,steptol,
     .                    mintol,maxstep,itnlimit,delta,printcode,
     .                    noscale,xf,termcode)
           print *,'SOLUTION For Broyden function IS'
           DO i = 1,n
              print *,'xf(i) = ',xf(i)
           ENDDO

           END SELECT



         ENDDO ! end factor (l) loop
         
      ENDDO 


      STOP
      END
