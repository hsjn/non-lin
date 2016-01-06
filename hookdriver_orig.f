                                                               
                                                                      
                                                                        
      subroutine hookdriver(n,xcur,fcur,FVEC,gradient,L,hessian,       
     &                 step,scale,noscale,Sf,steplength_max,step_tol,
     &                 iterat,epsilond,delta,mu,deltaprev,        
     &                 phi,phip,retcode,xnew,fnew,max_taken,iounit)      
!---------------------------------------------------------------------  
C INPUT 
c n
c xcur
c fcur
C FVEC
c gradient
c L(n,n)                



!---the following vectors are used for local storage
!  s(i)          i=1,2..n
!  xprev(i)      i=1,2..n  (f90 automatic arrays here)
!--------------------------------------------------------------------hsj
      USE copy

      implicit none
                              ! in: number of data points
      integer n,iterat,retcode,j,                                       &
     &      i,k,step_type,iounit
      real *8 xcur(*),L(n,n),hessian(*),gradient(*),step(*),xnew(*),      &
     &      s(n),xprev(n),temp,beta,epsilond,fnew,fcur,fprev
      real*8 scale(*),steplength_max,step_tol,delta,                    &
     &    deltaprev,phi,phip,alpha,sqrt,steplength,                     &
     &    dumy,phipinit,Sf(n)
      real*8 mu
      logical max_taken,noscale,firsthook,newton_step
      ! need a lower limit on steplength; this is arbitrary but a lot
      ! bigger than zero.                                               
      real*8, parameter :: steplength_min=1.0e-5
                                                                        
      external FVEC
                                   ! exception handling routine

                                                     
      retcode=4
                  
      step_type=1  !selects hookstep in trust region
      firsthook=.true.
      call euclidnorm(step,scale,n,noscale,steplength)
                                                                        
      ! purpose of this section is to calculate a new delta
      ! upon initialization                                             
      if(iterat .eq. 1 .or. delta .eq. -1.)then
          mu=0.0
          if(delta .eq. -1)then
              if(.not. noscale)then
                  do j=1,n
                      s(j)=gradient(j)/scale(j)
                  enddo
              else
                  call dcopy(n,gradient,1,s,1)
              endif
              call euclidnorm(s,scale,n,.true.,alpha)
              ! eliminate arithmetic problems from very large
              ! (and certainly pathological) alpha                      
              alpha = min(alpha, 1.0d8)
              alpha=alpha*alpha
              beta=0.0d0
              do i=1,n
                  temp=0.0d0
                  do j=i,n
                      dumy=1.0
                      if(.not. noscale)dumy=scale(j)
                      temp=temp+L(j,i)*gradient(j)/(dumy*dumy)
                  enddo
                  beta=beta+temp*temp
              enddo
              delta=alpha*sqrt(alpha)/beta
              if(delta .gt. steplength_max)delta=steplength_max
              ! if we came in with a gradient effectively at the
              ! solution, we can have an extremely small delta;         
              ! zero is even a possibility.  this will cause troubles
              ! in hookstep.  defend against that.
              if (delta .lt. steplength_min) delta=steplength_min
          endif
      endif
!                                                                       
!                                                                       
      do while (retcode .ge. 2)
          call hookstep(n,gradient,l,hessian,step,scale,noscale,        
     &              epsilond,delta,mu,deltaprev,phi,phip,firsthook,     
     &              phipinit,s,newton_step,steplength,                  
     &              iounit)
          deltaprev=delta
          call trust_region(n,xcur,fcur,FVEC,gradient,l,s,scale,        
     &        noscale,Sf,newton_step,steplength_max,step_tol,step_type,
     &        hessian,delta,retcode,xprev,fprev,xnew,fnew,max_taken)
                                    
c         if(phip .eq. 0.0)retcode=1    !hsj 5/10/93 accept solution
         print *,'phip,retcode =',phip,retcode
      enddo
!                                                                       
                                                  
      return
      end 
