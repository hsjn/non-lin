       SUBROUTINE nemodel(n,fc,jacf,g,Sf,Sx,macheps,global,
     .                  M,H,Sn)                               
c------------------------------------------------------------------------
c--- solve J*dx = -F(xc)
c    Factor model Jacobian,calculate Newton step
c    If Jacobian is singular modify it 
c---
c INPUT
c n              # unknowns
c fc(1..n)       Value of n dimensional function Fc(xc)  at xc
c jacf(n,n)      Jacobian of Fc
c g(1..n)        
c Sf(1..n)
c Sx(1..n)
c macheps
c global         solution method 

c global         = 1    line search 
c global         = 2    Hookstep trust region
c global         = 3    Dogleg trust region 

c OUPUT:
c M(n,n)        M is put into LOWER TRIANGULAR SYMMETRIC STORAGE
c               MODE on output.
c               (IE M(I,J),FOR J  .LE. I,IS ACCESSED AS M(K),WHERE
c               K= I*(I-1)/2 +J, for i .ge. j
c
c H             H IS IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE.
c               (IE H(I,J),FOR J  .GE. I,IS ACCESSED AS H(K),WHERE
c               K=((I-1)(2*N-I)+2*J)/2
c Sn(n)
c-----------------------------------------------------------HSJ-------------

       USE terminate

      IMPLICIT none
c  INPUT:
      INTEGER n,global
      REAL *8 fc(n),jacf(n,n),g(n),Sf(n),Sx(n),macheps
C  LOCAL:
      INTEGER
     .        i,j,k

      REAL*8
     .       M1(n),M2(n),est,sumd,Hnorm,temp,maxadd,zero
      LOGICAL sing
c  OUTPUT:
      REAL*8 M(n,n),H(*),Sn(n) !H has  n(n+1)/2 elements


c TEMPORARY
      REAL *8 work(n)
      INTEGER info 
      zero = 0.0d0


      DO i = 1,n
         DO j = 1,n
            M(i,j) = Sf(i)*jacf(i,j)
         ENDDO
      ENDDO


c------------------------------------------------------------
c      qrdecomp,coded from the Dennis- Schnabel algorithm
c      appears to be faulty in some way. So replace it with
c      Linpack routine dgeqrf and compensate for different
c      storage scheme:
c      CALL qrdecomp(n,M,M1,M2,sing)

      sing = .false.
      Call Dgeqrf(n,n,M,n,M1,work,n,info) ! lapack routine
c     returns qr factorization of jacobian in M, presumabley, if 
c     a diagonal element of r is zero then we dont crash
c     

      do i = 1,n
         M2(i) = M(i,i)  !M1 is ony used in qrsolve
         if(M2(i) .eq. 0.0d0)sing = .true.
      enddo              !which is also replaced so M1 becomes
                         !unused

      IF( .not. sing)THEN
         DO j =1 ,n
            DO i = 1,j-1
               M(i,j) = M(i,j)/Sx(j)
            ENDDO
            M2(j) = M2(j)/Sx(j)
         ENDDO
         CALL condest(n,M,M2,est)
      ELSE
         est = zero
      ENDIF  ! sing 
       print *,'sing,est  =',sing,est 
      IF(sing .or. est .gt. 1.0d0/SQRT(macheps))THEN
         print *,'JACOBIAN SINGULAR************************'
         DO i =1,n   !perturb jacobian so it isnt singular
            DO j =i,n
               sumd = 0.0d0
               DO k =1,n
                  sumd = sumd + jacf(k,i)*jacf(k,j)*Sf(k)**2
               ENDDO
               k=((i-1)*(2*n-i)+2*j)/2  !upper triangular storage for H
               H(k) = sumd 
            ENDDO
         ENDDO
         sumd =0.0d0
         DO j =1,n
            sumd =sumd + ABS(H(j))/Sx(j)
         ENDDO
         Hnorm = sumd/Sx(1)
         DO i = 2,n
            sumd =0.0D0
            DO  j= 1,i
                K=((j-1)*(2*N-j)+2*i)/2
                sumd = sumd + ABS(H(k))/Sx(j)
            ENDDO
            temp = sumd
            sumd =0.0D0
            DO j = i+1,n
                k=((i-1)*(2*n-i)+2*j)/2
                sumd = sumd + ABS(H(k))/Sx(j)
            ENDDO
            temp = temp +sumd
            temp = temp/Sx(i)
            Hnorm = MAX(Hnorm,temp)
         ENDDO
         DO i = 1,n
            k=((i-1)*(2*n-i)+2*i)/2
            H(k) = H(k)+ SQRT(n*macheps)*Hnorm*Sx(i)**2
         ENDDO
       
         
         CALL choldecomp(n,H,zero,macheps,M,maxadd)
         CALL cholsolve(n,g,M,Sn)


      ELSE  ! jacobian not singular, get normal Newton step:
         DO j = 1,n
            DO i =1 , j-1
               M(i,j) = M(i,j)*Sx(j)
            ENDDO
            M2(j) = M2(j)*Sx(j)
         ENDDO
         DO j= 1,n
            Sn(j) = -Sf(j)*Fc(j)      
         ENDDO


c-----------------------------------------------------------------------
c        now solve J*dx = -F(xc), where dx is the Newton step increment
c        (dx = xc new - xc prev). J is in factored form from above:
c     dont use qrsolve because qrdecomp doesn't produce the
c     right factors! (see above):
c         Call qrsolve(n,M,M1,M2,Sn)
c     qrsolve doesnt change M,M1,M2 . output is Sn (eg the Newton step)
c
c      REPLACEMENT for qrsolve:
c     Recall that in call to Dgeqrf we ended up with the matrix R 
c     in the upper triangle of M.
c      To solve the system
c           J*dx = Sn  (see Sn above,defined with minus sign included)
c     we have J = Q_1 * R   (Q = (Q_1,Q-2) , Q_1,Q_2 orthogonal, but
c     Q_2 drops out because it multiplies the lower half of the factored
c     matrix which is zero)
c      hence Q_1*R*dx  = Sn is just R*dx = Q_1(transpose)*Sn
c     this upper triangular system is easily solved.
C     The reason why all this is done is to control what happens
c     when J is near singular. 
c     First get Q_1 Transpose *Sn:


         CALL DORMQR('L','T',n,1,n,M,n,M1,Sn,n,work,n,info)



c     Now solve R * dx = Sn:

         CALL  DTRTRS('U','N','N',n,1,M,n,Sn,n,info)
c        if info = i>0  then element r(i,i) is zero
c        we checked for this above after call to Dgeqrf
c        so this situation shouldnt arise here:
         if(info .gt. 0)
     .   CALL STOP('subroutine NEMODEL: R singular', 0)

c---------------------------------------------------------------------


         IF(global .eq. 2 .or. global .eq. 3) THEN
            !put R transpose into lower triangle of M:
            DO i =1,n
               M(i,i) = M2(i)
               Do j= 1,i-1
                  M(i,j)=M(j,i)
               ENDDO
            ENDDO
         ENDIF     
         IF(global .eq. 2) THEN
            !store  R transpose* R in upper triangular symmetric storage
            !mode in H: 
               DO i =1,n                       
               sumd =0.0d0
               DO k =1,i
                  sumd =sumd + M(i,k)**2
               ENDDO
               K=((i-1)*(2*n-i)+2*i)/2         !picks out diagonal elements
               H(k) = sumd
               DO j = i+1,n
                  sumd =0.0d0
                  DO  k = 1,i
                     sumd =sumd + M(i,k)*M(j,k)
                  ENDDO
                  K=((I-1)*(2*N-I)+2*J)/2
                  H(k) = sumd
               ENDDO
             ENDDO 
          ENDIF
         
         
      ENDIF


      RETURN
      END


