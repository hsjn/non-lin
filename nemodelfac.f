      SUBROUTINE nemodelfac(n,fc,g,Sf,Sx,macheps,global,
     .                  restart,M,M2,jacf,H,Sn)
c------------------------------------------------------------------------
c--- Factor model Jacobian,calculate Newton step,If Jacobian is singular
c--- modify it
c INPUT
c n
c fc(n)
c g(n)
c Sf(n)
c Sx(n)
c macheps
c global
c restart

c INPUT/OUTPUT
c M(n,n) M2(n),jacf(n,n)
c              If this is a restart then M contains no useful
c              information on input
c              otherwise M contains R in its upper triangle, including the
c              diagonal elements.
c              On output M contains R in its upper triangle and
c              diagonal elements are in M2.

c OUTPUT:
c H(n)              upper triangular storage mode
c Sn(n)

c------------------------------------------------------------------------



      IMPLICIT none
c  INPUT:
      INTEGER n,global
      REAL *8 fc(n),g(n),Sf(n),Sx(n),macheps
      LOGICAL restart
C  LOCAL:
      INTEGER
     .        i,j,k,info

      REAL*8
     .       M1(n),work(n),est,sumd,Hnorm,temp,maxadd,zero
      LOGICAL sing
c  OUTPUT:
      REAL*8 M(n,n),M2(n),jacf(n,n),H(*),Sn(n)

c TEMPORAY DEBUG
         REAL*8 jac_copy(n,n),jactr(n,n)



       print *,'in nemodelfac'





         jac_copy(:,:) = jacf(:,:)

      zero = 0.0D0


      IF(restart) THEN                   !get new QR factorization
         !we start with a two dimensional M:

         DO i =1,n
            DO j=1,n
               M(i,j)= Sf(i)*jacf(i,j)
            ENDDO
         ENDDO

c        CALL qrdecomp(n,M,M1,M2,sing)
c        M is in factored form (qform does not modify M)
c        R is in upper triangle of M but diagonal elements
c        of R are in M2. Lower triangle of M contains
c        information to generate Q.
c        next factored form of M is used to get Q transpose
c        returned in jacf:
c        Call qform(n,M,M1,jacf)
c                   HOWEVER
c        qrdecomp,coded from the Dennis- Schnabel algorithm
c        appears to be faulty in some way. So replace it with
c        Linpack routine dgeqrf and compensate for different
c        storage scheme:

        sing = .false.  !we have to assume this here Dgeqr doesnt supply info
        jacf(:,:) = M(:,:)
        Call Dgeqrf(n,n,jacf,n,M1,work,n,info)      ! lapack routines
c       jacf now contains R in its upper triangle,including the diagonal
c       Put info into M,m2 as would have been done by qrdecomp:

        do i = 1,n
           M2(i) = jacf(i,i) 
        enddo   
        M(:,:) = jacf(:,:)              !puts M into factored form
                                        !as would have been done by qrdecomp
        Call Dorgqr(n,n,n,jacf,n,M1,work,n,info)    !jacf now contains Q transpose
        jactr(:,:)= jacf(:,:)
        do i=1,n
           do j=1,n
              jacf(i,j) =jactr(j,i)
           enddo
        enddo

c        Call Dorgql(n,n,n-1,jacf,n,M1,work,n,info)    !jacf now contains Q transpos



      ELSE  !restart = .false.
         sing = .false.
         DO i =1,n
            if(M(i,i) .eq. 0.0d0 )sing  = .TRUE.
         ENDDO
      ENDIF
      

      IF(.not. sing) THEN
         DO j = 1,n
            DO i = 1,j-1
                M(i,j) = M(i,j)/Sx(j)
            ENDDO
            M2(j) = M2(j)/Sx(j)
         ENDDO
         CALL condest(n,M,M2,est)
         DO j = 1,n
            DO i = 1,j-1
                M(i,j) = M(i,j)*Sx(j)
            ENDDO
            M2(j) = M2(j)*Sx(j)
         ENDDO
      ELSE
         est =0.0d0
      ENDIF


      If(sing .or. (est .gt. 1.D0/SQRT(macheps)))THEN
        !perturb Jacobian H = RT*R
         DO i =1,n
            DO j =i,n
               sumd =0.0D0
               DO k =1,i
                  sumd =sumd +M(k,i)*M(k,j)
               ENDDO
               k = ((i-1)*(2*n-i)+2*j)/2  !upper triangular storage for H
               H(k) = sumd 
            ENDDO
         ENDDO
         sumd =0.0D0
         DO j =1,n
c            sumd =sumd +ABS(H(1,j))/Sx(j)
            sumd =sumd +ABS(H(j))/Sx(j)
         ENDDO
         Hnorm = sumd/Sx(1)
         DO i =2,n
            sumd =0.0D0
            DO j =1,i
               k = ((j-1)*(2*n-j)+2*i)/2
               sumd =sumd +ABS(H(k))/Sx(j)
            ENDDO
            temp =sumd
            sumd =0.0d0
            DO j =i+1,n
               k = ((i-1)*(2*n-i)+2*j)/2
               sumd =sumd+ABS(H(k))/Sx(j)
            ENDDO
            temp = temp +sumd
            temp=temp/Sx(i)
            Hnorm = MAX(temp,Hnorm)
         ENDDO
         DO i =1,n
            k = ((i-1)*(2*n-i)+2*i)/2
            H(k) = H(k)+SQRT(n*macheps)*Hnorm*Sx(i)**2
         ENDDO


         !choldecomp takes H and outputs M (as a lower
         !triangular matrix)
         CALL choldecomp(n,H,zero,macheps,M,maxadd)
         !cholsolve takes lower triangular M and 
         !solves for Sn, the "newton' step. g is not changed:
         CALL cholsolve(n,g,M,Sn)
      ELSE
         !do normal Newton step
         !note that jacf  = Q transpose here
         !note that M contains R in upper triangle except  diagonal
         !which  is in M2
         DO i =1,n
            sumd =0.0D0
            DO j =1,n
               sumd =sumd+jacf(i,j)*Sf(j)*Fc(j)
            ENDDO
            Sn(i) = -sumd
         ENDDO
         !rsolve does not change M:
         CALL rsolve(n,M,M2,Sn)

         if(restart) call testing(n,jac_copy,sn,work,fc)
          if(restart) print *,"****************************restart****"


         IF(global .eq. 2 .or. global .eq. 3)THEN
c            copy R transpose into lower triangle of M:
            DO i =2,n
               Do j = 1,i-1
                  M(i,j) = M(j,i)
               ENDDO
            ENDDO
         ENDIF      

          IF(global .eq. 2)THEN
c             put L*Ltranspose into H in symmetric storage mode
             DO i =1,n
                DO j =i,n
                   sumd =0.0D0
                   DO k = 1,i
                      sumd =sumd +M(i,k)*M(j,k)
                   ENDDO
                   k = ((i-1)*(2*n-i)+2*j)/2
                   H(k) =sumd
                ENDDO
             ENDDO
          ENDIF
      ENDIF

      RETURN
      END
      


      subroutine  testing(n,jac,sn,work,fc)
      real *8 jac(n,n),sn(n),work(n),fc(n),ipiv(n)
      real *8 jac_c(n,n),sumd
      work(:) = -fc(:)
      jac_c(:,:) = jac(:,:)
      call dgetrf(n,n,jac,n,ipiv,info)
      call dgetrs('N',n,1,jac,n,ipiv,work,n,info)
      do i=1,n
c       print *,'i,work(i)',i,work(i)
      enddo
       do i=1,n
          sumd =0.0d0
           do j=1,n
             sumd =sumd + jac_c(i,j)*work(j)
           enddo
             sumd =sumd +fc(i)
c           print *,'i,ans =',i,sumd
        enddo



       do i=1,n
          sumd =sumd + (work(i) - sn(i))**2
       enddo
       if(sumd .gt. 1.e-12)then
          print *,'sumd =',sumd
          print *,'exit taken in testing'
           call exit(0)
       endif
      return
      end

