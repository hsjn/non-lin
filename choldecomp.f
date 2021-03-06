C
C
C
      SUBROUTINE CHOLDECOMP(N,H,POSDEF,EPSILOND,L,ADDMAX)
C----------------------------------------------------------------------
C  CALCULATE THE PERTURBED CHOLESKY DECOMPOSITION OF THE SYMMETRIC
C  MATRIX H. H MUST BE IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE.
C  (IE H(I,J),FOR J  .GE. I,IS ACCESSED AS H(K),WHERE
C       K=((I-1)(2*N-I)+2*J)/2
C   MATRIX H SHOULD BE POSITIVE DEFINITE. IF IT IS NOT THEN THE
C   DIAGONAL ELEMENS OF H ARE (IMPLICITELY) MODIFIED A MINIMAL AMOUNT UNTIL
C   H IS POSITIVE DEFINITE. NOTE THAT H IS ACTUALLY NOT CHANGED. INSTEAD
C   IT IS L THAT IS CHANGED TO ACCOUNT FOR THE NON POSITIVE DEFINITENESS
C   OF H. ON OUTPUT H WILL BE THE SAME AS ON INPUT.
C  THE LOWER TRIANGULAR MATRIX L,OF THE CHOLESKY DECOMPOSITION OF
C  H,(IE H=L*LT) IS RETURNED .
C  THIS SUBROUTINE FORMS L BY FIRST GETTING THE DIAGONAL ELEMENT
C  L(I,I) AND THEN CALCULATING THE ELEMENTS OF COLUMN I BELOW THIS
C  DIAGONAL ELEMENT (IE L(J,I),J=I+1,...N). 
C  INPUT
C  N            SQUARE SIZE OF MATRIX H
C  H            MATRIX WHICH MAY OR MAY NOT BE POSTIVE DEFINITE.
C               IF IT IS NOT THEN A DIAGONAL MATRIX D WILL BE ADDED
C               TO H TO MAKE IT POSITIVE DEFINITE. THE DIAGONAL
C               ELEMENTS OF D WILL NOT ALL BE THE SAME. THE LARGEST
C               ELEMENT OF D IS RETURNED AS ADDMAX.
C               NOTE THAT H MUST BE IN UPPER TRIANGULAR SYMMETRIC
C               STORAGE MODE AS EXPLAINED ABOVE. THE ADDITION OF
C               THE DIAGONAL ELEMENTS IS DONE IMPLICITELY,SO THAT H
C               IS ACTUALLY NEVER CHANGED.
C  POSDEF        A TOLERANCE PARAMETER.SET POSDEF=0.0 IF H IS KNOWN
C               TO BE POSITIVE DEFINITE A PRIORI.
C---OUTPUT
C  L(I,j)         LOWER TRIANGULAR CHOLESKY FACTOR OF H IN
C               Full STORAGE MODE
C  ADDMAX       THE MAXIMUM DIAGONAL ELEMENT OF MATRIX D
C               THAT WAS ADDED TO H.
C
C
C------------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      INTEGER N,N2,I,J,IHJJ,IJI,III,K
      REAL*8 MINMU,MINMU2,POSDEF,ADDMAX,
     &     HMAX,SQRTEPS,MINLJJ,SQRT,SQRT4
      REAL*8 SUMD,H(*),L(N,N),EPSILOND
      DATA SQRTEPS,SQRT4 /0.0D0,0.0D0/
C
C
      N2=2*N
      IF(SQRTEPS .EQ. 0.0D0)SQRTEPS=SQRT(EPSILOND)
      IF(SQRT4 .EQ. 0.0D0)SQRT4=SQRT(SQRTEPS)
      MINMU2=0.0
      MINMU=SQRT4*POSDEF
      IF(POSDEF .EQ. 0.0)THEN  !H is known to be postive definite
          HMAX=-1.
          DO I=1,N
              III=(I*(N2-I+3))/2-N
              HMAX=MAX(HMAX,ABS(H(III)))
          ENDDO
          POSDEF=SQRT(HMAX)
          MINMU2=SQRTEPS*POSDEF
      ENDIF
      ADDMAX=0.0D0
      DO J=1,N
C         FIRST GET THE DIAGONAL ELEMENT L(ILJJ). THE NORMALIZATION
C         OF THIS ELEMENT IS DONE AFTER WE DETERMINE IF A PERTURBATION
C         IS REQUIRED.
          SUMD=0.0D0
              DO I=1,J-1
                  SUMD=SUMD+L(j,i)**2
              ENDDO
          IHJJ=(J*(N2-J+3))/2-N
          L(J,J)=H(IHJJ)-SUMD
C         NEXT GET THE REMAINDER OF COLUMN J BELOW THE DIAGONAL.
C         MINLJJ KEEPS TRACK OF THE LARGEST (IN ABSOLUTE VALUE)
C         OFFDIAGONAL ELEMENT ENCOUNTERED IN THIS COLUMN:
C         THE NORMALIZATION OF L(I,J) IS DONE AFTER L(J,J) IS
C         NORMALIZED.
          MINLJJ=0.0
          DO I=J+1,N
              SUMD=0.0D0
                  DO K=1,J-1
                      SUMD=SUMD+L(i,k)*L(j,k)
                  ENDDO
              IJI=((J-1)*(N2-J)+2*I)/2
              L(i,j)=H(IJI)-SUMD
              MINLJJ=MAX(MINLJJ,ABS(L(I,J)))
          ENDDO
          MINLJJ=MAX(MINLJJ/POSDEF,MINMU)
C         NOW FINISH THE CALCULATION OF THE DIAGONAL TERM:
          IF(L(j,j) .GT. MINLJJ**2)THEN   !STANDARD CASE
              L(j,j)=SQRT(L(j,j))
          ELSE                            !PERTURBED CASE
              IF(MINLJJ .LT. MINMU2)MINLJJ=MINMU2
              ADDMAX=MAX(ADDMAX,MINLJJ**2-L(j,j))
              L(j,j)=MINLJJ
          ENDIF
C         FINISH THE CALCULATION OF THE REMAINDER OF COLUMN J:
          DO I=J+1,N
              L(I,J)=L(I,J)/L(J,J)
          ENDDO
      ENDDO
C
C
      RETURN
      END
