C
C
C
C
C
C
C
C
C
      SUBROUTINE MACHINE_EPSD(EPSILOND)
C-----------------------------------------------------------------------
C---GET THE (DOUBLE PRECISION) MACHINE EPSILON
C--------------------------------------------------------------------HSJ
C
      REAL *8 EPSILOND
      EPSILOND=1.0D0
      DO WHILE (1.D0+EPSILOND .GT. 1.0D0)
         EPSILOND=EPSILOND*0.5D0
      ENDDO
C
C
      RETURN
      END
