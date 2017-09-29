      REAL FUNCTION  SINSCA(EPSIL, LAYRU, NLYR, PHASE, OMEGA, TAU, UMU, &
     &                       UMU0, UTAU)

!                       SINGLY SCATTERED INTENSITY OF EQS. STW (65B,D,E)

!                I N P U T   V A R I A B L E S

!        EPSIL    10 TIMES MACHINE PRECISION
!        LAYRU    INDEX OF UTAU IN MULTI-LAYERED SYSTEM
!        NLYR     NUMBER OF SUBLAYERS
!        PHASE    PHASE FUNCTIONS OF SUBLAYERS
!        OMEGA    SINGLE SCATTERING ALBEDOS OF SUBLAYERS
!        TAU      OPTICAL THICKNESSES OF SUBLAYERS
!        UMU      COSINE OF EMERGENT ANGLE
!        UMU0     COSINE OF INCIDENT ZENITH ANGLE
!        UTAU     USER DEFINED OPTICAL DEPTH FOR OUTPUT INTENSITY

      REAL        PHASE(*), OMEGA(*), TAU(0:*)
!                                                         INITIALIZATION
      SINSCA=0.
      EXP0  =EXP(-UTAU/UMU0)
!                                      EQ. STW (65E): DIRECT TRANSMITTED

      IF (ABS(UMU+UMU0).LE.EPSIL)  THEN
         DO LYR=1, LAYRU-1
            SINSCA=SINSCA + OMEGA(LYR)*PHASE(LYR)*(TAU(LYR)-TAU(LYR-1))
         ENDDO
         SINSCA=EXP0 / UMU0 * (SINSCA + OMEGA(LAYRU)*PHASE(LAYRU)       &
     &                                   * (UTAU-TAU(LAYRU-1)))
         RETURN
      END IF
!                                               EQ. STW (65B): REFLECTED
      IF (UMU.GT.0.)  THEN
         DO LYR=LAYRU, NLYR
            EXP1=EXP(- ((TAU(LYR)-UTAU)/UMU + TAU(LYR)/UMU0))
            SINSCA=SINSCA + OMEGA(LYR) * PHASE(LYR) * (EXP0-EXP1)
            EXP0=EXP1
         ENDDO
      ELSE
!                                     EQ. STW (65D): DIFFUSE TRANSMITTED
         DO LYR=LAYRU, 1, -1
            EXP1=EXP(- ((TAU(LYR-1)-UTAU)/UMU + TAU(LYR-1)/UMU0))
            SINSCA=SINSCA + OMEGA(LYR) * PHASE(LYR) * (EXP0-EXP1)
            EXP0=EXP1
         ENDDO
      END IF

      SINSCA=SINSCA / (1. + UMU/UMU0)

      RETURN
      END
