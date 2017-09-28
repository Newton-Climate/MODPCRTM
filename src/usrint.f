      SUBROUTINE  USRINT(BPLANK, CMU, CWT, DELM0, EMU, EXPBEA, FBEAM,   &
     &                    FISOT, GC, GU, KK, LAMBER, LAYRU, LL, LYRCUT, &
     &                    MAZIM, MXCMU, MXUMU, NCUT, NLYR, NN, NSTR,    &
     &                    PLANK, NUMU, NTAU, PI, RMU, TAUCPR, TPLANK,   &
     &                    UMU, UMU0, UTAUPR, WK, ZBEAM, Z0U, Z1U, ZZ,   &
     &                    ZPLK0, ZPLK1, UUM)

!                     Compute intensity components at user output angles
!                             for azimuthal expansion terms in EQ. SD(2)

!                 I N P U T    V A R I A B L E S:

!       BPLANK :  integrated Planck function for emission from
!                 bottom boundary
!       CMU    :  abscissae for Gauss quadrature over angle cosine
!       CWT    :  weights for Gauss quadrature over angle cosine
!       DELM0  :  Kronecker delta, delta-sub-M0
!       EMU    :  surface directional emissivity (user angles)
!       EXPBEA :  transmission of incident beam, EXP(-TAUCPR/UMU0)
!       GC     :  eigenvectors at polar quadrature angles, SC(1)
!       GU     :  eigenvectors interpolated to user polar angles
!                 (i.e., g in EQ. SC(1))
!       KK     :  eigenvalues of coeff. matrix in EQ. SS(7)
!       LAYRU  :  layer number of user level -UTAU-
!       LL     :  constants of integration in EQ. SC(1), obtained
!                 by solving scaled version of EQ. SC(5);
!                 exponential term of EQ. SC(12) not included
!       LYRCUT :  logical flag for truncation of computational layer
!       MAZIM  :  order of azimuthal component
!       NCUT   :  total number of computational layers considered
!       NN     :  order of double-Gauss quadrature (NSTR/2)
!       RMU    :  surface bidirectional reflectivity (user angles)
!       TAUCPR :  cumulative optical depth (delta-M-Scaled)
!       TPLANK :  integrated Planck function for emission from
!                 top boundary
!       UTAUPR :  optical depths of user output levels in delta-M
!                 coordinates;  equal to -UTAU- if no delta-M
!       Z0U    :  Z-sub-zero in EQ. SS(16) interpolated to user
!                 angles from an equation derived from SS(16)
!       Z1U    :  Z-sub-one in EQ. SS(16) interpolated to user
!                 angles from an equation derived from SS(16)
!       ZZ     :  beam source vectors in EQ. SS(19)
!       ZPLK0  :  thermal source vectors -Z0-, by solving EQ. SS(16)
!       ZPLK1  :  thermal source vectors -Z1-, by solving EQ. SS(16)
!       ZBEAM  :  incident-beam source vectors
!       (Remainder are 'disort' input variables)

!                 O U T P U T    V A R I A B L E S:

!       UUM    :  azimuthal components of the intensity in EQ. STWJ(5)

!                 I N T E R N A L    V A R I A B L E S:

!       BNDDIR :  direct intensity down at the bottom boundary
!       BNDDFU :  diffuse intensity down at the bottom boundary
!       BNDINT :  intensity attenuated at both boundaries, STWJ(25-6)
!       DTAU   :  optical depth of a computational layer
!       LYREND :  end layer of integration
!       LYRSTR :  start layer of integration
!       PALINT :  intensity component from parallel beam
!       PLKINT :  intensity component from planck source
!       WK     :  scratch vector for saving 'EXP' evaluations

!       All the exponential factors (EXP1, EXPN,... etc.)
!       come from the substitution of constants of integration in
!       EQ. SC(12) into EQS. S1(8-9).  They all have negative
!       arguments so there should never be overflow problems.

!+---------------------------------------------------------------------+

      LOGICAL  LAMBER, LYRCUT, PLANK, NEGUMU
      INTEGER  LAYRU(*)
      REAL     CMU(*), CWT(*), EMU(*), EXPBEA(0:*),                     &
     &  GC(MXCMU,MXCMU,*),GU(0:MXUMU,MXCMU,*),KK(MXCMU,*),              &
     &         LL(MXCMU,*), RMU(MXUMU,0:*), TAUCPR(0:*),                &
     &         UUM(MXUMU,*), UMU(*), UTAUPR(*), WK(*),                  &
     &         Z0U(MXUMU,*), Z1U(MXUMU,*), ZBEAM(MXUMU,*),              &
     &         ZZ(MXCMU,*), ZPLK0(MXCMU,*), ZPLK1(MXCMU,*)

!    Incorporate constants of integration into interpolated eigenvectors

      DO LC=1, NCUT
         DO IQ=1, NSTR
            DO IU=1, NUMU
               GU(IU,IQ,LC)=GU(IU,IQ,LC) * LL(IQ,LC)
            ENDDO
         ENDDO
      ENDDO
!                                  Loop over levels at which intensities
!                                     are desired ('user output levels')
      DO LU=1, NTAU

!        MODTRAN FIX:  CHECK FOR POSITIVE FBEAM.
         IF(FBEAM.GT.0.)EXP0=EXP(- UTAUPR(LU) / UMU0)
         LYU=LAYRU(LU)

!                Loop over polar angles at which intensities are desired

         DO IU=1, NUMU
            IF (LYRCUT .AND. LYU.GT.NCUT)  GO TO 100
            NEGUMU=UMU(IU).LT.0.0
            IF(NEGUMU)  THEN
               LYRSTR=1
               LYREND=LYU - 1
               SGN=- 1.0
            ELSE
               LYRSTR=LYU + 1
               LYREND=NCUT
               SGN=1.0
            END IF

!               For downward intensity, integrate from top to 'LYU-1' in
!       EQ. S1(8); for upward, integrate from bottom to 'LYU+1' in S1(9)

            PALINT=0.0
            PLKINT=0.0
            DO LC=LYRSTR, LYREND

               DTAU=TAUCPR(LC) - TAUCPR(LC-1)
               EXP1= EXP((UTAUPR(LU) - TAUCPR(LC-1)) / UMU(IU))
               EXP2= EXP((UTAUPR(LU) - TAUCPR(LC)) / UMU(IU))

               IF (PLANK .AND. MAZIM.EQ.0)                              &
     &           PLKINT=PLKINT + SGN * (Z0U(IU,LC) * (EXP1 - EXP2) +    &
     &                    Z1U(IU,LC) * ((TAUCPR(LC-1) + UMU(IU))*EXP1 - &
     &                                   (TAUCPR(LC) + UMU(IU))*EXP2))

               IF (FBEAM.GT.0.0)  THEN
                  DENOM=1.0 + UMU(IU) / UMU0
                  IF (ABS(DENOM).LT.0.0001)  THEN
!                                                       L'Hospital limit
                     EXPN=(DTAU / UMU0) * EXP0
                  ELSE
                     EXPN=(EXP1 * EXPBEA(LC-1) - EXP2 * EXPBEA(LC))     &
     &                      * SGN / DENOM
                  END IF
                  PALINT=PALINT + ZBEAM(IU,LC) * EXPN
               END IF
!                                                          -KK- negative
               DO IQ=1, NN
                  WK(IQ)=EXP(KK(IQ,LC) * DTAU)
                  DENOM=1+UMU(IU)*KK(IQ,LC)
                  IF (ABS(DENOM).LT.0.0001)  THEN
!                                                       L'Hospital limit
                     EXPN=DTAU / UMU(IU) * EXP2
                  ELSE
                     EXPN=SGN * (EXP1 * WK(IQ) - EXP2) / DENOM
                  END IF
                  PALINT=PALINT + GU(IU,IQ,LC) * EXPN
               ENDDO
!                                                          -KK- positive
               DO IQ=NN+1, NSTR
                  DENOM=1.0 + UMU(IU) * KK(IQ,LC)
                  IF (ABS(DENOM).LT.0.0001)  THEN
!                                                       L'Hospital limit
                     EXPN=- DTAU / UMU(IU) * EXP1
                  ELSE
                     EXPN=SGN *(EXP1 - EXP2 * WK(NSTR+1-IQ)) / DENOM
                  END IF
                  PALINT=PALINT + GU(IU,IQ,LC) * EXPN
               ENDDO

            ENDDO
!                                Calculate contribution from user output
!                                      level to next computational level

            DTAU1=UTAUPR(LU) - TAUCPR(LYU-1)
            DTAU2=UTAUPR(LU) - TAUCPR(LYU)
            IF(ABS(DTAU1).LT.1.E-6 .AND. NEGUMU)  GO TO 50
            IF(ABS(DTAU2).LT.1.E-6 .AND. (.NOT.NEGUMU))  GO TO 50
            IF(NEGUMU) EXP1=EXP(DTAU1 / UMU(IU))
            IF(.NOT.NEGUMU) EXP2=EXP(DTAU2 / UMU(IU))

            IF (FBEAM.GT.0.0)  THEN
               DENOM=1.0 + UMU(IU) / UMU0
               IF (ABS(DENOM).LT.0.0001)  THEN
                  EXPN= (DTAU1 / UMU0) * EXP0
               ELSE IF (NEGUMU) THEN
                  EXPN=(EXP0 - EXPBEA(LYU-1) * EXP1) / DENOM
               ELSE
                  EXPN=(EXP0 - EXPBEA(LYU) * EXP2) / DENOM
               END IF
               PALINT=PALINT + ZBEAM(IU,LYU) * EXPN
            ENDIF
!                                                          -KK- negative
            DTAU=TAUCPR(LYU) - TAUCPR(LYU-1)
            DO IQ=1, NN
               DENOM=1.0 + UMU(IU) * KK(IQ,LYU)
               IF (ABS(DENOM).LT.0.0001)  THEN
                  EXPN=- DTAU2 / UMU(IU) * EXP2
               ELSE IF (NEGUMU) THEN
                  EXPN=(EXP(- KK(IQ,LYU) * DTAU2) -                     &
     &                     EXP(KK(IQ,LYU) * DTAU) * EXP1) / DENOM
               ELSE
                  EXPN=(EXP(- KK(IQ,LYU) * DTAU2) - EXP2) / DENOM
               END IF
               PALINT=PALINT + GU(IU,IQ,LYU) * EXPN
            ENDDO
!                                                          -KK- positive
            DO IQ=NN+1, NSTR
               DENOM=1.0 + UMU(IU) * KK(IQ,LYU)
               IF (ABS(DENOM).LT.0.0001)  THEN
                  EXPN=- DTAU1 / UMU(IU) * EXP1
               ELSE IF (NEGUMU) THEN
                  EXPN=(EXP(- KK(IQ,LYU) * DTAU1) - EXP1) / DENOM
               ELSE
                  EXPN=(EXP(- KK(IQ,LYU) * DTAU1) -                     &
     &                     EXP(- KK(IQ,LYU) * DTAU) * EXP2) / DENOM
               END IF
               PALINT=PALINT + GU(IU,IQ,LYU) * EXPN
            ENDDO

            IF (PLANK .AND. MAZIM.EQ.0)  THEN
              IF (NEGUMU)  THEN
                 EXPN=EXP1
                 FACT=TAUCPR(LYU-1) + UMU(IU)
              ELSE
                 EXPN=EXP2
                 FACT=TAUCPR(LYU) + UMU(IU)
              END IF
              PLKINT=PLKINT + Z0U(IU,LYU) * (1.- EXPN) +                &
     &                 Z1U(IU,LYU) *(UTAUPR(LU) + UMU(IU) - FACT*EXPN)
            END IF

!           Calculate intensity components attenuated at both boundaries
!           Note: no azimuthal intensity component for isotropic surface

50          BNDINT=0.0
            IF (NEGUMU .AND. MAZIM.EQ.0)  THEN
              BNDINT=(FISOT + TPLANK) * EXP(UTAUPR(LU) / UMU(IU))
            ELSE IF (.NOT.NEGUMU)  THEN
              IF (LYRCUT .OR. (LAMBER .AND. MAZIM.GT.0))  GO TO 90
              DO JQ=NN+1, NSTR
                 WK(JQ)=EXP(-KK(JQ,NLYR)*(TAUCPR(NLYR)-TAUCPR(NLYR-1)))
              ENDDO
              BNDDFU=0.0
              DO IQ=NN, 1, -1
                 DFUINT=0.0
                 DO JQ=1, NN
                    DFUINT=DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
                 ENDDO
                 DO JQ=NN+1, NSTR
                    DFUINT=DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)*WK(JQ)
                 ENDDO
                 IF (FBEAM.GT.0.0)                                      &
     &                DFUINT=DFUINT + ZZ(IQ,NLYR) * EXPBEA(NLYR)
                 DFUINT=DFUINT + DELM0 * (ZPLK0(IQ,NLYR)                &
     &                                + ZPLK1(IQ,NLYR) * TAUCPR(NLYR))
                 BNDDFU=BNDDFU + (1. + DELM0) * RMU(IU,NN+1-IQ)         &
     &                           * CMU(NN+1-IQ) * CWT(NN+1-IQ) * DFUINT
              ENDDO

              BNDDIR=0.0
!***************** VINCENT ROSS CHANGED FOR BRDF COUPLING *************
#ifndef BRDF_COUPLING
              IF (FBEAM.GT.0.0) BNDDIR=UMU0*FBEAM/PI * RMU(IU,0)        &
     &                                   * EXPBEA(NLYR)
#endif
!************************** END VINCENT ROSS **************************
!***************** VINCENT ROSS CHANGED FOR BRDF COUPLING *************
#ifdef BRDF_COUPLING
              BNDINT= BNDDFU
#else
              BNDINT=(BNDDFU + BNDDIR + DELM0 * EMU(IU) * BPLANK)       &
#endif
     &                 * EXP((UTAUPR(LU)-TAUCPR(NLYR)) / UMU(IU))
!************************** END VINCENT ROSS **************************
            END IF

90          UUM(IU,LU)=PALINT + PLKINT + BNDINT
100      CONTINUE
         ENDDO
      ENDDO

      RETURN
      END
