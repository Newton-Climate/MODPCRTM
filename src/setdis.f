      SUBROUTINE  SETDIS(CMU, CWT, DELTAM, DTAUC, EXPBEA, FBEAM, FLYR,  &
     &                    GL, IBCND, LAYRU, LYRCUT, MAXCOE,             &
     &                    MXUMU, MXCMU, NCUT, NLYR, NTAU, NN, NSTR,     &
     &                    PLANK, NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC,&
     &                    TAUCPR, UTAU, UTAUPR, NLOS, LC_OBS, OBSTAU,   &
     &                    OTAUPR, UMU, UMU0, USRTAU, USRANG)

!                            PERFORM MISCELLANEOUS SETTING-UP OPERATIONS
      IMPLICIT NONE

!       ROUTINES CALLED:  ERRMSG, QGAUSN

!       INPUT :  ALL ARE DISORT INPUT VARIABLES (SEE DOC FILE)
!                MXCMU       MAXIMUM NUMBER OF DISORT STREAMS.
!                MAXCOE      MAXIMUM NUMBER OF LEGENDRE COEFFICIENTS.
!                NLOS        NUMBER OF LINE-OF-SIGHT PATHS.
!                LC_OBS      LAYER INDEX OF OBSERVER
!                OBSTAU      OBSERVER OPTICAL DEPTH

!       OUTPUT:  NTAU,UTAU   IF USRTAU=FALSE
!                NUMU,UMU    IF USRANG=FALSE
!                CMU,CWT     COMPUTATIONAL POLAR ANGLES AND
!                            CORRESPONDING QUADRATURE WEIGHTS
!                EXPBEA      TRANSMISSION OF DIRECT BEAM
!                FLYR        TRUNCATED FRACTION IN DELTA-M METHOD
!                GL          PHASE FUNCTION LEGENDRE COEFFICIENTS MULTI-
!                            PLIED BY (2L+1) AND SINGLE-SCATTER ALBEDO
!                HLPR        LEGENDRE MOMENTS OF SURFACE BIDIRECTIONAL
!                            REFLECTIVITY, TIMES 2K+1
!                LAYRU       COMPUTATIONAL LAYER IN WHICH UTAU FALLS
!                LYRCUT      FLAG AS TO WHETHER RADIATION WILL BE ZEROED
!                            BELOW LAYER NCUT
!                NCUT        COMPUTATIONAL LAYER WHERE ABSORPTION
!                            OPTICAL DEPTH FIRST EXCEEDS  ABSCUT
!                NN          NSTR / 2
!                OPRIM       DELTA-M-SCALED SINGLE-SCATTER ALBEDO
!                TAUCPR      DELTA-M-SCALED OPTICAL DEPTH
!                UTAUPR      DELTA-M-SCALED VERSION OF  UTAU
!                OTAUPR      DELTA-M-SCALED VERSION OF OBSTAU

      LOGICAL  DELTAM, LYRCUT, PLANK, ONLYFL, USRTAU, USRANG
      INTEGER  IBCND, MAXCOE, MXUMU, MXCMU, NCUT, NLYR, NTAU, NN, NSTR, &
     &         NUMU, NLOS, LC_OBS(0:NLOS), LAYRU(*)
      REAL     FBEAM, UMU0
      REAL     CMU(*), CWT(*), DTAUC(*), EXPBEA(0:*),                   &
     &         FLYR(*), GL(0:MXCMU,*), OPRIM(*),                        &
     &         PMOM(0:MAXCOE,*), SSALB(*), TAUC(0:*), TAUCPR(0:*),      &
     &         UTAU(*), UTAUPR(*), UMU(*), OBSTAU(NLOS), OTAUPR(NLOS)

!     LOCAL VARIABLES:
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
      INTEGER  ILOS, LC, K, LU, IQ, IU
      REAL     ABSTAU, F
      REAL     ABSCUT
      DATA     ABSCUT / 20. /

!                    SET OUTPUT LEVELS AT COMPUTATIONAL LAYER BOUNDARIES

      IF (.NOT.USRTAU)  THEN
         NTAU=NLYR + 1
         DO LC=0, NTAU-1
            UTAU(LC+1)=TAUC(LC)
         ENDDO
      END IF
!                             APPLY DELTA-M SCALING AND MOVE DESCRIPTION
!                             OF COMPUTATIONAL LAYERS TO LOCAL VARIABLES
      EXPBEA(0)=1.0
      TAUCPR(0)=0.0
      ABSTAU=0.0
      DO LC=1, NLYR
         PMOM(0,LC)=1.0
         IF (ABSTAU.LT.ABSCUT)  NCUT=LC
         IF(.NOT.PLANK .OR. LC_OBS(0).LT.LC)                            &
     &      ABSTAU=ABSTAU+(1-SSALB(LC))*DTAUC(LC)

         IF (.NOT.DELTAM)  THEN
            OPRIM(LC)=SSALB(LC)
            TAUCPR(LC)=TAUC(LC)
            DO K=0, NSTR-1
               GL(K,LC)=(2*K+1) * OPRIM(LC) * PMOM(K,LC)
            ENDDO
            F=0.0
         ELSE
!                                                 DELTA-M TRANSFORMATION
            F=PMOM(NSTR,LC)
            OPRIM(LC)=SSALB(LC) * (1. - F) / (1. - F * SSALB(LC))
            TAUCPR(LC)=TAUCPR(LC-1) + (1. - F*SSALB(LC)) * DTAUC(LC)
            DO K=0, NSTR-1
               GL(K,LC)=(2*K+1) * OPRIM(LC) * (PMOM(K,LC)-F) / (1.-F)
            ENDDO
         END IF

         FLYR(LC)=F
         EXPBEA(LC)=0.0
         IF (FBEAM.GT.0.0)  EXPBEA(LC)=EXP(- TAUCPR(LC) / UMU0)
      ENDDO
!                IF NO THERMAL EMISSION, CUT OFF MEDIUM BELOW ABSORPTION
!               OPTICAL DEPTH=ABSCUT (NOTE THAT DELTA-M TRANSFORMATION
!                 LEAVES ABSORPTION OPTICAL DEPTH INVARIANT).  NOT WORTH
!                            THE TROUBLE FOR ONE-LAYER PROBLEMS, THOUGH.
      LYRCUT=ABSTAU.GE.ABSCUT .AND. IBCND.NE.1 .AND. NLYR.GT.1
      IF(.NOT.LYRCUT)NCUT=NLYR

!                     SET ARRAYS DEFINING LOCATION OF USER OUTPUT LEVELS
!                               WITHIN DELTA-M-SCALED COMPUTATIONAL MESH
      DO LU=1, NTAU
         DO LC=1, NLYR
            IF (UTAU(LU).GE.TAUC(LC-1) .AND. UTAU(LU).LE.TAUC(LC))      &
     &           GO TO 60
         ENDDO
         LC=NLYR

60       UTAUPR(LU)=UTAU(LU)
         IF(DELTAM) UTAUPR(LU)=TAUCPR(LC-1) + (1-SSALB(LC)*FLYR(LC))    &
     &                                        * (UTAU(LU) - TAUC(LC-1))
         LAYRU(LU)=LC
      ENDDO
      IF(DELTAM)THEN
          DO ILOS=1,NLOS
              IF(LC_OBS(ILOS).GT.0)OTAUPR(ILOS)=TAUCPR(LC_OBS(ILOS)-1)  &
     &          +(1-SSALB(LC_OBS(ILOS))*FLYR(LC_OBS(ILOS)))             &
     &          *(OBSTAU(ILOS)-TAUC(LC_OBS(ILOS)-1))
          ENDDO
      ELSE
          DO ILOS=1,NLOS
              IF(LC_OBS(ILOS).GT.0)OTAUPR(ILOS)=OBSTAU(ILOS)
          ENDDO
      ENDIF
!                        CALCULATE COMPUTATIONAL POLAR ANGLE COSINES AND
!                             ASSOCIATED QUADRATURE WEIGHTS FOR GAUSSIAN
!                              QUADRATURE ON THE INTERVAL (0,1) (UPWARD)
      NN=NSTR / 2
      CALL  QGAUSN(NN, CMU, CWT)
!                                      DOWNWARD (NEG) ANGLES AND WEIGHTS
      DO IQ=1, NN
         CMU(IQ+NN)=- CMU(IQ)
         CWT(IQ+NN)=  CWT(IQ)
      ENDDO
!                             COMPARE BEAM ANGLE TO COMPUTATIONAL ANGLES

      IF (FBEAM.GT.0.0)  THEN
         DO IQ=1, NN
            IF (ABS(UMU0-CMU(IQ))/UMU0 .LT. 1.E-4)  CALL ERRMSG         &
     &         ('SETDIS--beam angle=computational angle; change NSTR',  &
     &            .TRUE.)
         ENDDO
      END IF
!                  SET OUTPUT POLAR ANGLES TO COMPUTATIONAL POLAR ANGLES

      IF (.NOT.USRANG .OR. (ONLYFL .AND. MXUMU.GE.NSTR))  THEN
!M
!M       ALWAYS SKIP THIS BRANCH FOR MODTRAN:
!M       (THIS USE TO HAPPEN AUTOMATICALLY SINCE USRANG IS .TRUE. AND
!M       SINCE MXUMU EQUALED 1 AND NSTR WAS AT LEAST 2; TO ACCOMMODATE
!M       IBCND=1 CALCULATIONS, MXUMU HAS BEEN INCREASED TO 2).
!M       NUMU=NSTR
!M       DO IU=1, NN
!M          UMU(IU)=- CMU(NN+1-IU)
!M       ENDDO
!M       DO IU=NN+1, NSTR
!M          UMU(IU)=CMU(IU-NN)
!M       ENDDO
      END IF
!                             SHIFT POSITIVE USER ANGLE COSINES TO UPPER
!                         LOCATIONS AND PUT NEGATIVES IN LOWER LOCATIONS

      IF (USRANG .AND. IBCND.EQ.1)  THEN
         DO IU=1, NUMU
            UMU(IU+NUMU)=UMU(IU)
         ENDDO
         DO IU=1, NUMU
            UMU(IU)=- UMU(2*NUMU+1-IU)
         ENDDO
         NUMU=2*NUMU
      END IF

      RETURN
      END
