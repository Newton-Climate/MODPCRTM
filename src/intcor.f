      SUBROUTINE  INTCOR(EPSIL, FBEAM, FLYR, LAYRU, LYRCUT, MAXCOE,     &
     &                    MXCLY, MXUMU, NCOEF, NCUT, NPHI, NSTR,        &
     &                    NTAU, NUMU, OPRIM, PHASA, PHASE, PHASM,       &
     &                    PHIRAD, PI, PMOM, SSALB, TAUC, TAUCPR,        &
     &                    UMU, UMU0, UTAU, UTAUPR, UU)

!             CORRECT INTENSITY FIELD BY USING NAKAJIMA-TANAKA ALGORITHM
!         (1988).  FOR MORE DETAILS, SEE SECTION 3.6 OF STW NASA REPORT.

!                I N P U T   V A R I A B L E S

!        EPSIL   10 TIMES MACHINE PRECISION
!        FBEAM   INCIDENT BEAM RADIATION AT TOP
!        FLYR    TRUNCATED FRACTION IN DELTA-M METHOD.
!        LAYRU   INDEX OF UTAU IN MULTI-LAYERED SYSTEM
!        LYRCUT  LOGICAL FLAG FOR TRUNCATION OF COMPUTATIONAL LAYER
!        NCOEF   NUMBER OF PHASE FUNCTION LEGENDRE COEFFICIENTS SUPPLIED
!        NCUT    TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
!        NPHI    NUMBER OF USER AZIMUTHAL ANGLES
!        NSTR    NUMBER OF POLAR QUADRATURE ANGLES
!        NTAU    NUMBER OF USER-DEFINED OPTICAL DEPTHS
!        NUMU    NUMBER OF USER POLAR ANGLES
!        OPRIM   DELTA-M-SCALED SINGLE-SCATTER ALBEDO
!        PHIRAD  AZIMUTHAL ANGLES
!        PMOM    PHASE FUNCTION LEGENDRE COEFFICIENTS (K, LC)
!                    K=0 TO NCOEF, LC=1 TO NLYR WITH PMOM(0,LC)=1
!        SSALB   SINGLE SCATTERING ALBEDO AT COMPUTATIONAL LAYERS
!        TAUC    OPTICAL THICKNESS AT COMPUTATIONAL LEVELS
!        TAUCPR  DELTA-M-SCALED OPTICAL THICKNESS
!        UMU     COSINE OF EMERGENT ANGLE
!        UMU0    COSINE OF INCIDENT ZENITH ANGLE
!        UTAU    USER DEFINED OPTICAL DEPTHS
!        UTAUPR  DELTA-M-SCALED VERSION OF UTAU

!                O U T P U T   V A R I A B L E S

!        UU      CORRECTED INTENSITY FIELD; UU(IU,LU,J)
!                          IU=1,NUMU; LU=1,NTAU; J=1,NPHI

!                I N T E R N A L   V A R I A B L E S

!        CS      COSINE OF SCATTERING ANGLE
!        IMTHD   METHOD FLAG:  1: MS, 2: TMS, 3, IMS
!        PHASA   ACTUAL PHASE FUNCTION
!        PHASM   DELTA-M-SCALED PHASE FUNCTION
!        PHASE   SCRATCH ARRAY, PHASE FUNCTION MULTIPLIED BY
!                VARIOUS FACTORS

!   ROUTINES CALLED:  SECSCA, SINSCA

      LOGICAL  LYRCUT
      INTEGER  LAYRU(*)
      REAL     FLYR(*), OPRIM(*), PHASA(*), PHASE(*),                   &
     &         PHASM(*), PHIRAD(*), PMOM(0:MAXCOE,*),                   &
     &         SSALB(*), TAUC(0:*), TAUCPR(0:*), UMU(*),                &
     &         UTAU(*), UTAUPR(*), UU(MXUMU,MXCLY,*)

!                                   LOOP OVER POLAR AND AZIMUTHAL ANGLES

      DO IU=1, NUMU
!                          METHOD OF CORRECTION (STW SEC 3.6 FOR DETAIL)
         IMTHD=0
         IF (UMU(IU).GT.0.)  IMTHD=2
         IF (UMU(IU).LT.0.)  THEN
            IMTHD=3
!                                              SET AUREOLE +/- 10 DEGREE
            IF (ABS(UMU(IU)+UMU0) .LE. COS(10.))  THEN
               IF (TAUC(NCUT).GT.3. .AND. TAUC(NCUT).LT.8.)  IMTHD=1
               IF (TAUC(NCUT).GE.8.)  IMTHD=0
            END IF
         END IF
         IF (IMTHD.EQ.0)  GO TO 200

         DO JP=1, NPHI
!                                                  COS(SCATTERING ANGLE)

            CS=-UMU0*UMU(IU) + SQRT((1.-UMU0**2)*(1.-UMU(IU)**2)) *     &
     &                           COS(PHIRAD(JP))

!                                              INITIALIZE PHASE FUNCTION
            DO LC=1, NCUT
               PHASA(LC)=1.
               PHASM(LC)=1.
            ENDDO
!                                   INITIALIZE LEGENDRE POLY. RECURRENCE
            PL1=1.
            PL2=0.
            DO K=1, NCOEF
!                                              LEGENDRE POLY. RECURRENCE

               PL =((2*K-1) * CS * PL1 - (K-1) * PL2) / K
               PL2=PL1
               PL1=PL
!                                        CALCULATE ACTUAL PHASE FUNCTION
               DO LC=1, NCUT
                  PHASA(LC)=PHASA(LC) + (2*K+1) * PL * PMOM(K,LC)
               ENDDO

!                                       DELTA-M TRUNCATED PHASE FUNCTION
               IF(K.LT.NSTR)  THEN
                  DO LC=1, NCUT
                     PHASM(LC)=PHASM(LC) + (2*K+1) * PL *               &
     &                        (PMOM(K,LC)-FLYR(LC)) / (1.-FLYR(LC))
                  ENDDO
               END IF
            ENDDO
!                                                        ** MS METHOD **
!                                               US & UST OF EQ. STW (67)
            IF (IMTHD.EQ.1)  THEN
               DO LC=1, NCUT
                  PHASE(LC)=(1.-FLYR(LC)) * PHASM(LC)
               ENDDO
               DO LU=1, NTAU
                  IF(.NOT.LYRCUT .OR. LAYRU(LU).LE.NCUT)THEN
                     US =FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT, &
     &                  PHASA, SSALB, TAUC, UMU(IU), UMU0, UTAU(LU))
                     UTS=FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT, &
     &                  PHASE, SSALB, TAUC, UMU(IU), UMU0, UTAU(LU))

                     UU(IU,LU,JP)=UU(IU,LU,JP) + US - UTS
                  ENDIF
               ENDDO

               GO TO 100
            END IF
!                                                       ** TMS METHOD **
!                                             UPTS & UPS OF EQ. STW (68)
            DO LC=1, NCUT
               PHASE(LC)=PHASA(LC) / (1.-FLYR(LC)*SSALB(LC))
            ENDDO
            DO LU=1, NTAU
               IF(.NOT.LYRCUT .OR. LAYRU(LU).LE.NCUT)THEN
                  UPTS=FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT,   &
     &                PHASE, SSALB, TAUCPR, UMU(IU), UMU0, UTAUPR(LU))
                  UPS =FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT,   &
     &                PHASM, OPRIM, TAUCPR, UMU(IU), UMU0, UTAUPR(LU))

                  UU(IU,LU,JP)=UU(IU,LU,JP) + UPTS - UPS
               ENDIF
            ENDDO
!                                                       ** IMS METHOD **
!                     CORRECTION OF SECONDARY SCATTERING BELOW TOP LEVEL

            IF (IMTHD.EQ.3)  THEN
               LTAU=1
               IF (UTAU(1).LE.EPSIL)  LTAU=2
               DO LU=LTAU, NTAU
                  IF(.NOT.LYRCUT .OR. LAYRU(LU).LE.NCUT)THEN
                     UHAT=FBEAM/(4.*PI) * SECSCA(CS, FLYR, LAYRU(LU),   &
     &                   MAXCOE, NCOEF, NSTR, PMOM, SSALB, TAUC,        &
     &                   UMU(IU), UMU0, UTAU(LU))

                     UU(IU,LU,JP)=UU(IU,LU,JP) - UHAT
                  ENDIF
               ENDDO
            END IF
!                                                   END OF ANGULAR LOOPS
100      CONTINUE
         ENDDO
200   CONTINUE
      ENDDO

      RETURN
      END
