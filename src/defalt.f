      SUBROUTINE DEFALT(Z,P,T)

!     THIS SUBROUTINE INTERPOLATES PROFILES FROM THE 6 BUILT-IN MODEL
!     ATMOSPHERIC TO ALTITUDE Z.  THE JUNIT INDICES FROM /CARD1B/
!     INDICATE WHICH PROFILE SHOULD BE USED FOR EACH INTERPOLATION:

!                   JUNIT     MODEL ATMOSPHERE
!                     1          TROPICAL
!                     2          MID-LATITUDE SUMMER
!                     3          MID-LATITUDE WINTER
!                     4          HIGH-LAT SUMMER
!                     5          HIGH-LAT WINTER
!                     6          US STANDARD
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       Z        INPUT ALTITUDE [KM].
!       T        MODEL ATMOSPHERE TEMPERATURE CHANGE FROM DEFAULT [K].

!     OUTPUT ARGUMENTS:
!       P        PRESSURE FROM MODEL ATMOSPHERE [MBAR].
!       T        TEMPERATURE FROM MODEL ATMOSPHERE WITH CHANGE [K].
      DOUBLE PRECISION Z
      REAL P,T

!     COMMONS:
      INCLUDE 'YPROP.h'

!     /CARD1B/
      INTEGER JUNITP,JUNITT,JUNIT,JLOW
      REAL WMOL
      COMMON/CARD1B/JUNITP,JUNITT,JUNIT(13),WMOL(12),JLOW

!     /MLATM/
!       ALT      MODEL ATMOSPHERE ABOVE SEA LEVEL ALTITUDE LEVELS [KM].
!       PMLOG    NATURAL LOG OF PRESSURE PER MBAR FOR 6 MODEL PROFILES.
!       TMATM    6 MODEL ATMOSPHERE TEMPERATURE PROFILES [K].
!       AMOL     6 MODEL ATM DENSITY PROFILES, MOLECULES 1 TO 8 [ppmV].
      DOUBLE PRECISION ALT
      REAL PMLOG,TMATM,AMOL
      COMMON/MLATM/ALT(NLEVEL),PMLOG(NLEVEL,6),TMATM(NLEVEL,6),         &
     &  AMOL(NLEVEL,6,8)

!     /MLATMX/
      REAL AMOLX
      COMMON/MLATMX/AMOLX(NLEVEL,NMOLX)

!     /TRACE/
!       TRAC     TRACE MOLECULAR SPECIES DEFAULT PROFILES [PPMV].
      REAL TRAC
      COMMON/TRACE/TRAC(NLEVEL,21)

!     /CO2MIX/
      REAL CO2RAT
      COMMON/CO2MIX/CO2RAT

!     /CONSTN/
!       AMWT     MOLECULAR WEIGHTS [GM/MOLE].
!                  1  H2O   2  CO2   3  O3    4  N2O   5  CO    6  CH4
!                  7  O2    8  NO    9  SO2   10 NO2   11 NH3   12 HNO3
      REAL AMWT
      COMMON/CONSTN/AMWT(12)

!     /CD1BXY/
!       JUNITX   INPUT UNIT FLAG FOR CROSS-SECTION MOLECULES.
!       WMOLX    X-SECTION MOLECULAR DENSITY [PPMV AFTER JUBRAN CALL].
!       JUNITY   INPUT UNIT FLAG FOR AUXILIARY (Y) MOLECULES.
!       WMOLY    AUXILIARY MOLECULAR DENSITY [PPMV AFTER JUBRAN CALL].
      INTEGER JUNITX,JUNITY
      REAL WMOLX,WMOLY
      COMMON/CD1BXY/JUNITX,WMOLX(NMOLX),JUNITY,WMOLY(MMOLY)

!     /S_PROF/
!       S_UMIX   SCALE FACTORS FOR UNIFORMLY MIXED MOLECULAR PROFILES.
!       S_XSEC   SCALE FACTORS FOR CROSS-SECTION MOLECULAR PROFILES.
!       S_TRAC   SCALE FACTORS FOR TRACE MOLECULAR PROFILES.
!       L_UMIX   LOGICAL, .TRUE. IF S_UMIX VALUES ARE TO BE READ IN.
!       L_XSEC   LOGICAL, .TRUE. IF S_XSEC VALUES ARE TO BE READ IN.
!       L_TRAC   LOGICAL, .TRUE. IF S_TRAC VALUES ARE TO BE READ IN.
      REAL S_UMIX(4:NMOL+1),S_XSEC(NMOLX),S_TRAC(16)
      LOGICAL L_UMIX,L_XSEC,L_TRAC
      COMMON/S_PROF/S_UMIX,S_XSEC,S_TRAC,L_UMIX,L_XSEC,L_TRAC

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /MLATM/,/TRACE/,/MLATMX/,/CONSTN/
      EXTERNAL MLATMB,XMLATM

!     FUNCTIONS:
      REAL EXPINT

!     LOCAL VARIABLES:
!       KSPEC    SPECIES LOOP INDEX.
!       KSPEC7   KSPEC OFFSET BY THE NUMBER OF MODEL DEPENDENT SPECIES.
!       PMODEL   MODEL ATMOSPHERE PRESSURE AT Z [MBARS].
!       TMODEL   MODEL ATMOSPHERE TEMPERATURE AT Z [MBARS].
!       TRATIO   RATIO OF MODEL TEMPERATURE TO STANDARD TEMPERATURE.
!       MAP      MAPPING FROM ACTIVE Y-SPECIES TO DEFAULT PROFILE.
      INTEGER I0,I1,I2,I3,KSPEC,KSPEC7,MAP
      LOGICAL LINIT
      DOUBLE PRECISION Z0,Z1,Z2,Z3
      REAL FACTOR,PMODEL,TMODEL,TRATIO

!     FUNCTIONS:
!       ZLG4PT   4-POINT LAGRANGE INTERP W/ DOUBLE PRECISION ABSCISSA.
      REAL ZLG4PT

!     BEGIN CALCULATIONS:
      IF(Z.LE.ALT(2))THEN

!         INPUT ALTITUDE Z IS BELOW SECOND ALTITUDE.
          I1=1
          I2=2
      ELSEIF(Z.GE.ALT(NLEVM1))THEN

!         INPUT ALTITUDE Z IS ABOVE SECOND TO LAST ALTITUDE.
          I1=NLEVM1
          I2=NLEVEL
      ELSE

!         INTERMEDIATE ALTITUDE.
          I2=3
          DO I3=4,NLEVEL
              IF(Z.LE.ALT(I2))GOTO 10
              I2=I3
          ENDDO
   10     CONTINUE
          I1=I2-1
          I0=I2-2

!         LINEAR AND LAGRANGE 4-POINT INTERPOLATION COEFFICIENTS:
          Z0=ALT(I0)
          Z1=ALT(I1)
          Z2=ALT(I2)
          Z3=ALT(I3)
          LINIT=.TRUE.

!         TEST PRESSURE INPUT FLAG.
          IF(JUNITP.LE.6)THEN
              P=ZLG4PT(Z,Z0,Z1,Z2,Z3,PMLOG(I0,JUNITP),PMLOG(I1,JUNITP), &
     &          PMLOG(I2,JUNITP),PMLOG(I3,JUNITP),LINIT)
              LINIT=.FALSE.
              P=EXP(P)
          ENDIF

!         TEST TEMPERATURE INPUT FLAG.
          IF(JUNITT.LE.6)THEN
              TMODEL=ZLG4PT(Z,Z0,Z1,Z2,Z3,                              &
     &          TMATM(I0,JUNITT),TMATM(I1,JUNITT),                      &
     &          TMATM(I2,JUNITT),TMATM(I3,JUNITT),LINIT)
              T=TMODEL+T
              LINIT=.FALSE.
          ENDIF

!         IF JUNIT(1) IS NEGATIVE, DETERMINE DEFAULT RELATIVE HUMIDITY:
          IF(JUNIT(1).LT.0)THEN

!             DETERMINE MODEL PRESSURE AT Z:
              IF(JUNITP.EQ.-JUNIT(1))THEN
                  PMODEL=P
              ELSE
                  PMODEL=ZLG4PT(Z,Z0,Z1,Z2,Z3,                          &
     &              PMLOG(I0,-JUNIT(1)),PMLOG(I1,-JUNIT(1)),            &
     &              PMLOG(I2,-JUNIT(1)),PMLOG(I3,-JUNIT(1)),LINIT)
                  LINIT=.FALSE.
              ENDIF

!             DETERMINE MODEL TEMPERATURE AT Z:
              IF(JUNITT.NE.-JUNIT(1))THEN
                  TMODEL=ZLG4PT(Z,Z0,Z1,Z2,Z3,                          &
     &              TMATM(I0,-JUNIT(1)),TMATM(I1,-JUNIT(1)),            &
     &              TMATM(I2,-JUNIT(1)),TMATM(I3,-JUNIT(1)),LINIT)
                  LINIT=.FALSE.
              ENDIF

!             DETERMINE MODEL H2O VOLUME MIXING RATIO AT Z [PPMV]:
              WMOL(1)=ZLG4PT(Z,Z0,Z1,Z2,Z3,                             &
     &          AMOL(I0,-JUNIT(1),1),AMOL(I1,-JUNIT(1),1),              &
     &          AMOL(I2,-JUNIT(1),1),AMOL(I3,-JUNIT(1),1),LINIT)
              LINIT=.FALSE.

!             CONVERT VOLUME MIXING RATIO TO G/M3 (TRATIO IS EXCLUDED
!             FROM BOTH THE ACTUAL & SATURATED DENSITY CALCULATIONS)
              TRATIO=TZERO/TMODEL
              WMOL(1)=WMOL(1)*(PMODEL/PZERO)*AMWT(1)/STDVOL

!             COMPUTE RELATIVE HUMIDITY:
              WMOL(1)=100*WMOL(1)                                       &
     &          /EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
              JUNIT(1)=17
          ENDIF

!         PROFILE SCALING FOR UNIFORMLY MIXED MOLECULAR SPECIES?
          IF(L_UMIX)THEN

!             MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
              DO KSPEC=1,3
                  IF(JUNIT(KSPEC).LE.6)THEN
                      WMOL(KSPEC)=ZLG4PT(Z,Z0,Z1,Z2,Z3,                 &
     &                  AMOL(I0,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I1,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I2,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I3,JUNIT(KSPEC),KSPEC),LINIT)
                      LINIT=.FALSE.

!                     ADJUST CO2 MIXING RATIO BASED ON CARD1A INPUT.
                      IF(KSPEC.EQ.2)WMOL(2)=CO2RAT*WMOL(2)
                      JUNIT(KSPEC)=10
                  ENDIF
              ENDDO
              DO KSPEC=4,7
                  IF(JUNIT(KSPEC).LE.6)THEN
                      WMOL(KSPEC)=S_UMIX(KSPEC)*ZLG4PT(Z,Z0,Z1,Z2,Z3,   &
     &                  AMOL(I0,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I1,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I2,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I3,JUNIT(KSPEC),KSPEC),LINIT)
                      LINIT=.FALSE.
                      JUNIT(KSPEC)=10
                  ENDIF
              ENDDO

!             MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:
              DO KSPEC=8,NMOL
                  IF(JUNIT(KSPEC).LE.6)THEN
                      KSPEC7=KSPEC-7
                      WMOL(KSPEC)=S_UMIX(KSPEC)*ZLG4PT(Z,Z0,Z1,Z2,Z3,   &
     &                  TRAC(I0,KSPEC7),TRAC(I1,KSPEC7),                &
     &                  TRAC(I2,KSPEC7),TRAC(I3,KSPEC7),LINIT)
                      LINIT=.FALSE.
                      JUNIT(KSPEC)=10
                  ENDIF
              ENDDO
          ELSE

!             MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
              DO KSPEC=1,7
                  IF(JUNIT(KSPEC).LE.6)THEN
                      WMOL(KSPEC)=ZLG4PT(Z,Z0,Z1,Z2,Z3,                 &
     &                  AMOL(I0,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I1,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I2,JUNIT(KSPEC),KSPEC),                    &
     &                  AMOL(I3,JUNIT(KSPEC),KSPEC),LINIT)
                      LINIT=.FALSE.

!                     ADJUST CO2 MIXING RATIO BASED ON CARD1A INPUT.
                      IF(KSPEC.EQ.2)WMOL(2)=CO2RAT*WMOL(2)
                      JUNIT(KSPEC)=10
                  ENDIF
              ENDDO

!             MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:
              DO KSPEC=8,NMOL
                  IF(JUNIT(KSPEC).LE.6)THEN
                      KSPEC7=KSPEC-7
                      WMOL(KSPEC)=ZLG4PT(Z,Z0,Z1,Z2,Z3,                 &
     &                  TRAC(I0,KSPEC7),TRAC(I1,KSPEC7),                &
     &                  TRAC(I2,KSPEC7),TRAC(I3,KSPEC7),LINIT)
                      LINIT=.FALSE.
                      JUNIT(KSPEC)=10
                  ENDIF
              ENDDO
          ENDIF

!         MOLECULAR DENSITIES FOR CFC SPECIES:
          IF(JUNITX.LE.6)THEN
              FACTOR=SNGL((Z-Z1)/(Z2-Z1))
              IF(L_XSEC)THEN
                  DO KSPEC=1,NMOLX
                      WMOLX(KSPEC)=ZLG4PT(Z,Z0,Z1,Z2,Z3,                &
     &                  AMOLX(I0,KSPEC),AMOLX(I1,KSPEC),                &
     &                  AMOLX(I2,KSPEC),AMOLX(I3,KSPEC),LINIT)
                      IF(WMOLX(KSPEC).GT.AMOLX(I1,KSPEC) .OR.           &
     &                  WMOLX(KSPEC).LT.AMOLX(I2,KSPEC))WMOLX(KSPEC)    &
     &                  =EXPINT(AMOLX(I1,KSPEC),AMOLX(I2,KSPEC),FACTOR)
                      WMOLX(KSPEC)=S_XSEC(KSPEC)*WMOLX(KSPEC)
                  ENDDO
              ELSE
                  DO KSPEC=1,NMOLX
                      WMOLX(KSPEC)=ZLG4PT(Z,Z0,Z1,Z2,Z3,                &
     &                  AMOLX(I0,KSPEC),AMOLX(I1,KSPEC),                &
     &                  AMOLX(I2,KSPEC),AMOLX(I3,KSPEC),LINIT)
                      IF(WMOLX(KSPEC).GT.AMOLX(I1,KSPEC) .OR.           &
     &                  WMOLX(KSPEC).LT.AMOLX(I2,KSPEC))WMOLX(KSPEC)    &
     &                  =EXPINT(AMOLX(I1,KSPEC),AMOLX(I2,KSPEC),FACTOR)
                  ENDDO
              ENDIF
              LINIT=.FALSE.
          ENDIF

!         MOLECULAR DENSITIES FOR THE Y-SPECIES:
          IF(JUNITY.LE.6)THEN
              IF(L_TRAC)THEN
                  DO KSPEC=1,NMOLY
                      MAP=MAPY_D(KSPEC)
                      IF(MAP.GT.0)WMOLY(KSPEC)=S_TRAC(MAP)              &
     &                  *ZLG4PT(Z,Z0,Z1,Z2,Z3,TRAC(I0,MAP),             &
     &                  TRAC(I1,MAP),TRAC(I2,MAP),TRAC(I3,MAP),LINIT)
                  ENDDO
              ELSE
                  DO KSPEC=1,NMOLY
                      MAP=MAPY_D(KSPEC)
                      IF(MAP.GT.0)WMOLY(KSPEC)                          &
     &                  =ZLG4PT(Z,Z0,Z1,Z2,Z3,TRAC(I0,MAP),             &
     &                  TRAC(I1,MAP),TRAC(I2,MAP),TRAC(I3,MAP),LINIT)
                  ENDDO
              ENDIF
              LINIT=.FALSE.
          ENDIF
          RETURN
      ENDIF

!     USE 2-POINT INTERPOLATION/EXTRAPOLATION NEAR ALTITUDE END POINTS:
      FACTOR=SNGL((Z-ALT(I1))/(ALT(I2)-ALT(I1)))
      IF(JUNITP.LE.6)P=EXP(PMLOG(I1,JUNITP)                             &
     &  +FACTOR*(PMLOG(I2,JUNITP)-PMLOG(I1,JUNITP)))
      IF(JUNITT.LE.6)THEN
          TMODEL=TMATM(I1,JUNITT)
          TMODEL=TMODEL+FACTOR*(TMATM(I2,JUNITT)-TMODEL)
          T=TMODEL+T
      ENDIF

!     IF JUNIT(1) IS NEGATIVE, DETERMINE DEFAULT RELATIVE HUMIDITY:
      IF(JUNIT(1).LT.0)THEN

!         DETERMINE MODEL PRESSURE AT Z:
          IF(JUNITP.EQ.-JUNIT(1))THEN
              PMODEL=P
          ELSE
              PMODEL=PMLOG(I1,-JUNIT(1))
              PMODEL=EXP(PMODEL+FACTOR*(PMLOG(I2,-JUNIT(1))-PMODEL))
          ENDIF

!         DETERMINE MODEL TEMPERATURE AT Z:
          IF(JUNITT.NE.-JUNIT(1))TMODEL=TMATM(I1,-JUNIT(1))             &
     &      +FACTOR*(TMATM(I2,-JUNIT(1))-TMATM(I1,-JUNIT(1)))

!         DETERMINE MODEL H2O VOLUME MIXING RATIO AT Z [PPMV]:
          WMOL(1)=EXPINT(AMOL(I1,-JUNIT(1),1),                          &
     &                   AMOL(I2,-JUNIT(1),1),FACTOR)

!         CONVERT VOLUME MIXING RATIO TO G/M3 (TRATIO IS EXCLUDED
!         FROM BOTH THE ACTUAL & SATURATED DENSITY CALCULATIONS)
          TRATIO=TZERO/TMODEL
          WMOL(1)=WMOL(1)*(PMODEL/PZERO)*AMWT(1)/STDVOL

!         COMPUTE RELATIVE HUMIDITY:
          WMOL(1)=100*WMOL(1)                                           &
     &      /EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
          JUNIT(1)=17
      ENDIF

!     SCALE DEFAULT PROFILES OF UNIFORMLY MIXED MOLECULAR SPECIES?
      IF(L_UMIX)THEN

!         MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
          DO KSPEC=1,3
              IF(JUNIT(KSPEC).LE.6)THEN
                  WMOL(KSPEC)=EXPINT(AMOL(I1,JUNIT(KSPEC),KSPEC),       &
     &                               AMOL(I2,JUNIT(KSPEC),KSPEC),FACTOR)
                  IF(KSPEC.EQ.2)WMOL(2)=CO2RAT*WMOL(2)
                  JUNIT(KSPEC)=10
              ENDIF
          ENDDO
          DO KSPEC=4,7
              IF(JUNIT(KSPEC).LE.6)THEN
                  WMOL(KSPEC)=EXPINT(AMOL(I1,JUNIT(KSPEC),KSPEC),       &
     &              AMOL(I2,JUNIT(KSPEC),KSPEC),FACTOR)*S_UMIX(KSPEC)
                  JUNIT(KSPEC)=10
              ENDIF
          ENDDO

!         MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:
          DO KSPEC=8,NMOL
              IF(JUNIT(KSPEC).LE.6)THEN
                  WMOL(KSPEC)=S_UMIX(KSPEC)                             &
     &              *EXPINT(TRAC(I1,KSPEC-7),TRAC(I2,KSPEC-7),FACTOR)
                  JUNIT(KSPEC)=10
              ENDIF
          ENDDO
      ELSE

!         MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
          DO KSPEC=1,7
              IF(JUNIT(KSPEC).LE.6)THEN
                  WMOL(KSPEC)=EXPINT(AMOL(I1,JUNIT(KSPEC),KSPEC),       &
     &                               AMOL(I2,JUNIT(KSPEC),KSPEC),FACTOR)
                  IF(KSPEC.EQ.2)WMOL(2)=CO2RAT*WMOL(2)
                  JUNIT(KSPEC)=10
              ENDIF
          ENDDO

!         MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:
          DO KSPEC=8,NMOL
              IF(JUNIT(KSPEC).LE.6)THEN
                  WMOL(KSPEC)=EXPINT(TRAC(I1,KSPEC-7),                  &
     &                               TRAC(I2,KSPEC-7),FACTOR)
                  JUNIT(KSPEC)=10
              ENDIF
          ENDDO
      ENDIF

!     MOLECULAR DENSITIES FOR CFC SPECIES:
      IF(JUNITX.LE.6)THEN
          IF(L_XSEC)THEN
              DO KSPEC=1,NMOLX
                  WMOLX(KSPEC)=S_XSEC(KSPEC)                            &
     &              *EXPINT(AMOLX(I1,KSPEC),AMOLX(I2,KSPEC),FACTOR)
              ENDDO
          ELSE
              DO KSPEC=1,NMOLX
                  WMOLX(KSPEC)                                          &
     &              =EXPINT(AMOLX(I1,KSPEC),AMOLX(I2,KSPEC),FACTOR)
              ENDDO
          ENDIF
      ENDIF

!     MOLECULAR DENSITIES FOR THE Y-SPECIES:
      IF(JUNITY.LE.6)THEN
          IF(L_TRAC)THEN
              DO KSPEC=1,NMOLY
                  MAP=MAPY_D(KSPEC)
                  IF(MAP.GT.0)WMOLY(KSPEC)=S_TRAC(MAP)                  &
     &              *EXPINT(TRAC(I1,MAP),TRAC(I2,MAP),FACTOR)
              ENDDO
          ELSE
              DO KSPEC=1,NMOLY
                  MAP=MAPY_D(KSPEC)
                  IF(MAP.GT.0)WMOLY(KSPEC)                              &
     &              =EXPINT(TRAC(I1,MAP),TRAC(I2,MAP),FACTOR)
              ENDDO
          ENDIF
      ENDIF

!     INTERPOLATIONS COMPLETE.
      RETURN
      END

      REAL FUNCTION ZLG4PT(Z,Z1,Z2,Z3,Z4,F1,F2,F3,F4,LINIT)

!     FUNCTION ZLG4PT DOES A 4 POINT LAGRANGE INTERPOLATION FOR DOUBLE
!     PRECISION Z BETWEEN Z2 AND Z3, INCLUSIVE.  THE ABSCISSAE MUST BE
!     MONOTIC (Z1<Z2<Z3<Z4 OR Z1>Z2>Z3>Z4) AND IF F(Z) HAS AN
!     EXTREMUM BETWEEN Z2 AND Z3, THEN A LINEAR INTERPOLATION
!     IS USED INSTEAD.  THE 4-POINT INTERPOLATION FORMULA IS:

!     F(Z)  =  COEF1 P1     COEF2 P2     COEF3 P3     COEF4 P4

!     WHERE

!                 F1              F2              F3              F4
!          P1 = ------ ;   P2 = ------ ;   P4 = ------ ;   P4 = ------
!               DEMON1          DENOM2          DENOM3          DENOM4

!          COEF1  =  (Z - Z2) (Z - Z3) (Z - Z4)

!          COEF2  =  (Z - Z1) (Z - Z3) (Z - Z4)

!          COEF3  =  (Z - Z1) (Z - Z2) (Z - Z4)

!          COEF4  =  (Z - Z1) (Z - Z2) (Z - Z3)

!          DENOM1  =   (Z1 - Z2) (Z1 - Z3) (Z1 - Z4)

!          DENOM2  =   (Z2 - Z1) (Z2 - Z3) (Z2 - Z4)

!          DENOM3  =   (Z3 - Z1) (Z3 - Z2) (Z3 - Z4)

!          DENOM4  =   (Z4 - Z1) (Z4 - Z2) (Z4 - Z3)

!     THE EXTREMA IN F(Z) ARE FOUND BY SOLVING

!           d F(Z)        2
!           ------  =  A Z   -  2 B Z  +  C  =  0
!             dZ

!     WHERE

!        A  =  3 (P1 + P2 + P3 + P4)

!        B  =  (Z2 + Z3 + Z4) P1  +  (Z1 + Z3 + Z4) P2

!           +  (Z1 + Z2 + Z4) P3  +  (Z1 + Z2 + Z3) P4

!        C  =  (Z2 Z3 + Z2 Z4 + Z3 Z4) P1  +  (Z1 Z3 + Z1 Z4 + Z3 Z4) P2

!           +  (Z1 Z2 + Z1 Z4 + Z2 Z4) P3  +  (Z1 Z2 + Z1 Z3 + Z2 Z3) P4

!     THE DISCRIMINANT OF THE QUADRATIC EQUATION FOR Z
!     IS EXPANDED IN QUADRATIC TERMS OF P:

!                           _           _
!          2               \           \
!         B  - 4 A C   =    |  Pi       |       Pj  COEFij
!                          /_          /_
!                           i      j=i and j>i

!     WHERE

!                     1           2           2             2
!          COEFii  =  -  { (Zj-Zk)   + (Zj-Zl)   +   (Zk-Zl)  }
!                     2

!                              2
!          COEFij  =  2 (Zk-Zl)  +  (Zi-Zk) (Zj-Zl)  +  (Zi-Zl) (Zj-Zk)
      IMPLICIT NONE

!     PARAMETERS:
!       SMALL    SMALL NUMBER TOLERANCE
      REAL SMALL
      PARAMETER(SMALL=1.E-30)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       Z        ABSCISSA OF DESIRED ORDINATE.
!       Z1       FIRST ABSCISSA
!       Z2       SECOND ABSCISSA
!       Z3       THIRD ABSCISSA
!       Z4       FOURTH ABSCISSA
!       F1       FIRST ORDINATE
!       F2       SECOND ORDINATE
!       F3       THIRD ORDINATE
!       F4       FOURTH ORDINATE
!       LINIT    FLAG, FALSE IF ABSCISSAE ARE UNCHANGED.
      DOUBLE PRECISION Z,Z1,Z2,Z3,Z4
      REAL F1,F2,F3,F4
      LOGICAL LINIT

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     SAVED VARIABLES:
      REAL BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  &
     &  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, &
     &  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4
      SAVE BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  &
     &  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, &
     &  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4

!     LOCAL VARIABLES:
      DOUBLE PRECISION ZDIFF1,ZDIFF2,ZDIFF3,ZDIFF4,                     &
     &  ZDIF12,ZDIF13,ZDIF14,ZDIF23,ZDIF24,ZDIF34,                      &
     &  ZD12SQ,ZD13SQ,ZD14SQ,ZD23SQ,ZD24SQ,ZD34SQ,STORE,ZEXTRM
      REAL A,B,DSCRIM,ROOT,P1,P2,P3,P4

!     NEW ABSCISSAE CHECK
      IF(LINIT)THEN

!         CHECK INPUTS
          IF(Z2.EQ.Z3)GOTO 10
          IF(Z2.LT.Z3)THEN
             IF(Z1.GE.Z2 .OR. Z2.GT.Z .OR. Z.GT.Z3 .OR. Z3.GE.Z4)GOTO 10
          ELSE
             IF(Z1.LE.Z2 .OR. Z2.LT.Z .OR. Z.LT.Z3 .OR. Z3.LE.Z4)GOTO 10
          ENDIF

!         DEFINE CONSTANTS DEPENDENT ON Z, Z1, Z2, Z3 AND Z4 ONLY.
          ZDIFF1=Z-Z1
          ZDIFF2=Z-Z2
          ZDIFF3=Z-Z3
          ZDIFF4=Z-Z4
          COEF1=SNGL(ZDIFF2*ZDIFF3*ZDIFF4)
          COEF2=SNGL(ZDIFF1*ZDIFF3*ZDIFF4)
          COEF3=SNGL(ZDIFF1*ZDIFF2*ZDIFF4)
          COEF4=SNGL(ZDIFF1*ZDIFF2*ZDIFF3)
          ZDIF12=Z1-Z2
          ZDIF13=Z1-Z3
          ZDIF14=Z1-Z4
          ZDIF23=Z2-Z3
          ZDIF24=Z2-Z4
          ZDIF34=Z3-Z4
          DENOM1=SNGL( ZDIF12*ZDIF13*ZDIF14)
          DENOM2=SNGL(-ZDIF12*ZDIF23*ZDIF24)
          DENOM3=SNGL( ZDIF13*ZDIF23*ZDIF34)
          DENOM4=SNGL(-ZDIF14*ZDIF24*ZDIF34)

!         CONSTANTS USED IN EXTREMUM CHECK
          BCOEF1=SNGL(Z2+Z3+Z4)
          BCOEF2=SNGL(Z3+Z4+Z1)
          BCOEF3=SNGL(Z4+Z1+Z2)
          BCOEF4=SNGL(Z1+Z2+Z3)
          ZD12SQ=ZDIF12**2
          ZD13SQ=ZDIF13**2
          ZD14SQ=ZDIF14**2
          ZD23SQ=ZDIF23**2
          ZD24SQ=ZDIF24**2
          ZD34SQ=ZDIF34**2
          COEF11=SNGL((ZD23SQ+ZD24SQ+ZD34SQ)/2)
          COEF22=SNGL((ZD13SQ+ZD14SQ+ZD34SQ)/2)
          COEF33=SNGL((ZD12SQ+ZD14SQ+ZD24SQ)/2)
          COEF44=SNGL((ZD12SQ+ZD13SQ+ZD23SQ)/2)
          STORE=ZDIF13*ZDIF24+ZDIF14*ZDIF23
          COEF12=SNGL(2*ZD34SQ+STORE)
          COEF34=SNGL(2*ZD12SQ+STORE)
          STORE=ZDIF12*ZDIF34-ZDIF14*ZDIF23
          COEF13=SNGL(2*ZD24SQ+STORE)
          COEF24=SNGL(2*ZD13SQ+STORE)
          STORE=ZDIF12*ZDIF34+ZDIF13*ZDIF24
          COEF14=SNGL(2*ZD23SQ-STORE)
          COEF23=SNGL(2*ZD14SQ-STORE)
          CCOEF1=SNGL(Z2*Z3+Z2*Z4+Z3*Z4)
          CCOEF2=SNGL(Z1*Z3+Z1*Z4+Z3*Z4)
          CCOEF3=SNGL(Z1*Z2+Z1*Z4+Z2*Z4)
          CCOEF4=SNGL(Z1*Z2+Z1*Z3+Z2*Z3)

!         LINEAR INTERPOLATION FACTOR
          FACTOR=SNGL((Z-Z2)/(Z3-Z2))
      ENDIF

!     DETERMINE LOCATION OF EXTREMA
      P1=F1/DENOM1
      P2=F2/DENOM2
      P3=F3/DENOM3
      P4=F4/DENOM4
      A=3*(P1+P2+P3+P4)
      B=P1*BCOEF1+P2*BCOEF2+P3*BCOEF3+P4*BCOEF4
      IF(ABS(A).LT.SMALL)THEN
          IF(ABS(B).GT.SMALL)THEN

!             ONE EXTREMUM
              ZEXTRM                                                    &
     &          =DBLE((CCOEF1*P1+CCOEF2*P2+CCOEF3*P3+CCOEF4*P4)/(2*B))
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  ZLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
          ENDIF
      ELSE
          DSCRIM=P1*(P1*COEF11+P2*COEF12+P3*COEF13+P4*COEF14)           &
     &          +P2*(P2*COEF22+P3*COEF23+P4*COEF24)                     &
     &          +P3*(P3*COEF33+P4*COEF34)+P4**2*COEF44
          IF(DSCRIM.GE.0.)THEN

!             TWO EXTREMUM
              ROOT=SQRT(DSCRIM)
              ZEXTRM=DBLE((B+ROOT)/A)
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  ZLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
              ZEXTRM=DBLE((B-ROOT)/A)
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  ZLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
          ENDIF
      ENDIF

!     NO EXTREMA BETWEEN Z2 AND Z3.  USE 4-POINT LAGRANGE INTERPOLATION.
      ZLG4PT=COEF1*P1+COEF2*P2+COEF3*P3+COEF4*P4
      RETURN

!     INCORRECT INPUT ERROR.
   10 CONTINUE
      WRITE(IPR,'(/A,/(18X,A,F12.4))')                                  &
     &  ' Error in ZLG4PT:  ABSCISSAE OUT OF ORDER.',                   &
     &  ' Z1 =',Z1,' Z2 =',Z2,' Z  =',Z,' Z3 =',Z3,' Z4 =',Z4
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP ' Error in ZLG4PT:  ABSCISSAE OUT OF ORDER.'
      END
