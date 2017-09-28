      SUBROUTINE DFLTP(L,GNDALT)

!     THIS SUBROUTINE INTERPOLATES PROFILES FROM THE 6 BUILT-IN MODEL
!     ATMOSPHERES TO PRESSURE.  THE JUNIT INDICES FROM /CARD1B/
!     INDICATE WHICH PROFILE SHOULD BE USED FOR EACH INTERPOLATION:

!                   JUNIT     MODEL ATMOSPHERE
!                     1          TROPICAL
!                     2          MID-LATITUDE SUMMER
!                     3          MID-LATITUDE WINTER
!                     4          HIGH-LAT SUMMER
!                     5          HIGH-LAT WINTER
!                     6          US STANDARD

!     IF JUNITP IS 10 OR MORE, THE ALTITUDE IS DETERMINED
!     FROM THE HYDROSTATIC EQUATION.  FOR THE FIRST LEVEL,
!     THE ALTITUDE IS SET TO GNDALT.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       L        LAYER INDEX.
!       GNDALT   GROUND ALTITUDE [KM].
      INTEGER L
      DOUBLE PRECISION GNDALT

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
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
!       AMOLX    VERTICAL PROFILE OF CROSS-SECTION SPECIES [PPMV].
      REAL AMOLX
      COMMON/MLATMX/AMOLX(NLEVEL,NMOLX)

!     /TRACE/
!       TRAC     TRACE MOLECULAR SPECIES DEFAULT PROFILES [PPMV].
      REAL TRAC
      COMMON/TRACE/TRAC(NLEVEL,21)

!     /CO2MIX/
      REAL CO2RAT
      COMMON/CO2MIX/CO2RAT

!     /CD1BXY/
!       JUNITX   INPUT UNIT FLAG FOR CROSS-SECTION MOLECULES.
!       WMOLX    X-SECTION MOLECULAR DENSITY [PPMV AFTER JUBRAN CALL].
!       JUNITY   INPUT UNIT FLAG FOR AUXILIARY (Y) MOLECULES.
!       WMOLY    AUXILIARY MOLECULAR DENSITY [PPMV AFTER JUBRAN CALL].
      INTEGER JUNITX,JUNITY
      REAL WMOLX,WMOLY
      COMMON/CD1BXY/JUNITX,WMOLX(NMOLX),JUNITY,WMOLY(MMOLY)

!     /M_PTWO/
!       PPROF    PRESSURE PROFILE [MB].
!       TPROF    TEMPERATURE PROFILE [K].
!       WH2O     H2O VOLUME MIXING RATIO PROFILE [PPMV].
!       WO3      O3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL PPROF,TPROF,WH2O,WO3
      COMMON/M_PTWO/PPROF(LAYDIM),TPROF(LAYDIM),WH2O(LAYDIM),WO3(LAYDIM)

!     /MPROF/
!       ZM       PROFILE LEVEL ALTITUDES [KM].
!       PM       PROFILE LEVEL PRESSURES [MBAR].
!       TM       PROFILE LEVEL TEMPERATURES [K].
!       RFNDX    PROFILE LEVEL REFRACTIVITIES.
!       LRHSET   FLAG, .TRUE. IF RELATIVE HUMIDITY IS NOT TO BE SCALED.
      DOUBLE PRECISION ZM
      REAL PM,TM,RFNDX
      LOGICAL LRHSET
      COMMON/MPROF/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    &
     &  RFNDX(LAYDIM),LRHSET(LAYDIM)

!     /CONSTN/
!       AMWT     MOLECULAR WEIGHTS [GM/MOLE].
!                  1  H2O   2  CO2   3  O3    4  N2O   5  CO    6  CH4
!                  7  O2    8  NO    9  SO2   10 NO2   11 NH3   12 HNO3
      REAL AMWT
      COMMON/CONSTN/AMWT(12)

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
      EXTERNAL MLATMB,XMLATM,DEVCBD

!     FUNCTIONS:
!       HYSTAT   INTEGRATES HYDROSTATIC EQUATION TO FIND ALTITUDE.
      REAL FLG4PT
      DOUBLE PRECISION HYSTAT

!     LOCAL VARIABLES:
!       KSPEC    SPECIES LOOP INDEX.
!       KSPEC7   KSPEC LESS 7.
!       TMODEL   MODEL ATMOSPHERE TEMPERATURE [K].
!       TRATIO   RATIO OF MODEL TEMPERATURE TO STANDARD TEMPERATURE.
!       MAP      MAPPING FROM ACTIVE Y-SPECIES TO DEFAULT PROFILE.
      INTEGER I0,I1,I2,I3,J0,J1,J2,J3,KSPEC,KSPEC7,LM1,MAP
      LOGICAL LINIT
      REAL PLOG,PBNDLO,PBNDHI,P0,P1,P2,P3,TMODEL,TRATIO

!     DEFINE LOG(PRESSURE/MBAR)
      IF(PPROF(L).LE.0.)THEN
          WRITE(IPR,'(2A)')' Error:  Pressure must be speified',        &
     &      ' for each atmospheric level when MODEL is 8.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' Input pressure must be positive when MODEL is 8.'
      ELSEIF(L.GT.1)THEN
          LM1=L-1
          IF(PPROF(L).GT.PPROF(LM1))THEN
              WRITE(IPR,'(2A)')                                         &
     &          ' Error:  Pressure inversions are not allowed.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' Pressure inversions are not allowed.'
          ENDIF
      ENDIF
      PLOG=LOG(PPROF(L))

!     TEMPERATURE:
      LINIT=.TRUE.
      IF(JUNITT.LE.6)THEN

!         TEMPERATURE DETERMINED FROM MODEL ATMOSPHERE:
          PBNDLO=PMLOG(2,JUNITT)
          PBNDHI=PMLOG(NLEVM1,JUNITT)
          IF(PLOG.GE.PBNDLO)THEN

!             INPUT PRESSURE IS ABOVE SECOND ATMOSPHERIC LEVEL VALUE:
              TMODEL=TMATM(2,JUNITT)
              TMODEL=TMODEL+(TMATM(1,JUNITT)-TMODEL)                    &
     &          *(PLOG-PBNDLO)/(PMLOG(1,JUNITT)-PBNDLO)
          ELSEIF(PLOG.LE.PBNDHI)THEN

!             INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!             SECOND-TO-LAST ATMOSPHERIC LEVEL:
              TMODEL=TMATM(NLEVM1,JUNITT)
              TMODEL=TMODEL+(TMATM(NLEVEL,JUNITT)-TMODEL)               &
     &          *(PLOG-PBNDHI)/(PMLOG(NLEVEL,JUNITT)-PBNDHI)
          ELSE

!             INTERMEDIATE PRESSURE:
              DO I2=2,NLEVM1
                  IF(PLOG.GE.PMLOG(I2,JUNITT))GOTO 10
              ENDDO
   10         CONTINUE
              I0=I2-2
              I1=I2-1
              I3=I2+1
              P0=PMLOG(I0,JUNITT)
              P1=PMLOG(I1,JUNITT)
              P2=PMLOG(I2,JUNITT)
              P3=PMLOG(I3,JUNITT)
              TMODEL=FLG4PT(PLOG,P0,P1,P2,P3,                           &
     &          TMATM(I0,JUNITT),TMATM(I1,JUNITT),                      &
     &          TMATM(I2,JUNITT),TMATM(I3,JUNITT),LINIT)
              IF(JUNITP.EQ.JUNITT)LINIT=.FALSE.
          ENDIF
          TPROF(L)=TMODEL+TPROF(L)
      ENDIF

!     ALTITUDE:
      IF(JUNITP.GE.1 .AND. JUNITP.LE.6)THEN

!         ALTITUDE DETERMINED FROM MODEL ATMOSPHERE:
          PBNDLO=PMLOG(2,JUNITP)
          PBNDHI=PMLOG(NLEVM1,JUNITP)
          IF(PLOG.GE.PBNDLO)THEN

!             INPUT PRESSURE IS ABOVE PRESSURE OF
!             SECOND ATMOSPHERIC LEVEL:
              ZM(L)=ALT(2)+(ALT(1)-ALT(2))                              &
     &          *DBLE((PLOG-PBNDLO)/(PMLOG(1,JUNITP)-PBNDLO))
          ELSEIF(PLOG.LE.PBNDHI)THEN

!             INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!             SECOND-TO-LAST ATMOSPHERIC LEVEL:
              ZM(L)=ALT(NLEVM1)+(ALT(NLEVEL)-ALT(NLEVM1))               &
     &          *DBLE((PLOG-PBNDHI)/(PMLOG(NLEVEL,JUNITP)-PBNDHI))
          ELSE
              IF(LINIT)THEN

!                 P AND T FROM DIFFERENT MODEL ATMOSPHERES:
                  DO I2=2,NLEVM1
                      IF(PLOG.GE.PMLOG(I2,JUNITP))GOTO 20
                  ENDDO
   20             CONTINUE
                  I0=I2-2
                  I1=I2-1
                  I3=I2+1
                  P0=PMLOG(I0,JUNITP)
                  P1=PMLOG(I1,JUNITP)
                  P2=PMLOG(I2,JUNITP)
                  P3=PMLOG(I3,JUNITP)
              ENDIF
              ZM(L)=DBLE(FLG4PT(PLOG,P0,P1,P2,P3,SNGL(ALT(I0)),         &
     &          SNGL(ALT(I1)),SNGL(ALT(I2)),SNGL(ALT(I3)),LINIT))
              LINIT=.FALSE.
          ENDIF
      ELSEIF(L.EQ.1)THEN

!         SET ALTITUDE AT PPROF(1) TO INPUT GNDALT:
          ZM(1)=GNDALT
      ELSE

!         INTEGRATE THE HYDROSTATIC EQUATION ASSUMING
!         TEMPERATURE VARIES LINEARLY WITH LOG(PRESSURE):
          ZM(L)=HYSTAT(ZM(LM1),TPROF(LM1),PPROF(LM1),TPROF(L),PPROF(L))
      ENDIF

!     IF JUNIT(1) IS NEGATIVE, DETERMINE DEFAULT RELATIVE HUMIDITY:
      IF(JUNIT(1).LT.0)THEN

!         COMPUTE MODEL H2O VOLUME MIXING RATIO AND TEMPERATURE [PPMV]:
          PBNDLO=PMLOG(2,-JUNIT(1))
          PBNDHI=PMLOG(NLEVM1,-JUNIT(1))
          IF(PLOG.GE.PBNDLO)THEN

!             INPUT PRESSURE IS ABOVE SECOND ATMOSPHERIC LEVEL VALUE:
              WMOL(1)=AMOL(2,-JUNIT(1),1)
              WMOL(1)=WMOL(1)*(AMOL(1,-JUNIT(1),1)/WMOL(1))             &
     &          **((PLOG-PBNDLO)/(PMLOG(1,-JUNIT(1))-PBNDLO))
              IF(JUNITT.NE.-JUNIT(1))THEN

!                 TEMPERATURE DETERMINED FROM SAME MODEL ATMOSPHERE:
                  TMODEL=TMATM(2,-JUNIT(1))
                  TMODEL=TMODEL+(TMATM(1,-JUNIT(1))-TMODEL)             &
     &              *(PLOG-PBNDLO)/(PMLOG(1,-JUNIT(1))-PBNDLO)
              ENDIF
          ELSEIF(PLOG.LE.PBNDHI)THEN

!             INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!             SECOND-TO-LAST ATMOSPHERIC LEVEL:
              WMOL(1)=AMOL(NLEVM1,-JUNIT(1),1)
              WMOL(1)=WMOL(1)*(AMOL(NLEVEL,-JUNIT(1),1)/WMOL(1))        &
     &          **((PLOG-PBNDHI)/(PMLOG(NLEVEL,-JUNIT(1))-PBNDHI))
              IF(JUNITT.NE.-JUNIT(1))THEN

!                 TEMPERATURE DETERMINED FROM SAME MODEL ATMOSPHERE:
                  TMODEL=TMATM(NLEVM1,-JUNIT(1))
                  TMODEL=TMODEL+(TMATM(NLEVEL,-JUNIT(1))-TMODEL)        &
     &              *(PLOG-PBNDHI)/(PMLOG(NLEVEL,-JUNIT(1))-PBNDHI)
              ENDIF
          ELSEIF(JUNITP.EQ.-JUNIT(1))THEN

!             P AND MOL FROM THE SAME MODEL ATMOSPHERE:
              WMOL(1)=FLG4PT(PLOG,P0,P1,P2,P3,                          &
     &          AMOL(I0,JUNITP,1),AMOL(I1,JUNITP,1),                    &
     &          AMOL(I2,JUNITP,1),AMOL(I3,JUNITP,1),LINIT)
              LINIT=.FALSE.

!             TEMPERATURE DETERMINED FROM SAME MODEL ATMOSPHERE:
              IF(JUNITT.NE.-JUNIT(1))TMODEL=FLG4PT(PLOG,P0,P1,P2,P3,    &
     &          TMATM(I0,-JUNIT(1)),TMATM(I1,-JUNIT(1)),                &
     &          TMATM(I2,-JUNIT(1)),TMATM(I3,-JUNIT(1)),LINIT)
          ELSE

!             P AND MOL FROM DIFFERENT MODEL ATMOSPHERES:
              DO J2=2,NLEVM1
                  IF(PLOG.GE.PMLOG(J2,-JUNIT(1)))GOTO 30
              ENDDO
   30         CONTINUE
              J0=J2-2
              J1=J2-1
              J3=J2+1
              LINIT=.TRUE.
              WMOL(1)=FLG4PT(PLOG,                                      &
     &          PMLOG(J0,-JUNIT(1)),PMLOG(J1,-JUNIT(1)),                &
     &          PMLOG(J2,-JUNIT(1)),PMLOG(J3,-JUNIT(1)),                &
     &          AMOL(J0,-JUNIT(1),1),AMOL(J1,-JUNIT(1),1),              &
     &          AMOL(J2,-JUNIT(1),1),AMOL(J3,-JUNIT(1),1),LINIT)
              IF(JUNITT.NE.-JUNIT(1))THEN

!                 TEMPERATURE DETERMINED FROM SAME MODEL ATMOSPHERE:
                  LINIT=.FALSE.
                  TMODEL=FLG4PT(PLOG,                                   &
     &              PMLOG(J0,-JUNIT(1)),PMLOG(J1,-JUNIT(1)),            &
     &              PMLOG(J2,-JUNIT(1)),PMLOG(J3,-JUNIT(1)),            &
     &              TMATM(J0,-JUNIT(1)),TMATM(J1,-JUNIT(1)),            &
     &              TMATM(J2,-JUNIT(1)),TMATM(J3,-JUNIT(1)),LINIT)
                  LINIT=.TRUE.
              ENDIF
          ENDIF

!         CONVERT VOLUME MIXING RATIO TO G/M3 (TRATIO IS EXCLUDED
!         FROM BOTH THE ACTUAL & SATURATED DENSITY CALCULATIONS):
          TRATIO=TZERO/TMODEL
          WMOL(1)=WMOL(1)*(PPROF(L)/PZERO)*AMWT(1)/STDVOL

!         COMPUTE RELATIVE HUMIDITY:
          WMOL(1)=100*WMOL(1)                                           &
     &      /EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
          JUNIT(1)=17
      ENDIF

!     MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
      DO KSPEC=1,7
          IF(JUNIT(KSPEC).LE.6)THEN
              PBNDLO=PMLOG(2,JUNIT(KSPEC))
              PBNDHI=PMLOG(NLEVM1,JUNIT(KSPEC))
              IF(PLOG.GE.PBNDLO)THEN

!                 INPUT PRESSURE IS ABOVE PRESSURE OF
!                 SECOND ATMOSPHERIC LEVEL:
                  WMOL(KSPEC)=AMOL(2,JUNIT(KSPEC),KSPEC)
                  WMOL(KSPEC)=WMOL(KSPEC)*                              &
     &              (AMOL(1,JUNIT(KSPEC),KSPEC)/WMOL(KSPEC))**          &
     &              ((PLOG-PBNDLO)/(PMLOG(1,JUNIT(KSPEC))-PBNDLO))
              ELSEIF(PLOG.LE.PBNDHI)THEN

!                 INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!                 SECOND-TO-LAST ATMOSPHERIC LEVEL:
                  WMOL(KSPEC)=AMOL(NLEVM1,JUNIT(KSPEC),KSPEC)
                  WMOL(KSPEC)=WMOL(KSPEC)*                              &
     &              (AMOL(NLEVEL,JUNIT(KSPEC),KSPEC)/WMOL(KSPEC))**     &
     &              ((PLOG-PBNDHI)/(PMLOG(NLEVEL,JUNIT(KSPEC))-PBNDHI))
              ELSEIF(JUNITP.EQ.JUNIT(KSPEC))THEN

!                 P AND MOL FROM THE SAME MODEL ATMOSPHERE:
                  WMOL(KSPEC)=FLG4PT(PLOG,P0,P1,P2,P3,                  &
     &              AMOL(I0,JUNITP,KSPEC),AMOL(I1,JUNITP,KSPEC),        &
     &              AMOL(I2,JUNITP,KSPEC),AMOL(I3,JUNITP,KSPEC),LINIT)
                  LINIT=.FALSE.
              ELSE

!                 P AND MOL FROM DIFFERENT MODEL ATMOSPHERES:
                  DO J2=2,NLEVM1
                      IF(PLOG.GE.PMLOG(J2,JUNIT(KSPEC)))GOTO 40
                  ENDDO
   40             CONTINUE
                  J0=J2-2
                  J1=J2-1
                  J3=J2+1
                  LINIT=.TRUE.
                  WMOL(KSPEC)=FLG4PT(PLOG,                              &
     &              PMLOG(J0,JUNIT(KSPEC)),PMLOG(J1,JUNIT(KSPEC)),      &
     &              PMLOG(J2,JUNIT(KSPEC)),PMLOG(J3,JUNIT(KSPEC)),      &
     &              AMOL(J0,JUNIT(KSPEC),KSPEC),                        &
     &              AMOL(J1,JUNIT(KSPEC),KSPEC),                        &
     &              AMOL(J2,JUNIT(KSPEC),KSPEC),                        &
     &              AMOL(J3,JUNIT(KSPEC),KSPEC),LINIT)
              ENDIF
              JUNIT(KSPEC)=10

!             ADJUST CO2 MIXING RATIO BASED ON CARD1A INPUT.
              IF(KSPEC.EQ.2)THEN
                  WMOL(2)=CO2RAT*WMOL(2)
              ELSEIF(L_UMIX .AND. KSPEC.GE.4)THEN
                  WMOL(KSPEC)=S_UMIX(KSPEC)*WMOL(KSPEC)
              ENDIF
          ENDIF
      ENDDO

!     MOLECULES CONSTANT (IN ALTITUDE) WITH MODEL ATMOSPHERE:
      DO KSPEC=8,NMOL
          IF(JUNIT(KSPEC).LE.6)THEN
              KSPEC7=KSPEC-7
              PBNDLO=PMLOG(2,JUNIT(KSPEC))
              PBNDHI=PMLOG(NLEVM1,JUNIT(KSPEC))
              IF(PLOG.GE.PBNDLO)THEN

!                 INPUT PRESSURE IS ABOVE PRESSURE OF
!                 SECOND ATMOSPHERIC LEVEL:
                  WMOL(KSPEC)=TRAC(2,KSPEC7)
                  WMOL(KSPEC)=WMOL(KSPEC)*(TRAC(1,KSPEC7)/WMOL(KSPEC))  &
     &              **((PLOG-PBNDLO)/(PMLOG(1,JUNIT(KSPEC))-PBNDLO))
              ELSEIF(PLOG.LE.PBNDHI)THEN

!                 INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!                 SECOND-TO-LAST ATMOSPHERIC LEVEL:
                  WMOL(KSPEC)=TRAC(NLEVM1,KSPEC7)
                  WMOL(KSPEC)=WMOL(KSPEC)*                              &
     &              (TRAC(NLEVEL,KSPEC7)/WMOL(KSPEC))**                 &
     &              ((PLOG-PBNDHI)/(PMLOG(NLEVEL,JUNIT(KSPEC))-PBNDHI))
              ELSEIF(JUNITP.EQ.JUNIT(KSPEC))THEN

!                 P AND MOL FROM THE SAME MODEL ATMOSPHERE:
                  WMOL(KSPEC)=FLG4PT(PLOG,P0,P1,P2,P3,                  &
     &              TRAC(I0,KSPEC7),TRAC(I1,KSPEC7),                    &
     &              TRAC(I2,KSPEC7),TRAC(I3,KSPEC7),LINIT)
                  LINIT=.FALSE.
              ELSE

!                 P AND MOL FROM DIFFERENT MODEL ATMOSPHERES:
                  DO J2=2,NLEVM1
                      IF(PLOG.GE.PMLOG(J2,JUNIT(KSPEC)))GOTO 50
                  ENDDO
   50             CONTINUE
                  J0=J2-2
                  J1=J2-1
                  J3=J2+1
                  LINIT=.TRUE.
                  WMOL(KSPEC)=FLG4PT(PLOG,                              &
     &              PMLOG(J0,JUNIT(KSPEC)),PMLOG(J1,JUNIT(KSPEC)),      &
     &              PMLOG(J2,JUNIT(KSPEC)),PMLOG(J3,JUNIT(KSPEC)),      &
     &              TRAC(J0,KSPEC7),TRAC(J1,KSPEC7),                    &
     &              TRAC(J2,KSPEC7),TRAC(J3,KSPEC7),LINIT)
              ENDIF
              JUNIT(KSPEC)=10
              IF(L_UMIX)WMOL(KSPEC)=S_UMIX(KSPEC)*WMOL(KSPEC)
          ENDIF
      ENDDO

!     MOLECULAR DENSITIES FOR CFC SPECIES:
      IF(JUNITX.LE.6)THEN
          LINIT=JUNITP.NE.JUNITX
          DO KSPEC=1,NMOLX
              PBNDLO=PMLOG(2,JUNITX)
              PBNDHI=PMLOG(NLEVM1,JUNITX)
              IF(PLOG.GE.PBNDLO)THEN

!                 INPUT PRESSURE IS ABOVE PRESSURE OF
!                 SECOND ATMOSPHERIC LEVEL:
                  WMOLX(KSPEC)=AMOLX(2,KSPEC)
                  WMOLX(KSPEC)=WMOLX(KSPEC)                             &
     &              *(AMOLX(1,KSPEC)/WMOLX(KSPEC))                      &
     &              **((PLOG-PBNDLO)/(PMLOG(1,JUNITX)-PBNDLO))
              ELSEIF(PLOG.LE.PBNDHI)THEN

!                 INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!                 SECOND-TO-LAST ATMOSPHERIC LEVEL:
                  WMOLX(KSPEC)=AMOLX(NLEVM1,KSPEC)
                  WMOLX(KSPEC)=WMOLX(KSPEC)                             &
     &              *(AMOLX(NLEVEL,KSPEC)/WMOLX(KSPEC))                 &
     &              **((PLOG-PBNDHI)/(PMLOG(NLEVEL,JUNITX)-PBNDHI))
              ELSE
                  IF(LINIT)THEN

!                     P AND MOL FROM THE SAME MODEL ATMOSPHERE:
                      DO I2=2,NLEVM1
                          IF(PLOG.GE.PMLOG(I2,JUNITX))GOTO 60
                      ENDDO
   60                 CONTINUE
                      I0=I2-2
                      I1=I2-1
                      I3=I2+1
                      P0=PMLOG(I0,JUNITX)
                      P1=PMLOG(I1,JUNITX)
                      P2=PMLOG(I2,JUNITX)
                      P3=PMLOG(I3,JUNITX)
                      LINIT=.TRUE.
                  ENDIF
                  WMOLX(KSPEC)=FLG4PT(PLOG,P0,P1,P2,P3,                 &
     &              AMOLX(I0,KSPEC),AMOLX(I1,KSPEC),                    &
     &              AMOLX(I2,KSPEC),AMOLX(I3,KSPEC),LINIT)
                  LINIT=.FALSE.
              ENDIF
              IF(L_XSEC)WMOLX(KSPEC)=S_XSEC(KSPEC)*WMOLX(KSPEC)
          ENDDO
      ENDIF

!     MOLECULAR DENSITIES FOR THE Y-SPECIES:
      IF(JUNITY.LE.6)THEN
          LINIT=.TRUE.
          DO KSPEC=1,NMOLY
              MAP=MAPY_D(KSPEC)
              IF(MAP.GT.0)THEN
                  PBNDLO=PMLOG(2,JUNITY)
                  PBNDHI=PMLOG(NLEVM1,JUNITY)
                  IF(PLOG.GE.PBNDLO)THEN

!                     INPUT PRESSURE IS ABOVE PRESSURE OF
!                     SECOND ATMOSPHERIC LEVEL:
                      WMOLY(KSPEC)=TRAC(2,MAP)*(TRAC(1,MAP)/TRAC(2,MAP))&
     &                  **((PLOG-PBNDLO)/(PMLOG(1,JUNITY)-PBNDLO))
                  ELSEIF(PLOG.LE.PBNDHI)THEN

!                     INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!                     SECOND-TO-LAST ATMOSPHERIC LEVEL:
                      WMOLY(KSPEC)=TRAC(NLEVM1,MAP)
                      WMOLY(KSPEC)=WMOLY(KSPEC)                         &
     &                  *(TRAC(NLEVEL,MAP)/WMOLY(KSPEC))                &
     &                  **((PLOG-PBNDHI)/(PMLOG(NLEVEL,JUNITY)-PBNDHI))
                  ELSEIF(.NOT.LINIT)THEN

!                     P AND MOL FROM THE SAME MODEL ATMOSPHERE:
                      WMOLY(KSPEC)=FLG4PT(PLOG,P0,P1,P2,P3,             &
     &                  TRAC(I0,MAP),TRAC(I1,MAP),                      &
     &                  TRAC(I2,MAP),TRAC(I3,MAP),LINIT)
                  ELSE

!                     P AND MOL FROM DIFFERENT MODEL ATMOSPHERES:
                      DO J2=2,NLEVM1
                          IF(PLOG.GE.PMLOG(J2,JUNITY))GOTO 70
                      ENDDO
   70                 CONTINUE
                      J0=J2-2
                      J1=J2-1
                      J3=J2+1
                      P0=PMLOG(J0,JUNITY)
                      P1=PMLOG(J1,JUNITY)
                      P2=PMLOG(J2,JUNITY)
                      P3=PMLOG(J3,JUNITY)
                      WMOLY(KSPEC)=FLG4PT(PLOG,P0,P1,P2,P3,TRAC(J0,MAP),&
     &                  TRAC(J1,MAP),TRAC(J2,MAP),TRAC(J3,MAP),LINIT)
                      LINIT=.FALSE.
                  ENDIF
                  IF(L_TRAC)WMOLY(KSPEC)=S_TRAC(MAP)*WMOLY(KSPEC)
              ENDIF
          ENDDO
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION HYSTAT(ZREF,TREF,PREF,T,P)

!     INTEGRATES HYDROSTATIC EQUATION ASSUMING TEMPERATURE
!     VARIES LINEARLY WITH LOG(PRESSURE):

!     INPUT ARGUMENTS:
!       ZREF     REFERENCE ALTITUDE [KM].
!       TREF     REFERENCE TEMPERATURE [K].
!       PREF     REFERENCE PRESSURE [MB].
!       T        TEMPERATURE AT UNKNOWN ALTITUDE [K].
!       P        PRESSURE AT UNKNOWN ALTITUDE [MB].
      DOUBLE PRECISION ZREF
      REAL TREF,PREF,T,P

!     PARAMETERS:
!       GRAV     GRAVITATIONAL CONSTANT OF EARTH [KM3/SEC2].
!       GASCON   GAS CONSTANT [GM KM2 / SEC2 MOL K].
!       HSTAT    CONSTANT USED IN HYDROSTATIC EQUATION [KM K].
      INCLUDE 'PARAMS.h'
      DOUBLE PRECISION GRAV,GASCON,HSTAT
      PARAMETER(GRAV=398600.4418D0,GASCON=.00831441D0,                  &
     &  HSTAT=2*AIRMWT*GRAV/GASCON)

!     /CARD3/
!       H1ALT    OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    FINAL (TANGENT FOR LIMB PATH) ALTITUDE [KM].
!       OBSZEN   OBSERVER ZENITH ANGLE (H1ALT TO H2ALT) [DEG].
!       HRANGE   DISTANCE FROM H1ALT TO H2ALT [KM].
!       BETA     EARTH CENTER ANGLE BETWEEN H1ALT AND H2ALT [DEG].
!       REARTH   RADIUS OF THE EARTH [KM].
!       HMIN     PATH MINIMUM ALTITUDE [KM].
!       HMAX     PATH MAXIMUM ALTITUDE [KM].
!       CKRANG   MAXIMUM PATH RANGE FOR K-DISTRIBUTION OUTPUT
!                (=0. FOR TOTAL PATH ONLY; <0. FOR ALL RANGES).
!       BCKZEN   ZENITH ANGLE FOR BACKWARD (H2ALT TO H1ALT) PATH [DEG].
!       REARTH   EARTH RADIUS [KM].
!       ANGLEM   LUNAR PHASE ANGLE [0 TO 180 DEG].
!       LENN     PATH LENGTH SWITCH (0=SHORT, 1=LONG).
!       IDAY     DAY OF YEAR [0-366, DEFAULT (0) IS DAY 91].
!       ISOURC   SOURCE FLAG [0=SUN AND 1=MOON].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
      INTEGER LENN,IDAY,ISOURC,NLOS
      DOUBLE PRECISION H1ALT,H2ALT,OBSZEN,HRANGE,BETA,REARTH,           &
     &  HMIN,HMAX,CKRANG,BCKZEN,ANGLEM
      COMMON/CARD3/H1ALT(MLOS),H2ALT(MLOS),OBSZEN(MLOS),HRANGE(MLOS),   &
     &  BETA(MLOS),HMIN(MLOS),HMAX(MLOS),CKRANG(MLOS),BCKZEN(MLOS),     &
     &  REARTH,ANGLEM,LENN(MLOS),IDAY,ISOURC,NLOS

!     LOCAL VARIABLES:
!       RADIUS   EARTH CENTER DISTANCE [KM].
!       PRATLN   LOGARITHM OF THE PRESSURE RATIO.
      DOUBLE PRECISION RADIUS,PRATLN

!     INTEGRATE:
      RADIUS=REARTH+ZREF
      PRATLN=LOG(DBLE(PREF/P))
      HYSTAT=ZREF+PRATLN*RADIUS/(HSTAT/(RADIUS*DBLE(TREF+T))-PRATLN)
      RETURN
      END

      REAL FUNCTION FLG4PT(Z,Z1,Z2,Z3,Z4,F1,F2,F3,F4,LINIT)

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     FUNCTION FLG4PT DOES A 4 POINT LAGRANGE INTERPOLATION FOR
!     Z BETWEEN Z2 AND Z3, INCLUSIVE.  THE ABSCISSAE MUST BE
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

!     DECLARE INPUTS
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
      REAL Z,Z1,Z2,Z3,Z4,F1,F2,F3,F4
      LOGICAL LINIT

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     SAVED VARIABLES.
      REAL BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  &
     &  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, &
     &  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4
      SAVE BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  &
     &  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, &
     &  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4

!     LOCAL VARIABLES.
      REAL A,B,ZEXTRM,DSCRIM,ROOT,ZDIFF1,ZDIFF2,ZDIFF3,ZDIFF4,          &
     &  ZDIF12,ZDIF13,ZDIF14,ZDIF23,ZDIF24,ZDIF34,P1,P2,P3,P4,          &
     &  ZD12SQ,ZD13SQ,ZD14SQ,ZD23SQ,ZD24SQ,ZD34SQ,STORE

!     DATA:
!       SMALL    SMALL NUMBER TOLERANCE
      REAL SMALL
      DATA SMALL/1.E-30/

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
          COEF1=ZDIFF2*ZDIFF3*ZDIFF4
          COEF2=ZDIFF1*ZDIFF3*ZDIFF4
          COEF3=ZDIFF1*ZDIFF2*ZDIFF4
          COEF4=ZDIFF1*ZDIFF2*ZDIFF3
          ZDIF12=Z1-Z2
          ZDIF13=Z1-Z3
          ZDIF14=Z1-Z4
          ZDIF23=Z2-Z3
          ZDIF24=Z2-Z4
          ZDIF34=Z3-Z4
          DENOM1= ZDIF12*ZDIF13*ZDIF14
          DENOM2=-ZDIF12*ZDIF23*ZDIF24
          DENOM3= ZDIF13*ZDIF23*ZDIF34
          DENOM4=-ZDIF14*ZDIF24*ZDIF34

!         CONSTANTS USED IN EXTREMUM CHECK
          BCOEF1=Z2+Z3+Z4
          BCOEF2=Z3+Z4+Z1
          BCOEF3=Z4+Z1+Z2
          BCOEF4=Z1+Z2+Z3
          ZD12SQ=ZDIF12**2
          ZD13SQ=ZDIF13**2
          ZD14SQ=ZDIF14**2
          ZD23SQ=ZDIF23**2
          ZD24SQ=ZDIF24**2
          ZD34SQ=ZDIF34**2
          COEF11=(ZD23SQ+ZD24SQ+ZD34SQ)/2
          COEF22=(ZD13SQ+ZD14SQ+ZD34SQ)/2
          COEF33=(ZD12SQ+ZD14SQ+ZD24SQ)/2
          COEF44=(ZD12SQ+ZD13SQ+ZD23SQ)/2
          STORE=ZDIF13*ZDIF24+ZDIF14*ZDIF23
          COEF12=2*ZD34SQ+STORE
          COEF34=2*ZD12SQ+STORE
          STORE=ZDIF12*ZDIF34-ZDIF14*ZDIF23
          COEF13=2*ZD24SQ+STORE
          COEF24=2*ZD13SQ+STORE
          STORE=ZDIF12*ZDIF34+ZDIF13*ZDIF24
          COEF14=2*ZD23SQ-STORE
          COEF23=2*ZD14SQ-STORE
          CCOEF1=Z2*Z3+Z2*Z4+Z3*Z4
          CCOEF2=Z1*Z3+Z1*Z4+Z3*Z4
          CCOEF3=Z1*Z2+Z1*Z4+Z2*Z4
          CCOEF4=Z1*Z2+Z1*Z3+Z2*Z3

!         LINEAR INTERPOLATION FACTOR
          FACTOR=(Z-Z2)/(Z3-Z2)
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
              ZEXTRM=(CCOEF1*P1+CCOEF2*P2+CCOEF3*P3+CCOEF4*P4)/(2*B)
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  FLG4PT=F2+FACTOR*(F3-F2)
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
              ZEXTRM=(B+ROOT)/A
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  FLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
              ZEXTRM=(B-ROOT)/A
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  FLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
          ENDIF
      ENDIF

!     NO EXTREMA BETWEEN Z2 AND Z3.  USE 4-POINT LAGRANGE INTERPOLATION.
      FLG4PT=COEF1*P1+COEF2*P2+COEF3*P3+COEF4*P4
      RETURN

!     INCORRECT INPUT ERROR.
   10 CONTINUE
      WRITE(IPR,'(/A,/(18X,A,F12.4))')                                  &
     &  ' Error in FLG4PT:  ABSCISSAE OUT OF ORDER.',                   &
     &  ' Z1 =',Z1,' Z2 =',Z2,' Z  =',Z,' Z3 =',Z3,' Z4 =',Z4
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP ' Error in FLG4PT:  ABSCISSAE OUT OF ORDER.'
      END
