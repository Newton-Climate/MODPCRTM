      SUBROUTINE AERNSM(JPRT,GNDALT,MARIC1,MARK,ICH,LMODEL)

!     ROUTINE AERNSM DEFINES ALTITUDE, PRESSURE, TEMPERATURE,
!     MOLECULAR, AEROSOL, CLOUD, AND RAIN PROFILES.
      IMPLICIT NONE

!     ARGUMENTS:
      INTEGER JPRT,MARIC1,MARK,ICH(4)
      LOGICAL LMODEL
      DOUBLE PRECISION GNDALT

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'
      INCLUDE 'YPROP.h'

!     /M_PTWO/
!       PPROF    PRESSURE PROFILE [MB].
!       TPROF    TEMPERATURE PROFILE [K].
!       WH2O     H2O VOLUME MIXING RATIO PROFILE [PPMV].
!       WO3      O3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL PPROF,TPROF,WH2O,WO3
      COMMON/M_PTWO/PPROF(LAYDIM),TPROF(LAYDIM),WH2O(LAYDIM),WO3(LAYDIM)

!     /MDATXY/
!       WMOLXT   CROSS-SECTION MOLECULE DENSITY PROFILE [PPMV].
!       WMOLYT   AUXILIARY (Y) MOLECULE DENSITY PROFILE [PPMV].
      REAL WMOLXT,WMOLYT
      COMMON/MDATXY/WMOLXT(NMOLX,LAYDIM),WMOLYT(MMOLY,LAYDIM)

!     /M_UMIX/
!       WCO2     CO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WN2O     N2O VOLUME MIXING RATIO PROFILE [PPMV].
!       WCO      CO VOLUME MIXING RATIO PROFILE [PPMV].
!       WCH4     CH4 VOLUME MIXING RATIO PROFILE [PPMV].
!       WO2      O2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WNO      NO VOLUME MIXING RATIO PROFILE [PPMV].
!       WSO2     SO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WNO2     NO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WHNO3    HNO3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL WCO2,WN2O,WCO,WCH4,WO2,WNO,WSO2,WNO2,WNH3,WHNO3
      COMMON/M_UMIX/WCO2(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),              &
     &  WCH4(LAYDIM),WO2(LAYDIM),WNO(LAYDIM),WSO2(LAYDIM),              &
     &  WNO2(LAYDIM),WNH3(LAYDIM),WHNO3(LAYDIM)

!     /CARD1/
!       MODEL    MODEL ATMOSPHERE INDEX.
!       ITYPE    SLANT PATH TYPE.
!       IEMSCT   RADIATIVE TRANSFER MODE.
!                  0 FOR TRANSMITTANCE
!                  1 FOR THERMAL EMISSION ONLY
!                  2 FOR THERMAL EMISSION PLUS SOLAR SCATTER
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!                  4 FOR SOLAR SCATTER ONLY
!       M1       MODEL ATMOSPHERE FOR PRESSURE & TEMPERATURE PROFILES.
!       M2       MODEL ATMOSPHERE FOR H2O PROFILE.
!       M3       MODEL ATMOSPHERE FOR O3 PROFILE.
!       I_RD2C   READ CARD 2C, 2C1, ... IF EQUAL 1; SKIP IF EQUAL TO 0.
!       NOPRNT   PRINT FLAG.
!       MODTRN   MODTRAN BAND MODEL FLAG.
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT,MODTRN

!     /CARD1A/
      INTEGER M4,M5,M6,MDEF,IRD1,IRD2
      COMMON/CARD1A/M4,M5,M6,MDEF,IRD1,IRD2

!     /CARD1B/
      INTEGER JUNITP,JUNITT,JUNIT,JLOW
      REAL WMOL
      COMMON/CARD1B/JUNITP,JUNITT,JUNIT(13),WMOL(12),JLOW

!     /JM1B/
      CHARACTER JCHAR*17
      COMMON/JM1B/JCHAR

!     /CARD2/
!       IHAZE    BOUNDARY LAYER AEROSOL MODEL NUMBER.
!       ISEASN   SEASON NUMBER (1=SPRING-SUMMER, 2=FALL-WINTER).
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
!       ICSTL    COASTAL AIRMASS MODEL NUMBER.
!       ICLD     CLOUD MODEL NUMBER.
!       IVSA     VERTICAL STRUCTURE ALGORITHM (0=OFF, 1=ON).
!       VIS      SURFACE VISIBILITY (GROUND METEOROLOGICAL RANGE) [KM].
!       WSS      CURRENT WIND SPEED (M/S).
!       WHH      24-HOUR WIND SPEED (M/S).
!       RAINRT   RAIN RATE (MM/HR)
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      LOGICAL LSAP
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT,LSAP

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      DOUBLE PRECISION CTHIK,CALT
      REAL CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,        &
     &  CCOLIP,CHUMID,ASYMWD,ASYMIP
      INTEGER IREG,IREGC
      DOUBLE PRECISION ALTB
      COMMON/CARD2D/IREG(4),ALTB(4),IREGC(4)

!     /CNTRL/
!       NSEG     NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       THERML   FLAG TO CALCULATE THERMAL SCATTER.
      INTEGER NSEG,ML,MLFLX,IMULT
      LOGICAL THERML
      COMMON/CNTRL/NSEG(0:MLOSP1),ML,MLFLX,IMULT,THERML

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

!     /DEN/
!       DENSTY   PROFILE LEVEL DENSITIES [ATM CM / KM FOR MOST SPECIES].
      REAL DENSTY
      COMMON/DEN/DENSTY(0:MEXTXY,1:LAYDIM)

!     /ZVSALY/
      DOUBLE PRECISION ZVSA
      REAL RHVSA,AHVSA
      INTEGER IHVSA
      COMMON/ZVSALY/ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)

!     /MDLZ/
      REAL HMDLZ
      COMMON/MDLZ/HMDLZ(8)

!     /TITL/
      CHARACTER HHAZE(16)*20,HSEASN(2)*20,HVULCN(8)*20,                 &
     &  HMET(2)*20,HMODEL(0:8)*20
      COMMON/TITL/HHAZE,HSEASN,HVULCN,HMET,HMODEL

!     /CD1BXY/
!       JUNITX   INPUT UNIT FLAG FOR CROSS-SECTION MOLECULES.
!       WMOLX    X-SECTION MOLECULAR DENSITY [PPMV AFTER JUBRAN CALL].
!       JUNITY   INPUT UNIT FLAG FOR AUXILIARY (Y) MOLECULES.
!       WMOLY    AUXILIARY MOLECULAR DENSITY [PPMV AFTER JUBRAN CALL].
      INTEGER JUNITX,JUNITY
      REAL WMOLX,WMOLY
      COMMON/CD1BXY/JUNITX,WMOLX(NMOLX),JUNITY,WMOLY(MMOLY)

!     /JM2C3/
      REAL AHAZE(4),EQLWCZ,RRATZ
      INTEGER IHA1,ICLD1,IVUL1,ISEA1,ICHR
      COMMON/JM2C3/AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,                       &
     &   IVUL1,ISEA1,ICHR

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /TITL/
      EXTERNAL DEVCBD,TITLE

!     FUNCTIONS:
!       JOU      TRANSLATE UNIT SPECIFIER CHARACTER INTO INTEGER LABEL.
!       LRDSAP   TRUE FOR ATMOSPHERIC LEVELS WITH SPECTRAL PROFILE DATA.
      INTEGER JOU
      REAL CLDPRF,AERPRF
      LOGICAL LRDSAP

!     LOCAL VARIABLES:
!       JDEF     VALUE OF JUNIT FOR MDEF SPECIES.
!       LCIRZ    FLAG, TRUE IF CIRRUS ALTITUDES NEED TO BE DEFINED.
!       WMOLYS   TEMPORARY PLACE HOLDER FOR WMOLY ARRAY VALUES.
!       RD_2C1   READ CARD 2C1 LOGICAL SWITCH.
      CHARACTER HHOL*20,AHOL1*20,AHOL2*20
      INTEGER JDEF,ICLDS,ITYAER,IC1,NUMAER,K,I
      LOGICAL LCIRZ,LDESRT,RD_2C1
      REAL WH100,RH,VIS1,HAZ1,HAZ2
      DOUBLE PRECISION FAC,CLDD,CLD0,CLD1,CLD2,CLD3

!     LOCAL ARRAYS:
      INTEGER ITY1(LAYDIM),IH1(LAYDIM),IS1(LAYDIM),IVL1(LAYDIM)
      REAL AHAST(LAYDIM),WMOLYS(MMOLY)
      DOUBLE PRECISION ZGN(LAYDIM)

!     DATA:
      CHARACTER AHAHOL(13)*20
      DOUBLE PRECISION CLDTOP(10)
      DATA AHAHOL/             'CUMULUS             ',                  &
     &  'ALTOSTRATUS         ','STRATUS             ',                  &
     &  'STRATUS STRATO CUM  ','NIMBOSTRATUS        ',                  &
     &  'DRIZZLE 2.0 MM/HR   ','LT RAIN 5.0 MM/HR   ',                  &
     &  'MOD RAIN 12.5 MM/HR ','HEAVY RAIN 25 MM/HR ',                  &
     &  'EXTREME RAIN 75MM/HR','USER ATMOSPHERE     ',                  &
     &  'USER RAIN NO CLOUD  ','CIRRUS CLOUD        '/
      DATA CLDTOP/3.D0,3.D0,1.D0,2.D0,.66D0,1.D0,.66D0,.66D0,3.D0,3.D0/

!     INITIALIZATIONS:
      IC1=1
      NUMAER=7
      IF(LMODEL)THEN
          RD_2C1=.FALSE.
      ELSEIF(I_RD2C.NE.1)THEN
          RETURN
      ELSEIF(IVSA.EQ.1)THEN

!         VERTICAL STRUCTURE ALGORITHM:
          IF(MODEL.EQ.0)THEN
              WRITE(IPR,'(/2A)')' ERROR: ',                             &
     &          ' MODEL equals 0 and army (VSA) model cannot mix.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' MODEL equals 0 and army (VSA) model cannot mix.'
          ENDIF
          RD_2C1=.FALSE.
          IRD1=0
          IRD2=0
          ML=ML+10-JLOW
          IF(ML.GT.LAYDIM)THEN
              WRITE(IPR,'(/2A)')' WARNING:  ML exceeds parameter',      &
     &          ' LAYDIM and army (VSA) model top layer was truncated.'
              ML=LAYDIM
          ENDIF
          ZVSA(10)=ZVSA(9)+.01D0
          RHVSA(10)=0.
          AHVSA(10)=0.
          IHVSA(10)=0
      ELSE
          RD_2C1=.TRUE.
      ENDIF

!     INITIALIZE COMMON /CARD2D/:
      IREGC(1)=0
      IREGC(2)=0
      IREGC(3)=0
      IREGC(4)=0
      ALTB(1)=0.D0
      ALTB(2)=0.D0
      ALTB(3)=0.D0
      ALTB(4)=0.D0
      LDESRT=.TRUE.
      IF(ICLD.EQ.18 .OR. ICLD.EQ.19)CALL CIRR18(ICLD,LCIRZ)

!     DEFAULTS:
      IF(IVULCN.LE.0)IVULCN=1
      IF(ISEASN.LE.0)ISEASN=1
      IF(.NOT.LJMASS)WRITE(IPR,'(/A,I3)')' MODEL ATMOSPHERE NO.',MODEL
      IF(LMODEL)THEN
          CALL FLAYZ(ML,ICLD,GNDALT,IVSA)
      ELSE
          IF(MDEF.EQ.-1)THEN
              JDEF=10
          ELSE
              JDEF=6
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,'(/10X,A)')                          &
     &      ' MODEL 0, 7, or 8 USER INPUT DATA:'
      ENDIF

!     LOOP OVER LAYERS:
      IF(ML.GT.LAYDIM)THEN
          WRITE(IPR,'(/2A,I5,A,/19X,A,I5,A)')' Error in AERNSM: ',      &
     &      ' Number of atmospheric levels (ML =',ML,') exceeds',       &
     &      'parameter LAYDIM (=',LAYDIM,').  LAYDIM must be increased.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' Error:  Parameter LAYDIM must be increased.'
      ENDIF
      DO K=1,ML
          LRHSET(K)=.FALSE.
          RH=0.
          WH2O(K)=0.
          WO3(K)=0.
          IHA1=0
          ICLD1=0
          ISEA1=0
          IVUL1=0
          AHAZE(1)=0.
          AHAZE(2)=0.
          AHAZE(3)=0.
          AHAZE(4)=0.
          EQLWCZ=0.
          RRATZ=0.
          ICHR=0
          DENSTY(16,K)=0.
          WCO2(K)=0.
          WCO(K)=0.
          WCH4(K)=0.
          WN2O(K)=0.
          WO2(K) =0.
          WNH3(K)=0.
          WNO (K)=0.
          WNO2(K)=0.
          WSO2(K)=0.
          WHNO3(K)= 0.
          WMOL(1)=0.
          WMOL(2)=0.
          WMOL(3)=0.
          WMOL(4)=0.
          WMOL(5)=0.
          WMOL(6)=0.
          WMOL(7)=0.
          WMOL(8)=0.
          WMOL(9)=0.
          WMOL(10)=0.
          WMOL(11)=0.
          WMOL(12)=0.
          JCHAR='                '
          IF(RD_2C1)THEN

!             READ IN USER-SPECIFIED ATMOSPHERE.
!             FOR MOLECULAR SPECIES, JCHAR IS DEFINED AS FOLLOWS:
!               JCHAR   JUNIT
!               -----   -----
!               " ",A     10    VOLUME MIXING RATIO (PPMV)
!                 B       11    NUMBER DENSITY (CM-3)
!                 C       12    MASS MIXING RATIO (GM(K)/KG(AIR))
!                 D       13    MASS DENSITY (GM M-3)
!                 E       14    PARTIAL PRESSURE (MB)
!                 F       15    DEW POINT TEMPERATURE (K) - H2O ONLY
!                 G       16    DEW POINT TEMPERATURE (C) - H2O ONLY
!                 H       17    RELATIVE HUMIDITY (%) - H2O ONLY
!                1-6     1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE

!             OTHER 'JCHAR' SPECIFICATIONS -
!               JCHAR   JUNIT
!               -----   -----
!               " ",A     10    PRESSURE IN (MB)
!                 B       11    PRESSURE IN (ATM)
!                 C       12    PRESSURE IN (TORR)
!                1-6     1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE

!               " ",A     10    AMBIENT TEMPERATURE (K)
!                 B       11    AMBIENT TEMPERATURE (C)
!                 C       12    DELTA TEMPERATURE (K)
!                1-6     1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
              IF(LJMASS)THEN
                  CALL INITCD('CARD2C1')
              ELSE
                  READ(IRD,'(6F10.0,A17)')ZM(K),                        &
     &              PPROF(K),TPROF(K),WMOL(1),WMOL(2),WMOL(3),JCHAR
                  WRITE(IPR,'(F10.5,1P,5E10.3,10X,A17)')ZM(K),          &
     &              PPROF(K),TPROF(K),WMOL(1),WMOL(2),WMOL(3),JCHAR
              ENDIF
              IF(IRD1.EQ.1)THEN
                  IF(LJMASS)THEN
                      CALL INITCD('CARD2C2')
                  ELSE
                      READ(IRD,'((8F10.0))')(WMOL(I),I=4,12)
                      WRITE(IPR,'((1P,8E10.3))')(WMOL(I),I=4,12)
                  ENDIF
                  IF(MDEF.EQ.2)THEN

!                     THE EXTRA SPECIES (I.E. THE WMOLX SPECIES) ARE
!                     READ IF MDEF=2 AND IRD1=1, 8 SPECIES PER LINE:
                      IF(LJMASS)THEN
                          CALL INITCD('CARD2C2X')
                      ELSE
                          READ(IRD,'((8F10.0))')(WMOLX(I),I=1,NMOLX)
                          WRITE(IPR,'(1P,(8E10.3))')(WMOLX(I),I=1,NMOLX)
                      ENDIF
                  ELSE
                      DO I=1,NMOLX
                          WMOLX(I)=0.
                      ENDDO
                  ENDIF
                  IF(NMOLYC.GT.0)THEN

!                     THE NSPECY Y-SPECIES ARE READ IN, 8 PER LINE.
                      READ(IRD,'((8F10.0))')(WMOLYS(I),I=1,NMOLYC)
                      WRITE(IPR,'((1P,8E10.3))')(WMOLYS(I),I=1,NMOLYC)
                      DO I=1,NMOLY
                          WMOLY(I)=WMOLYS(MAPY_F(I))
                      ENDDO
                  ENDIF
              ELSE
                  DO I=1,NMOLX
                      WMOLX(I)=0.
                  ENDDO
              ENDIF
              IF(IRD2.EQ.1)THEN

!                 CARD2C3:  ONLY ONE OF IHA1, ICLD1  OR IVUL1 IS
!                           ALLOWED. IF IHA1>0, OTHERS IGNORED; IF
!                           IHA1=0 AND ICLD1=0, USE ICLD1.  IF AHAZE
!                           AND EQLWCZ ARE BOTH ZERO, DEFAULT PROFILES
!                           ARE LOADED FROM IHA1, ICLD1, OR IVUL1.

!                   AHAZE    AEROSOL EXTINCTION AT 550 NM [KM-1].
!                   EQLWCZ   LIQUID WATER CONTENT AT ALTITUDE Z FOR
!                            AEROSOL, CLOUD OR FOG MODELS [PPMV].
!                   RRATZ    RAIN RATE AT ALTITUDE Z [MM/HR].
!                   IHA1     BOUNDARY LAYER AEROSOL MODEL USED FOR
!                            SPECTRAL EXTINCTION.
!                   IVUL1    STRATOSPHERIC AEROSOL MODEL USED FOR
!                            SPECTRAL EXTINCTION.
!                   ICLD1    CLOUD MODEL USED FOR SPECTRAL EXTINCTION.
!                   ISEA1    AEROSOL SEASON CONTROL FOR ALTITUDE Z.
!                   ICHR     AEROSOL PROFILE REGION SWITCH FOR IHA=7.
                  IF(LJMASS)THEN
                      CALL INITCD('CARD2C3')
                  ELSE
                      READ(IRD,'(10X,3F10.0,5I5)')AHAZE(1),             &
     &                  EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1,ICHR
                      WRITE(IPR,'(10X,3F10.3,5I5)')AHAZE(1),            &
     &                  EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1,ICHR
                  ENDIF
                  AHAZE(2)=0.
                  AHAZE(3)=0.
                  AHAZE(4)=0.
              ELSEIF(IRD2.EQ.2)THEN

!                 READ 4 AEROSOL PROFILES:
                  IF(LJMASS)THEN
                      CALL INITCD('CARD2C3')
                  ELSE
                      READ(IRD,'(10X,F10.0,10X,4F10.0)')                &
     &                  AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
                      WRITE(IPR,'(10X,F10.3,10X,4F10.3)')               &
     &                  AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
                  ENDIF
              ENDIF
          ENDIF
          IF((IHA1.EQ.0 .AND. ICLD1.NE.11) .OR. IHA1.NE.7)ICHR=0
          IF(MODEL.EQ.0)THEN

!             HORIZONTAL (CONSTANT PRESSURE) PATH.
              HMDLZ(1)=SNGL(ZM(K))
              HMDLZ(2)=PPROF(K)
              HMDLZ(3)=TPROF(K)
              HMDLZ(4)=WMOL(1)
              HMDLZ(5)=WMOL(2)
              HMDLZ(6)=WMOL(3)
              HMDLZ(7)=AHAZE(1)
          ENDIF

!         JUNITP, JUNITT, JUNIT() AND JUNITX:
          IF(RD_2C1)THEN
              JUNITP=JOU(JCHAR(1:1),M1)
              JUNITT=JOU(JCHAR(2:2),M1)
              JUNIT(1)=JOU(JCHAR(3:3),M2)
              JUNIT(2)=JOU(JCHAR(4:4),JDEF)
              JUNIT(3)=JOU(JCHAR(5:5),M3)
              JUNIT(4)=JOU(JCHAR(6:6),M5)
              JUNIT(5)=JOU(JCHAR(7:7),M6)
              JUNIT(6)=JOU(JCHAR(8:8),M4)
              JUNIT(7)=JOU(JCHAR(9:9),JDEF)
              JUNIT(8)=JOU(JCHAR(10:10),JDEF)
              JUNIT(9)=JOU(JCHAR(11:11),JDEF)
              JUNIT(10)=JOU(JCHAR(12:12),JDEF)
              JUNIT(11)=JOU(JCHAR(13:13),JDEF)
              JUNIT(12)=JOU(JCHAR(14:14),JDEF)
              JUNIT(13)=JOU(JCHAR(15:15),JDEF)
              JUNITX=JOU(JCHAR(16:16),JDEF)
              JUNITY=JOU(JCHAR(17:17),JDEF)
          ELSE
              JUNITP=M1
              JUNITT=M1
              JUNIT(1)=M2
              JUNIT(2)=6
              JUNIT(3)=M3
              JUNIT(4)=M5
              JUNIT(5)=M6
              JUNIT(6)=M4
              JUNIT(7)=6
              JUNIT(8)=6
              JUNIT(9)=6
              JUNIT(10)=6
              JUNIT(11)=6
              JUNIT(12)=6
              JUNIT(13)=6
              JUNITX=6
              JUNITY=6
          ENDIF
          IF(IVSA.EQ.1 .AND. .NOT.LMODEL)THEN
              IF(MODEL.EQ.8)THEN
                  WRITE(IPR,'(/2A,/(18X,A))')' Error in AERNSM:  The',  &
     &              ' pressure-dependent radiosonde data (MODEL=8)',    &
     &              ' option cannot be used with the Army Vertical',    &
     &              ' Structure Algorithm (IVSA=1).'
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP 'Error: MODEL=8 cannot be coupled with IVSA=1.'
              ENDIF
              CALL VSANSM(K,AHAZE(1),IHA1,ZM(K))
          ELSE
              CALL CHECKP(PPROF(K),JUNITP)
              CALL CHECKT(TPROF(K),JUNITT,M1)
              IF(MODEL.EQ.8)THEN
                  CALL DFLTP(K,GNDALT)
              ELSE
                  CALL DEFALT(ZM(K),PPROF(K),TPROF(K))
              ENDIF
              CALL JUBRAN(PPROF(K),TPROF(K))
              WH2O(K)=WMOL(1)
              WCO2(K)=WMOL(2)
              WO3(K)=WMOL(3)
              WN2O(K)=WMOL(4)
              WCO(K)=WMOL(5)
              WCH4(K)=WMOL(6)
              WO2(K)=WMOL(7)
              WNO(K)=WMOL(8)
              WSO2(K)=WMOL(9)
              WNO2(K)=WMOL(10)
              WNH3(K)=WMOL(11)
              WHNO3(K)=WMOL(12)
              DO I=1,NMOLX
                  WMOLXT(I,K)=WMOLX(I)
              ENDDO
              DO I=1,NMOLY
                  WMOLYT(I,K)=WMOLY(I)
              ENDDO
          ENDIF

!         WITH ALTITUDE SET, RAIN RATE CAN BE DEFINED:
          IF(IRD2.EQ.0 .AND. ZM(K).LE.6.D0)RRATZ=RAINRT

!         SANITY CHECKS ON PRESSURE AND TEMPERATURE:
          IF(PPROF(K).LE.0.)THEN
              WRITE(IPR,'(/A,I3,A,1P,E12.4,A)')' Error in routine AER'//&
     &          'NSM:  The pressure at level',K,' is',PPROF(K),' mbars.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' Non-positive pressure encountered.'
          ELSEIF(K.GT.1)THEN
              IF(PPROF(K).GE.PPROF(K-1))                                &
     &          WRITE(IPR,'(/2A,I3,/30X,2(A,1P,E12.4,A,I3))')           &
     &          ' Warning from routine AERNSM:  A pressure inversion',  &
     &          ' occurs between level',K-1,' (P=',PPROF(K-1),          &
     &          ' mbar) and level',K,' (P=',PPROF(K),' mbar).'
          ENDIF
          IF(TPROF(K).LE.0.)THEN
              WRITE(IPR,'(/A,I3,A,F9.3,A)')' Error in routine AER'//    &
     &          'NSM:  The temperature at level',K,' is',TPROF(K),' K.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' Non-positive temperature encountered.'
          ELSEIF(TPROF(K).LT.30. .OR. TPROF(K).GT.400.)THEN
              WRITE(IPR,'(/A,I3,A,F9.3,A)')' Warning from routine AER'//&
     &          'NSM:  The temperature at level',K,' is',TPROF(K),' K.'
          ENDIF

!         CHECK GROUND ALTITUDE:
          IF(GNDALT.NE.ZM(1))THEN

!             FIX GROUND ALTITUDE TO THE BOTTOM OF ATMOSPHERIC PROFILE.
              IF(ABS(GNDALT-ZM(1)).GT..0001 .AND. .NOT.LJMASS)          &
     &          WRITE(IPR,'(2A,F8.4,A,/10X,A,F8.4,A)')' WARNING:  THE', &
     &          ' INPUT GROUND ALTITUDE (',GNDALT,'KM) IS BEING SET',   &
     &          ' TO THE BOTTOM OF THE ATMOSPHERIC PROFILE,',ZM(1),'KM.'
              GNDALT=ZM(1)
          ENDIF

!         PLACE ORIGINAL ALTITUDES IN ZGN
          IF(ZM(K).LT.6.D0)THEN
              ZGN(K)=6*(ZM(K)-GNDALT)/(6-GNDALT)
          ELSE
              ZGN(K)=ZM(K)
          ENDIF
          ICLDS=ICLD1
          IF(ICLD1.EQ.0)ICLD1=ICLD
          IF(ICLD1.GT.11)ICLD1=0
          IF(IHA1.NE.0)THEN
              IVUL1=0
              ICLD1=0
          ELSEIF(ICLD1.NE.0)THEN
              IVUL1=0
          ENDIF
          IF(AHAZE(1).EQ.0. .AND. EQLWCZ.EQ.0.)THEN
              IF(IVSA.EQ.1 .AND. ICLD1.EQ.0)THEN
                  IF(MODEL.LT.7)CALL LAYVSA(K,RH,AHAZE(1),IHA1)
              ELSE
                  CALL LAYCLD(K,EQLWCZ,RRATZ,ICLD1,GNDALT)
                  IF(RAINRT.GT.0. .AND. ZM(K).LT.6.D0)RRATZ=RAINRT
                  IF(ICLD1.GE.1 .AND. ICLD1.LE.10)THEN
                      IF(ZM(K).GT.CLDTOP(ICLD1)+GNDALT)RRATZ=0.
                  ENDIF
              ENDIF
          ENDIF

!         DENSTY(16,K):
          IF(ICLDS.EQ.18 .OR. ICLDS.EQ.19)THEN
              IF(AHAZE(1).GT.0.)THEN
                  DENSTY(16,K)=AHAZE(1)
                  AHAZE(1)=0.
              ELSEIF(EQLWCZ.GT.0)THEN
                  IF(ICLDS.EQ.18)THEN
                      DENSTY(16,K)=EQLWCZ/5.811E-2
                  ELSE
                      DENSTY(16,K)=EQLWCZ/3.446E-3
                  ENDIF
                  EQLWCZ=0.
              ENDIF
          ELSEIF(ICLDS.EQ.0 .AND. (ICLD.EQ.18 .OR. ICLD.EQ.19))THEN
              IF(LCIRZ)THEN
                  CLDD=CTHIK/10
                  CLD0=CALT-CLDD/2
                  IF(CLD0.LE.GNDALT)CLD0=GNDALT
                  CLD1=CLD0+CLDD
                  CLD2=CLD0+CTHIK
                  CLD3=CLD1+CTHIK
                  LCIRZ=.FALSE.
              ENDIF
              IF(ZM(K).LE.CLD0 .OR. ZM(K).GE.CLD3)THEN
                  DENSTY(16,K)=0.
              ELSEIF(ZM(K).LT.CLD1)THEN
                  DENSTY(16,K)=CEXT*SNGL((ZM(K)-CLD0)/CLDD)
              ELSEIF(ZM(K).LE.CLD2)THEN
                  DENSTY(16,K)=CEXT
              ELSE
                  DENSTY(16,K)=CEXT*SNGL((CLD3-ZM(K))/CLDD)
              ENDIF
          ENDIF
          DENSTY(66,K)=EQLWCZ
          DENSTY(67,K)=0.
          IF(ICLDS.EQ.0 .AND. EQLWCZ.EQ.0.)ICLD1=0
          DENSTY(3,K)=RRATZ
          IF(LMODEL .AND. (EQLWCZ.GT.0. .OR. RRATZ.GT.0.))RH=100.
          AHAST(K)=AHAZE(1)

!           IHA1    IHAZE FOR THIS LAYER
!           ISEA1   ISEASN FOR THIS LAYER
!           IVUL1   IVULCN FOR THE LAYER
          IF(ISEA1.EQ.0)ISEA1=ISEASN
          IF(IHA1.GT.0)THEN
              ITYAER=IHA1
          ELSE
              ITYAER=IHAZE
          ENDIF
          IF(IVUL1.GT.0)THEN
              IVULCN=IVUL1
          ELSE
              IVUL1=IVULCN
          ENDIF
          IF(K.GT.1)THEN
              IF(ICHR.NE.1)THEN
                  IF(ICLD1.EQ.IREGC(IC1))THEN
                      IF(IHA1.EQ.0 .AND. ICLD1.EQ.0)THEN
                          IF(ZGN(K).GT.30.00001D0)THEN
                              ITYAER=19
                          ELSEIF(ZGN(K).GT.10.00001D0)THEN
                              ITYAER=IVULCN+10
                          ELSEIF(ZGN(K).GT.2.00001D0)THEN
                              ITYAER=6
                          ENDIF
                          IF(ITYAER.EQ.ICH(IC1))GOTO 10
                      ELSE
                          NUMAER=7
                          IF(IC1.GT.1)NUMAER=IC1+10
                          IF(IHA1.EQ.0 .OR. IHA1.EQ.ICH(IC1))GOTO 10
                      ENDIF
                  ELSEIF(ICLD1.NE.0)THEN
                      IF(ICLD1.EQ.IREGC(1))THEN
                          NUMAER=7
                          ALTB(1)=ZM(K)
                          GOTO 20
                      ELSEIF(IC1.GT.1 .AND. ICLD1.EQ.IREGC(2))THEN
                          NUMAER=12
                          ALTB(2)=ZM(K)
                          GOTO 20
                      ELSEIF(IC1.GT.2 .AND. ICLD1.EQ.IREGC(3))THEN
                          NUMAER=13
                          ALTB(3)=ZM(K)
                          GOTO 20
                      ENDIF
                  ELSE
                      IF(IHA1.EQ.0 .AND. ICLD1.EQ.0)THEN
                          IF(ZGN(K).GT.30.00001D0)THEN
                              ITYAER=19
                          ELSEIF(ZGN(K).GT.10.00001D0)THEN
                              ITYAER=IVULCN+10
                          ELSEIF(ZGN(K).GT.2.00001D0)THEN
                              ITYAER=6
                          ENDIF
                      ENDIF
                      IF(ITYAER.EQ.ICH(1))THEN
                          NUMAER=7
                          ALTB(1)=ZM(K)
                          GOTO 20
                      ELSEIF(IC1.GT.1 .AND. ITYAER.EQ.ICH(2))THEN
                          NUMAER=12
                          ALTB(2)=ZM(K)
                          GOTO 20
                      ELSEIF(IC1.GT.2 .AND. ITYAER.EQ.ICH(3))THEN
                          NUMAER=13
                          ALTB(3)=ZM(K)
                          GOTO 20
                      ENDIF
                  ENDIF
              ENDIF
              IF(IC1.LT.4)THEN
                  IC1=IC1+1
                  NUMAER=IC1+10
              ELSE
                  IC1=4
                  NUMAER=14
                  ITYAER=ICH(IC1)
              ENDIF
          ENDIF
   10     CONTINUE
          ICH(IC1)=ITYAER
          IREGC(IC1)=ICLD1
          ALTB(IC1)=ZM(K)
   20     CONTINUE

!         SATURATED WATER VAPOR DENSITY [GM / M3]
          WH100=TZERO/TPROF(K)
          WH100=WH100*EXP(18.9766-(14.9595+2.43882*WH100)*WH100)
          IF(RH.GT.0.)WH2O(K)=.01*RH*WH100
          DENSTY(7,K)=0.
          DENSTY(12,K)=0.
          DENSTY(13,K)=0.
          DENSTY(14,K)=0.
          DENSTY(15,K)=0.
          RELHUM(K)=0.
          IF(WH2O(K).GT.0.)THEN
              RELHUM(K)=100.*WH2O(K)/WH100
              IF(RELHUM(K).GT.100.)THEN
                  IF(RELHUM(K).GT.100.05)WRITE(IPR,'(/A,2(F12.5,A))')   &
     &              ' WARNING:  The relative humidity at',ZM(K),        &
     &              'km was reset from',RELHUM(K),'% to 100%'
                  RELHUM(K)=100.
                  WH2O(K)=WH100
              ELSEIF(RELHUM(K).LT.0.)THEN
                  WRITE(IPR,'(/A,2(F12.5,A))')                          &
     &              ' WARNING:  The relative humidity at',ZM(K),        &
     &              'km was reset from',RELHUM(K),'% to 0%'
                  RELHUM(K)=0.
                  WH2O(K)=0.
              ENDIF
          ENDIF
          VIS1=VIS
          IF(IRD2.EQ.2)THEN
              DENSTY(7,K)=AHAZE(1)
              DENSTY(12,K)=AHAZE(2)
              DENSTY(13,K)=AHAZE(3)
              DENSTY(14,K)=AHAZE(4)
          ELSEIF(AHAZE(1).GT.0.)THEN
              DENSTY(NUMAER,K)=AHAZE(1)
              IF(ITYAER.NE.3 .AND. ITYAER.NE.10)GOTO 30
          ENDIF
          IF(ITYAER.EQ.3 .AND. MARIC1.EQ.0)THEN
              CALL MARINE(VIS1,MODEL,RELHUM(K),                         &
     &          WSS,WHH,ICSTL,EXTC,ABSC,IC1)
              IREG(IC1)=1
              VIS=VIS1
              MARIC1=IC1
              MARK=K
          ELSEIF(ITYAER.EQ.10 .AND. LDESRT)THEN
              CALL DESATT(WSS,VIS1)
              IREG(IC1)=1
              VIS=VIS1
              LDESRT=.FALSE.
          ENDIF
          IF(IRD2.NE.2 .AND. AHAZE(1).LE.0.)THEN
              IF(IHA1.LE.0)IHA1=IHAZE
              IF(EQLWCZ.GT.0.)THEN
                  DENSTY(NUMAER,K)                                      &
     &              =CLDPRF(EQLWCZ,RELHUM(K),ICLD1,IHA1,IC1,ICH(IC1))
              ELSE
                  IF(ZM(K).LT.6.D0)THEN
                      I=INT(ZGN(K)+1.D-6)+1
                      FAC=ZGN(K)-(I-1)
                  ELSEIF(ZM(K).GE.70.D0)THEN
                      I=32
                      FAC=(ZM(K)-70)/30
                  ELSEIF(ZM(K).GE.50.D0)THEN
                      I=31
                      FAC=(ZM(K)-50)/20
                  ELSEIF(ZM(K).GE.25.D0)THEN
                      FAC=(ZM(K)-25)/5
                      I=INT(FAC)
                      FAC=FAC-I
                      I=I+26
                  ELSE
                      I=INT(ZM(K))
                      FAC=ZM(K)-I
                      I=I+1
                      IF(I.LE.0)THEN

!                         THIS FIX IS FOR THE ALTITUDES -1.0 KM.  IN
!                         THIS CASE, I.LE.0 AND FAC SHOULD BE ZM(K)-0.
!                         SET I T0 1 AND FAC TO THE APPROPRIATE VALUE.
                          FAC=ZM(K)
                          I=1
                      ENDIF
                  ENDIF
                  HAZ1=AERPRF(I,VIS1,IHA1,ISEA1,IVUL1)
                  DENSTY(NUMAER,K)=0.
                  IF(HAZ1.GT.0.)THEN
                      I=I+1
                      HAZ2=AERPRF(I,VIS1,IHA1,ISEA1,IVUL1)
                      IF(HAZ2.GT.0.)                                    &
     &                  DENSTY(NUMAER,K)=HAZ1*(HAZ2/HAZ1)**SNGL(FAC)
                  ENDIF
              ENDIF
          ENDIF
   30     CONTINUE
          ITY1(K)=ITYAER
          IF(AHAZE(1).EQ.0.)THEN
              IH1(K)=IHA1
          ELSE
              IH1(K)=-99
          ENDIF
          IS1(K)=ISEA1
          IVL1(K)=IVUL1

!         CHECK FOR SPECTRAL AEROSOL PROFILES (SAP):
          IF(LSAP)THEN
              IF(LRDSAP(K))THEN

!                 AEROSOL DATA OVERWRITTEN BY SPECTRAL AEROSOL PROFILE:
                  DENSTY(7,K)=0.
                  DENSTY(12,K)=0.
                  DENSTY(13,K)=0.
                  DENSTY(14,K)=0.
              ENDIF
          ENDIF
      ENDDO
      JPRT=1

!     CHECK GROUND ALTITUDE:
      IF(GNDALT.LT.ZM(1))THEN

!         RAISE GROUND ALTITUDE TO THE BOTTOM OF ATMOSPHERE:
          IF(GNDALT.LT.ZM(1)-.00005D0)                                  &
     &      WRITE(IPR,'(2A,F8.4,A,/10X,A,F8.4,A)')' WARNING:  The',     &
     &      ' input ground altitude (',GNDALT,' km) is being raised',   &
     &      ' to the bottom of the atmosphere,',ZM(1),' km.'
          GNDALT=ZM(1)
      ENDIF
      IF(LJMASS .OR. (LMODEL .AND. IVSA.EQ.0 .AND. ICLD.EQ.0            &
     &  .AND. RAINRT.EQ.0. .AND. GNDALT.EQ.0.D0))RETURN
      JPRT=0
      IF(IVSA.EQ.1)THEN
          HHOL='VSA DEFINED         '
      ELSEIF(ICLD.GE.18)THEN
          HHOL=AHAHOL(13)
      ELSEIF(ICLD.LE.0 .OR. ICLD.GT.12)THEN
          HHOL=AHAHOL(12)
      ELSE
          HHOL=AHAHOL(ICLD)
      ENDIF
      IF(ICLD.NE.0)THEN
          WRITE(IPR,'(/2A)')' CLOUD AND/OR RAIN TYPE CHOSEN IS ',HHOL
          WRITE(IPR,'(//(2A))')                                         &
     &      '      Z         P        T     REL H    H2O   ',           &
     &      '  CLD AMT   RAIN RATE                      AEROSOL',       &
     &      '     (KM)      (MB)     (K)     (%)  (GM / M3)',           &
     &      ' (GM / M3)  (MM / HR) TYPE                 PROFILE',       &
     &      '                              [Before scaling]'
      ELSE
          WRITE(IPR,'(//(2A))')                                         &
     &      '      Z         P        T     REL H    H2O   ',           &
     &      '                       AEROSOL',                           &
     &      '     (KM)      (MB)     (K)     (%)  (GM / M3)',           &
     &      '  TYPE                 PROFILE',                           &
     &      '                              [Before scaling]'
      ENDIF
      DO K=1,ML
          IF(ITY1(K).LE.0)THEN
              ITYAER=1
          ELSEIF(ITY1(K).EQ.18)THEN
              ITYAER=13
          ELSEIF(ITY1(K).GE.16 .AND. ITY1(K).LE.19)THEN
              ITYAER=11
          ELSE
              ITYAER=ITY1(K)
          ENDIF
          IHA1=IH1(K)
          ISEA1=IS1(K)
          IVUL1=IVL1(K)
          IF((IVSA.EQ.1 .AND. K.LE.9) .OR. DENSTY(66,K).GT.0. .OR.      &
     &      DENSTY(3,K).GT.0. .OR. IHAZE.EQ.0)THEN
              AHOL1=HHOL
          ELSE
              AHOL1=HHAZE(ITYAER)
          ENDIF
          IF(AHAST(K).EQ.0.)THEN
              AHOL2=AHOL1
          ELSEIF(DENSTY(66,K).GT.0. .OR. DENSTY(3,K).GT.0.)THEN
              AHOL2=HHOL
          ELSE
              AHOL2='USER DEFINED        '
          ENDIF
          IF(ICLD.NE.0)THEN
              IF(ZGN(K).GT.2.00001D0)THEN
                  IF(PPROF(K).GT.9.9995)THEN
                      WRITE(IPR,'(F10.5,F10.3,2F8.2,1P,3E10.3,1X,3A)')  &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),WH2O(K),      &
     &                  DENSTY(66,K),DENSTY(3,K),AHOL1,AHOL2,           &
     &                  HSEASN(ISEA1)(1:13)
                  ELSE
                      WRITE(IPR,'(F10.5,F10.6,2F8.2,1P,3E10.3,1X,3A)')  &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),WH2O(K),      &
     &                  DENSTY(66,K),DENSTY(3,K),AHOL1,AHOL2,           &
     &                  HSEASN(ISEA1)(1:13)
                  ENDIF
              ELSE
                  IF(PPROF(K).GT.9.9995)THEN
                      WRITE(IPR,'(F10.5,F10.3,2F8.2,1P,3E10.3,1X,2A)')  &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),WH2O(K),      &
     &                  DENSTY(66,K),DENSTY(3,K),AHOL1,AHOL2
                  ELSE
                      WRITE(IPR,'(F10.5,F10.6,2F8.2,1P,3E10.3,1X,2A)')  &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),WH2O(K),      &
     &                  DENSTY(66,K),DENSTY(3,K),AHOL1,AHOL2
                  ENDIF
              ENDIF
          ELSE
              IF(ZGN(K).GT.2.00001D0)THEN
                  IF(PPROF(K).GT.9.9995)THEN
                      WRITE(IPR,'(F10.5,F10.3,2F8.2,1P,E10.3,1X,3A)')   &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),WH2O(K),      &
     &                  AHOL1,AHOL2,HSEASN(ISEA1)(1:13)
                  ELSE
                      WRITE(IPR,'(F10.5,F10.6,2F8.2,1P,E10.3,1X,3A)')   &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),WH2O(K),      &
     &                  AHOL1,AHOL2,HSEASN(ISEA1)(1:13)
                  ENDIF
              ELSE
                  IF(PPROF(K).GT.9.9995)THEN
                      WRITE(IPR,'(F10.5,F10.3,2F8.2,1P,E10.3,1X,2A)')   &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),              &
     &                  WH2O(K),AHOL1,AHOL2
                  ELSE
                      WRITE(IPR,'(F10.5,F10.6,2F8.2,1P,E10.3,1X,2A)')   &
     &                  ZM(K),PPROF(K),TPROF(K),RELHUM(K),              &
     &                  WH2O(K),AHOL1,AHOL2
                  ENDIF
              ENDIF
          ENDIF
      ENDDO
      RETURN
      END
