      SUBROUTINE CD2(LAPLUS,LMODEL,ICH,LNFLRT,FLRT)

!     PROCESS CARD2 INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
!       NL       NUMBER OF LEVELS IN MODEL ATMOSPHERES.
      INCLUDE 'PARAMS.h'
      INTEGER NL
      PARAMETER(NL=36)

!     ARGUMENTS:
!       LAPLUS   LOGICAL FLAG FOR THE AEROSOL A+ OPTION.
!       LMODEL   FLAG, .TRUE. IF MODEL ATMOSPHERE IS USED.
!       ICH      NUMERIC LABEL FOR AEROSOL MODEL
!       LNFLRT   LENGTH OF FILE ROOT NAME.
!       FLRT     FILE ROOT NAME.
      LOGICAL LAPLUS,LMODEL
      INTEGER ICH(4),LNFLRT
      CHARACTER FLRT*(NAMLEN-4)

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BASE.h'

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

!     /JM2/
!       CNOVAM   CHARACTER STRING FLAG FOR NOVAM OPTION (IF 'N' or 'n').
      CHARACTER APLUS*2,CNOVAM*1,ARUSS*3
      COMMON/JM2/APLUS,CNOVAM,ARUSS

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

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     /JM2APLUS/
      DOUBLE PRECISION ZAER11,ZAER12,ZAER21,ZAER22,                     &
     &  ZAER31,ZAER32,ZAER41,ZAER42
      REAL SCALE1,SCALE2,SCALE3,SCALE4
      COMMON/JM2APLUS/ZAER11,ZAER12,ZAER21,ZAER22,ZAER31,ZAER32,        &
     &  ZAER41,ZAER42,SCALE1,SCALE2,SCALE3,SCALE4

!     /USSPC/
      INTEGER NARSPC
      REAL VARSPC
      LOGICAL LARUSS
      COMMON/USSPC/NARSPC(4),VARSPC(4,NWAVLN),LARUSS

!     /COMNOV/
!       LNOVAM   LOGICAL FLAG, .TRUE. IF NOVAM AEROSOLS ARE USED.
      LOGICAL LNOVAM
      REAL EXTNOV(MNOV,MXWVLN),ABSNOV(MNOV,MXWVLN),                     &
     &  ASMNOV(MNOV,MXWVLN),WLNOV(MXWVLN)
      INTEGER NNOV,NWLNOV
      COMMON/COMNOV/LNOVAM,EXTNOV,ABSNOV,ASMNOV,WLNOV,NNOV,NWLNOV

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

!     /NOVDAT/
      DOUBLE PRECISION ALTNOV
      INTEGER NLNOV
      REAL RHNOV,PNOV,TNOV,DENNOV
      COMMON/NOVDAT/ALTNOV(MLNOV),NLNOV,RHNOV(MLNOV),                   &
     &  PNOV(MLNOV),TNOV(MLNOV),DENNOV(MNOV,MLNOV)

!     /VSBD/
      REAL VSB
      COMMON/VSBD/VSB(10)

!     /DISRT/
!       UMU      MONOTONICALLY INCREASING LIST OF DISTINCT USER-PATH
!                COSINE POLAR ANGLES.
!       PHI      MONOTONICALLY INCREASING LIST OF DISTINCT RELATIVE
!                SOLAR AZIMUTH ANGLES [0 TO 180 DEG].
!       NSTR     NUMBER OF DISCRETE ORDINATE STREAMS.
!       NAZ      NUMBER OF DISORT AZIMUTH COMPONENTS.
!       N2GAUS   ORDER OF DOUBLE-GAUSS QUADRATURES.
!       NUMU     NUMBER OF DISTINCT USER LINE-OF-SIGHT POLAR ANGLES.
!       MAPUMU   MAPPING FROM LINE-OF-SIGHT INDEX TO UMU ARRAY ENTRY.
!       NPHI     NUMBER OF DISTINCT RELATIVE SOLAR AZIMUTH ANGLES.
!       MAPPHI   MAPPING FROM LINE-OF-SIGHT INDEX TO PHI ARRAY ENTRY.
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       DISAZM   LOGICAL FLAG, TRUE FOR DISORT WITH AZIMUTH DEPENDENCE.
!       DISALB   LOGICAL FLAG, TRUE FOR DISORT SPHERICAL ALBEDO OPTION.
!       LDISCL   LOGICAL FLAG, TRUE FOR ISAACS SCALED TO DISORT.
      REAL UMU,PHI
      INTEGER NSTR,NAZ,N2GAUS,NUMU,MAPUMU,NPHI,MAPPHI
      LOGICAL DIS,DISAZM,DISALB,LDISCL
      COMMON/DISRT/UMU(MXUMU),PHI(MAXPHI),NSTR,NAZ,N2GAUS,NUMU,         &
     &  MAPUMU(MLOS),NPHI,MAPPHI(MLOS),DIS,DISAZM,DISALB,LDISCL
      SAVE /DISRT/

!     /SAP/
!       NWVSAP   NUMBER OF AEROSOL SPECTRAL GRID POINTS.
!       NLGSAP   HIGHEST AEROSOL PHASE FUNCTION LEGENDRE MOMENT.
!       NANSAP   NUMBER OF AEROSOL PHASE FUNCTION ANGULAR POINTS.
!       LEVSAP   NUMBER OF AEROSOL ATMOSPHERIC LEVELS.
!       ANGSAP   AEROSOL PHASE FUNCTION ANGULAR GRID [DEGREES].
!       COSSAP   COSINE OF THE AEROSOL PHASE FUNCTION ANGLES.
!       WAVSAP   AEROSOL SPECTRAL GRID [MICRONS].
!       LEGSAP   AEROSOL LEGENDRE COEFFICIENTS.
!       PFSAP    AEROSOL PHASE FUNCTION [SR-1].
!       LOSSAP   LOGICAL FLAG, TRUE FOR LINE-OF-SIGHT PATH.
      INTEGER NWVSAP,NLGSAP,NANSAP,LEVSAP
      DOUBLE PRECISION ANGSAP,COSSAP,PFSAP
      REAL WAVSAP,LEGSAP
      LOGICAL LOSSAP
      COMMON/SAP/NWVSAP,NLGSAP,NANSAP,LEVSAP,ANGSAP(MANSAP),            &
     &  COSSAP(MANSAP),WAVSAP(MWVSAP),LEGSAP(MLGSAP,MWVSAP,LAYDIM),     &
     &  PFSAP(MANSAP,MWVSAP,LAYDIM),LOSSAP

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     FUNCTIONS:
!       GETVIS   CONVERTS VERTICAL AEROSOL OPTICAL DEPTH TO VISIBILITY.
      REAL GETVIS

!     LOCAL VARIABLES:
!       LOPEN    LOGICAL FLAG, .TRUE. IF FILE IS OPEN.
!       LEXIST   LOGICAL FLAG, .TRUE. IF FILE EXISTS.
!       IOS      RESULT OF IOSTAT TEST.
!       IANSAP   AEROSOL PHASE FUNCTION ANGULAR GRID INDEX.
!       IANGM1   IANSAP MINUS 1.
!       AOD      VERTICAL AEROSOL OPTICAL DEPTH.
      LOGICAL LOPEN,LEXIST
      INTEGER IOS,IANSAP,IANGM1
      REAL AOD

!     SAVED COMMONS:
      SAVE /VSBD/

!     READ CARD2:
      IF(LJMASS)THEN
          CALL INITCD('CARD2')
      ELSE
          READ(IRD,'(A2,I3,A1,I4,A3,I2,3I5,5F10.0)')                    &
     &      APLUS,IHAZE,CNOVAM,ISEASN,ARUSS,                            &
     &      IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT
          WRITE(IPR,'(/A,A2,I3,A1,I4,A3,I2,3I5,5F10.5)')                &
     &      ' CARD 2  *****',APLUS,IHAZE,CNOVAM,ISEASN,ARUSS,           &
     &      IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT

          !WRITE(*,*) 'APLUS = ',APLUS
          !WRITE(*,*) 'IHAZE = ',IHAZE
          !WRITE(*,*) 'CNOVAM = ',CNOVAM
          !WRITE(*,*) 'ISEASN = ',ISEASN
          !WRITE(*,*) 'ARUSS = ',ARUSS
          !WRITE(*,*) 'IVULCN = ',IVULCN
          !WRITE(*,*) 'ICSTL = ',ICSTL
          !WRITE(*,*) 'ICLD = ',ICLD
          !WRITE(*,*) 'IVSA = ',IVSA
          !WRITE(*,*) 'VIS = ',VIS
          !WRITE(*,*) 'WSS = ',WSS
          !WRITE(*,*) 'WHH = ',WHH
          !WRITE(*,*) 'RAINRT = ',RAINRT
          !WRITE(*,*) 'GNDALT = ',GNDALT
      ENDIF

!     ARUSS='USS' MEANS USER-SUPPLED SPECTRAL DATA
!     (THIS IS NOT THE SWITCH FOR PHASE FUNCTIONS)
!     ARUSS='SAP' MEANS SPECTRAL AEROSOL PROFILES
      CALL UPCASE(ARUSS)
      LARUSS=ARUSS.EQ.'USS'
      LSAP=ARUSS.EQ.'SAP'

!     SAP OPTION AND MODEL ATMOSPHERE INCOMPATIBILITY CHECK
      IF(LSAP .AND. LMODEL)THEN
          WRITE(IPR,'(/A,/20X,A)')' Warning from CD2:  The Spectral'    &
     &      //' Aerosol Profiles (SAP) option cannot be used in',       &
     &      'conjunction with a model atmosphere.'                      &
     &      //'  The SAP option is being turned off.'
          LSAP=.FALSE.
      ENDIF

!     CHECK FOR SPECTRAL AEROSOL PROFILES DATA:
      IF(LSAP)THEN

!         OPEN SPECTRAL AEROSOL PROFILES FILE AND READ HEADER:
          INQUIRE(UNIT=ISAP,OPENED=LOPEN)

          IF(.NOT.LOPEN)THEN
              IF(LNFLRT.GT.0)THEN
                  INQUIRE(FILE=FLRT(1:LNFLRT)//'.sap',EXIST=LEXIST)
                  IF(LEXIST)THEN
                      CALL OPNFL(ISAP,0,FLRT(1:LNFLRT)//'.sap',         &
     &                  'OLD','FORMATTED','CD2')
                  ELSE
                      WRITE(IPR,'(/A,/(20X,A))')                        &
     &                  ' Warning from CD2:  The Spectral Aerosol'      &
     &                  //' Profiles (SAP) option is on, but file',     &
     &                  FLRT(1:LNFLRT)//'.sap','does not exist.'        &
     &                  //'  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ENDIF
              ELSE
                  INQUIRE(FILE='SpecAerProf.dat',EXIST=LEXIST)
                  IF(LEXIST)THEN
                      CALL OPNFL(ISAP,0,'SpecAerProf.dat',              &
     &                  'OLD','FORMATTED','CD2')
                  ELSE
                      WRITE(IPR,'(/A,/20X,A)')' Warning from CD2:  The' &
     &                  //' Spectral Aerosol Profiles (SAP) option is', &
     &                  'on, but file SpecAerProf.dat does not exist.'  &
     &                  //'  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ENDIF
              ENDIF
          ENDIF

!         READ AND CHECK HEADER:
          IF(LSAP)THEN
              READ(ISAP,*,IOSTAT=IOS)NWVSAP,NLGSAP,NANSAP
              IF(IOS.NE.0)THEN
                  WRITE(IPR,'(/A,/(20X,A))')' Warning from CD2:  The'   &
     &              //' Spectral Aerosol Profiles (SAP) option is on,', &
     &              'but header of aerosol data file could not read.'   &
     &              //'  The SAP option is being turned off.'
                  LSAP=.FALSE.
              ELSE
                  IF(NWVSAP.LT.1)THEN
                      WRITE(IPR,'(/A,/(20X,A))')                        &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the number', &
     &                  'of aerosol spectral points is not positive.'   &
     &                  //'  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ELSEIF(NWVSAP.GT.MWVSAP)THEN
                      WRITE(IPR,'(/A,2(/20X,A),I4,A))')                 &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the number', &
     &                  'of aerosol spectral points is too large;'      &
     &                  //' parameter MWVSAP in file PARAMS.h',         &
     &                  'must be increased to',NWVSAP,                  &
     &                  '  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ENDIF
                  IF(NLGSAP.LT.1)THEN
                      WRITE(IPR,'(/A,/(20X,A))')                        &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the number', &
     &                  'of Legendre moments is not positive.'          &
     &                  //'  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ELSEIF(DIS .AND. NLGSAP.LT.NSTR)THEN
                      WRITE(IPR,'(/A,/(20X,A))')                        &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the number', &
     &                  'of Legendre moments is less than'//            &
     &                  ' the number of DISORT streams.',               &
     &                  'The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ELSEIF(NLGSAP.GT.MLGSAP)THEN
                      WRITE(IPR,'(/A,/20X,A,/20X,I4,A)')                &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the number', &
     &                  'of Legendre moments is too large; parameter'// &
     &                  ' MLGSAP in file PARAMS.h must be increased to',&
     &                  NLGSAP,'  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ENDIF
                  IF(NANSAP.LT.2)THEN
                      WRITE(IPR,'(/A,/20X,A)')                          &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the number', &
     &                  'of aerosol phase function angles is less than' &
     &                  //' 2.  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ELSEIF(NANSAP.GT.MANSAP)THEN
                      WRITE(IPR,'(/A,2(/20X,A),I4,A)')                  &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the number', &
     &                  'of aerosol phase function angles is'//         &
     &                  ' too large; parameter MANSAP in',              &
     &                  'PARAMS.h must be increased to',NANSAP,         &
     &                  '  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ENDIF
              ENDIF
              IF(LSAP)THEN
                  READ(ISAP,*,IOSTAT=IOS)                               &
     &              (ANGSAP(IANSAP),IANSAP=1,NANSAP)
                  IF(IOS.NE.0)THEN
                      WRITE(IPR,'(/A,/(20X,A))')' Warning from CD2:  '//&
     &                  'The Spectral Aerosol Profiles (SAP) option is',&
     &                  'on, but phase function angles could not read.',&
     &                  'The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ELSEIF(ANGSAP(1).NE.0.D0                              &
     &              .OR. ANGSAP(NANSAP).NE.180.D0)THEN
                      WRITE(IPR,'(/A,/(20X,A))')                        &
     &                  ' Warning from CD2:  The Spectral Aerosol'//    &
     &                  ' Profiles (SAP) option is on, but the first',  &
     &                  'and last phase function angles were not 0 and' &
     &                  //' 180.  The SAP option is being turned off.'
                      LSAP=.FALSE.
                  ELSE
                      IANGM1=1
                      COSSAP(1)=1.D0
                      DO IANSAP=2,NANSAP
                          IF(ANGSAP(IANSAP).LE.ANGSAP(IANGM1))THEN
                              WRITE(IPR,'(2A,/(20X,A,I5,A,F15.10))')    &
     &                          ' Warning from CD2:  The Spectral'//    &
     &                          ' Aerosol Profiles (SAP) option is on,',&
     &                          ' but the aerosol phase angles'//       &
     &                          ' are not monotonically increasing',    &
     &                          'Angle',IANGM1,' (deg):',ANGSAP(IANGM1),&
     &                          'Angle',IANSAP,' (deg):',ANGSAP(IANSAP),&
     &                          'The SAP option is being turned off.'
                              LSAP=.FALSE.
                          ENDIF
                          COSSAP(IANSAP)=COS(ANGSAP(IANSAP)/DPDEG)
                          IANGM1=IANSAP
                      ENDDO
                      COSSAP(NANSAP)=-1.D0
                  ENDIF
              ENDIF
          ENDIF
      ENDIF

!     *****THE NEW AEROSOL IMPROVEMENTS (A+ SCHEME) (BEGIN)*****
!     THE A+ SCHEME GENERALIZES THE MODTRAN CURRENT AEROSOL MODELS.
!     THIS IS INDEPENDENT OF NOVAM RELATED MODIFICATIONS/IMPROVEMENTS.
!     IF IHAZE > 0 (I.E. AEROSOLS ARE INVOKED) AND APLUS='A+',
!     FOUR PAIRS OF ALTITUDE ARE READ IN THE FOLLOWING CARD.
!     SEE THE ROUTINE APRFNU.F.
      LAPLUS=APLUS.EQ.'A+'
      IF(LAPLUS)THEN
          IF(IHAZE.LT.0 .OR. IHAZE.EQ.7 .OR. LSAP .OR. ICLD.EQ.11)THEN
              WRITE(IPR,'(//A,/(23X,A))')                               &
     &          ' WARNING from CD2:  Inputs IHAZE<0, IHAZE=7,'          &
     &          //' LSAP, and ICLD=11 are incompatible',                &
     &          'with the "A+" option and "A+" is being turned off.',   &
     &          'An error will occur if tape5 includes CARD 2A+.'
              LAPLUS=.FALSE.
          ENDIF

!         IF(LAPLUS)READ CARD 2A+
          IF(LAPLUS)THEN
              IF(LJMASS)THEN
                  !CALL INITCD('CARD2APLUS')
              ELSE
                  READ(IRD,'((3(1X,F9.0),20X,3(1X,F9.0)))')             &
     &              ZAER11,ZAER12,SCALE1,ZAER21,ZAER22,SCALE2,          &
     &              ZAER31,ZAER32,SCALE3,ZAER41,ZAER42,SCALE4
              ENDIF

!             CHECK FOR TROPOSPHERIC AEROSOL:
              IF(IHAZE.EQ.6)THEN
                  WRITE(IPR,'(/2A,/(15X,A))')' ***  WARNING: ',         &
     &              ' Aerosol regions 1 (nominally, 0-2 km) and',       &
     &              ' 2 (nominally, 2-10 km) are combined when',        &
     &              ' the Tropospheric aerosol model (IHAZE=6)',        &
     &              ' is selected.  From CARD 2A+, the following',      &
     &              ' assignments define the tropospheric region:',     &
     &              '   LOWER BOUNDING ALTITUDE = MIN(ZAER11,ZAER21)',  &
     &              '   UPPER BOUNDING ALTITUDE = MAX(ZAER12,ZAER22)',  &
     &              '   AEROSOL SCALE FACTOR    = MAX(SCALE1,SCALE2)'
                  ZAER11=MIN(ZAER11,ZAER21)
                  ZAER12=MAX(ZAER12,ZAER22)
                  SCALE1=MAX(SCALE1,SCALE2)
                  ZAER21=0.D0
                  ZAER22=0.D0
                  SCALE2=1.
              ENDIF
              IF(.NOT.LJMASS)WRITE(IPR,'(A,12(1X,F9.4))')' CARD 2A+ **',&
     &          ZAER11,ZAER12,SCALE1,ZAER21,ZAER22,SCALE2,              &
     &          ZAER31,ZAER32,SCALE3,ZAER41,ZAER42,SCALE4
          ENDIF
      ENDIF

!     *****THE NEW AEROSOL IMPROVEMENTS (A+ SCHEME) (END)*****

!     CNOVAM = "N" TRIGGERS NOVAM
      LNOVAM=CNOVAM.EQ.'N' .OR. CNOVAM.EQ.'n'

!     NOVAM WITH A+ INCOMPATIBILITY CHECK:
      IF(LNOVAM .AND. LAPLUS)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' Error in CD2:  NOVAM AND A+ OPTIONS ARE INCOMPATIBLE.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' ERROR:  NOVAM AND A+ OPTIONS ARE INCOMPATIBLE.'
      ENDIF

!     NOVAM WITH SAP INCOMPATIBILITY CHECK:
      IF(LNOVAM .AND. LSAP)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' Error in CD2:  NOVAM AND SAP OPTIONS ARE INCOMPATIBLE.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' ERROR:  NOVAM AND SAP OPTIONS ARE INCOMPATIBLE.'
      ENDIF

!     CHECK IF IHAZE OR ICLD NEED TO BE RESET.
!       IF IHAZE < 0, THEN NO AEROSOLS BUT CLOUDS IF 0 < ICLD < 11
!       IF IHAZE = 0 AND CNOVAM .NE. "N", NO AEROSOLS AND NO CLOUDS
!       IF IHAZE > 0, THEN AEROSOLS AND, IF ICLD > 0, CLOUDS
      IF(IHAZE.EQ.0 .AND. ICLD.NE.0)THEN

!         RESET ICLD TO ZERO - NO CLOUDS
          WRITE(IPR,'(/2A,I3,A)')' WARNING:  INPUT ICLD IS BEING',      &
     &      ' RESET FROM',ICLD,' TO 0 SINCE IHAZE EQUALS 0.'
          ICLD=0
      ELSEIF(IHAZE.LT.0)THEN

!         FOR INTERNAL USE, SET IHAZE TO ZERO (CLOUDS WILL BE
!         INCLUDED IF ICLD IS BETWEEN 1 AND 10, INCLUSIVE).
          IHAZE=0
      ENDIF

!     CHECK GROUND ALTITUDE.
      IF(.NOT.LJMASS .AND. GNDALT.NE.0.D0)                              &
     &  WRITE(IPR,'(/A,F10.5)')'   GNDALT =',GNDALT
      IF(GNDALT.GE.6.D0)THEN
          WRITE(IPR,'(/2A,F10.5)')' WARNING from CD2:  GNDALT',         &
     &      ' (> 6 KM) IS BEING RESET TO 0 KM; GNDALT WAS',GNDALT
          GNDALT=0.D0
      ENDIF

!     CHECK FOR WINTER:
      IF((MODEL.EQ.3 .OR. MODEL.EQ.5) .AND. ISEASN.EQ.0)ISEASN=2

!     CHECK FOR UNDEFINED VISIBILITY:
      IF(VIS.LT.0.)THEN

!         CONVERT VERTICAL AEROSOL OPTICAL DEPTH (-VIS) TO VISIBILITY:
          AOD=-VIS
          VIS=GETVIS(AOD,GNDALT,ISEASN.NE.2,IVULCN,IPR)
          WRITE(IPR,'(/(A,F10.5,A))')' The 550nm vertical'//            &
     &      ' optical depth (',AOD,') corresponds to a surface',        &
     &      ' meteorological range (visibility) of',VIS,'km.'
      ENDIF
      IF(VIS.LE.0. .AND. IHAZE.GT.0)THEN

!         AEROSOL MODEL DEPENDENT VISIBILITY:
          VIS=VSB(IHAZE)
      ENDIF

!     DEFINE AEROSOL BOUNDARY LEVEL ARRAY, ICH:
      IF(LMODEL)THEN
          ML=NL
          IF(IVSA.EQ.1 .AND. IHAZE.EQ.3)                                &
     &      CALL MARINE(VIS,MODEL,0.,WSS,WHH,ICSTL,EXTC,ABSC,1)
          ICH(1)=MIN(IHAZE,1)
          ICH(2)=6
          ICH(3)=9+MIN(IVULCN,1)
      ELSE
          ICH(1)=1
          ICH(2)=0
          ICH(3)=10
      ENDIF
      ICH(4)=18
      IF(ICLD.EQ.11)THEN
          ICH(4)=ICH(3)
          ICH(3)=ICH(2)
          ICH(2)=ICLD
      ENDIF

!     NOVAM AEROSOL MODEL:
      IF(LNOVAM)CALL NOVAER(IHAZE,VIS,ALTNOV,RHNOV,PNOV,TNOV,DENNOV)
      IF(.NOT.LJMASS .AND. RAINRT.NE.0.)WRITE(IPR,'(/A,F9.3,A)')        &
     &  '   RAIN MODEL CALLED, RAIN RATE = ',RAINRT,' MM/HR'

!     RETURN TO DRIVER:
      RETURN
      END
