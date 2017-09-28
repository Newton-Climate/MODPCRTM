      SUBROUTINE CD1A(LENN,ASTMC,ASTMX,ASTMO,AERRH,LNFILT,LENDAT,       &
     &  KNTRVL,LASTM,LDASTM,NSSALB,RHASYM,KPRINT,LABEL)

!     PROCESS CARD1A INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       LENN     PATH LENGTH SWITCH (0=SHORT, 1=LONG).
!       ASTMC    AEROSOL ANGSTROM LAW COEFFICIENT @ 550 NM.
!       ASTMX    AEROSOL ANGSTROM LAW EXPONENT.
!       ASTMO    AEROSOL ANGSTROM LAW OFFSET.
!       AERRH    BOUNDARY LAYER/TROPOSPHERIC AEROSOL RELATIVE HUMID [%].
!       LNFILT   LENGTH OF FILTER RESPONSE FUNCTION FILE NAME.
!       LENDAT   LENGTH OF DATDIR STRING.
!       KNTRVL   NUMBER OF SUB-INTERVALS USED IN CORRELATED-K ALGORITHM.
!       LASTM    LOGICAL FLAG, TRUE IS ANGSTROM PARAMETERS ARE INPUT.
!       LDASTM   LOGICAL IMPLYING THAT ASTMX IS A DELTA VALUE.
!       KPRINT   LOGICAL FLAG DICTATING K-DISTRIBUTION OUTPUT.
      REAL ASTMC,ASTMX,ASTMO,AERRH,RHASYM
      INTEGER LENN,LNFILT,LENDAT,KNTRVL,NSSALB
      LOGICAL LASTM,LDASTM(2),KPRINT
      CHARACTER LABEL*6

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     /JM1A1/
!       RESCHR   BAND MODEL RESOLUTION CHARACTER STRING.
!       DISSTR   CHARACTER STRING USED TO READ IN DISORT LOGICALS.
!       H2OSTR   VERTICAL WATER COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE WATER
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE WATER COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS A WATER COLUMN SCALING FACTOR).
!       O3STR    VERTICAL OZONE COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE OZONE
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE OZONE COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS AN OZONE COLUMN SCALING FACTOR).
!       USRSUN   USER-SPECIFIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       BMROOT   PREFIX OF MOLECULAR BAND MODEL PARAMETERS FILE.
!       FILTNM   NAME OF FILTER RESPONSE FUNCTION FILE.
!       H2OAER   FLAG, TRUE IF DEFAULT AEROSOL PROPERTIES ARE REVISED
!                BASED ON WATER COLUMN SCALING.
!       DATDIR   NAME OF THE MODTRAN DATA DIRECTORY.
      CHARACTER RESCHR*2,DISSTR*3,H2OSTR*10,O3STR*10,USRSUN*(NAMLEN),   &
     &  FILTNM*(NAMLEN),BMROOT*(NAMLEN),H2OAER*1,DATDIR*(NAMLEN-LENSUN)
      COMMON/JM1A1/RESCHR,DISSTR,H2OSTR,O3STR,USRSUN,                   &
     &  FILTNM,BMROOT,H2OAER,DATDIR

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

!     /JM1/
!       CODE     CODE INPUT FLAG: 'F' OR 'L' FOR LOWTRAN BAND MODEL.
!                                 'T' OR 'M' FOR MODTRAN BAND MODEL.
!                                 'C' OR 'K' FOR CORRELATED-K MODEL.
!       SPEED    COMPUTATIONAL SPEED FLAG USED WITH CK MODEL ['F' OR '2'
!                FOR FAST (12 K'S), 'M' OR '1' FOR MODERATE (17 K'S),
!                AND 'S' OR '0' FOR SLOW (33 K'S); 'S' IS DEFAULT]
!       SURREF   SURFACE REFLECTANCE CHARACTER STRING.
!                IF FIRST NON-BLANK CHARACTER IS "B",
!                  BI-DIRECTIONAL REFLECTANCE DISTRIBUTION
!                  FUNCTION (BRDF) DATA IS READ IN.
!                IF FIRST NON-BLANK CHARACTER IS "L" OR "-"
!                  SURFACE IS MODELED AS A LAMBERTIAN REFLECTOR AND
!                  SPECTRAL ALBEDO IS READ FROM FILE "DATA/spec_alb.dat"
!                OTHERWISE, THE CHARACTER STRING IS ASSUMED TO CONTAIN
!                  A SPECTRALLY INDEPENDENT VALUE FOR SURFACE ALBEDO.
      CHARACTER CODE*1,SPEED*1,SURREF*7
      COMMON/JM1/CODE,SPEED,SURREF

!     /JM1A2/
!       SOLCON   SOLAR CONSTANT (ACTUAL VALUE IF POSITIVE,
!                                SCALE FACTOR IF NEGATIVE).
      LOGICAL LSUN
      REAL SFWHM,CO2MX,SOLCON
      COMMON/JM1A2/LSUN,SFWHM,CO2MX,SOLCON

!     /CO2MIX/
!       CO2RAT   RATIO OF USER DEFINED TO DEFAULT CO2 PROFILE.
      REAL CO2RAT
      COMMON/CO2MIX/CO2RAT

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

!     /SUNFLG/
!       LNFLRT   LENGTH OF I/O FILE ROOT NAME, 0 IF NO mod5root.in FILE.
!       LRDSUN   SOLAR DATA FLAG, TRUE IF IRRADIANCES IN COMMON BLOCK
!                /SOL01/ HAVE BEEN MODIFIED FROM THE BLOCK DATA.
      INTEGER LNFLRT
      LOGICAL LRDSUN
      COMMON/SUNFLG/LNFLRT,LRDSUN

!     /SUNNAM/
!       SUNFIL   NAME OF DEFAULT FILE CONTAINING SOLAR IRRADIANCE DATA.
!       FLRT     ROOT NAME FOR ALL I/O FILES.
      CHARACTER FLRT*(NAMLEN-4),SUNFIL*(LENSUN)
      COMMON/SUNNAM/FLRT,SUNFIL

!     /PATH/
!       PTHCOS   COSINE OF PATH ZENITH AT PATH BOUNDARIES.
!       PTHZEN   PATH ZENITH AT PATH BOUNDARIES [DEG].
!       PTHECA   SENSOR TO PATH EARTH CENTER ANGLE [DEG].
!       PTHALT   ALTITUDES AT PATH BOUNDARIES [KM].
!       PTH_MS   ALTITUDES AT PATH BOUNDARIES FOR THE MS PATH.
!       PTHSEG   PATH SEGMENT LENGTH [KM].
!       PTHRNG   SENSOR TO PATH BOUNDARY RANGE [KM].
!       JMAX     NUMBER OF DISTINCT LOS PATH SEGMENT ENDPOINT ALTITUDES.
!       IKHMIN   PATH BOUNDARY INDEX OF PATH MINIMUM ALTITUDE.
!       IKHMAX   PATH BOUNDARY INDEX OF PATH MAXIMUM ALTITUDE.
!       IKOUT    NUMBER OF PATH BOUNDARIES K DATA IS OUTPUT.
!       NTKDIS   RECORD NUMBER FOR K-DISTRIBUTION TRANSMITTANCE FILE.
!       NRKDIS   RECORD NUMBER FOR K-DISTRIBUTION RADIANCE FILE.
!       MAPPTH   MAPPING FROM PATH SEGMENT MIDPOINT TO VERTICAL LAYER.
!       IPTHHT   ALTITUDES (HEIGHTS) AT PATH BOUNDARIES [M].
!       LOWALT   VERTICAL LAYER BOUNDARY INDEX AT OR JUST BELOW PTHALT.
!       FACALT   ALTITUDE INTERPOLATION FRACTION FOR PTHALT
!       PATH_T   TEMPERATURE AT PATH BOUNDARIES [K].
!       PATH_P   PRESSURE AT PATH BOUNDARIES [ATM].
!       PTHRH    RELATIVE HUMIDITY AT PATH BOUNDARIES [K].
!       LSSGEO   LOGICAL FLAG, .TRUE. FOR SOLAR PATHS.
!       LTANMX   LOGICAL FLAG, .TRUE. IF PATH HAS A TANGENT MAXIMUM.
      DOUBLE PRECISION PTHCOS,PTHZEN,PTHECA,PTHALT,PTH_MS,PTHSEG,PTHRNG
      INTEGER JMAX,IKHMIN,IKHMAX,IKOUT,NTKDIS,NRKDIS,MAPPTH,            &
     &  IPTHHT,LOWALT
      REAL FACALT,PATH_T,PATH_P,PTHRH
      LOGICAL LSSGEO,LTANMX
      COMMON/PATH/PTHCOS(0:LAYTWO),PTHZEN(0:LAYTWO),PTHECA(0:LAYTWO),   &
     &  PTHALT(0:LAYTWO,1:MLOS),PTH_MS(0:LAYDIM),PTHSEG(LAYTWO),        &
     &  PTHRNG(0:LAYTWO,1:MLOS),JMAX,IKHMIN(MLOS),IKHMAX(MLOS),         &
     &  IKOUT(MLOS),MAPPTH(LAYTWO,1:MLOS),IPTHHT(0:LAYTWO),NTKDIS,      &
     &  NRKDIS,LOWALT(0:LAYTWO,1:MLOS),FACALT(0:LAYTWO,1:MLOS),         &
     &  PATH_T(0:LAYTWO,1:MLOS),PATH_P(0:LAYTWO,1:MLOS),                &
     &  PTHRH(0:LAYTWO,1:MLOS),LSSGEO,LTANMX

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

!     /NAMEX/
!       CNAMEX   NAME OF CROSS-SECTION (X) SPECIES.
      CHARACTER CNAMEX*8
      COMMON/NAMEX/CNAMEX(NMOLX)

!     SAVED COMMONS:
      SAVE /SUNFLG/,/SUNNAM/,/NAMEX/

!     FUNCTIONS:
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
!       IRECLN   RETURNS RECORD LENGTH FOR A DIRECT ACCESS FILE.
      INTEGER LENSTR,IRECLN

!     LOCAL VARIABLES:
!       CSTR     THIRD DIGIT OF NSTR INPUT (BACKWARD COMPATIBLE).
!       C_PROF   CHARACTER STRING USED FOR SCALING DEFAULT PROFILES.
!       CDTDIR   CHARACTER SRING USED FOR INPUTTING DATA DIRECTORY NAME.
!       CDASTM   CHARACTER STRING USED TO DEFINE LOGICAL FLAGS LDASTM.
!       CKNAME   CORRELATED-K DISTRIBUTIONS FILE FULL PATH NAME.
!       DFTSUN   DEFAULT TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       NSPEED   INTEGER CORRESPONDING TO SPEED.
!       L_USUN   LENGTH OF FILE NAME USRSUN.
!       IOS      I/O STATUS OF A FORMATTED READ.
!       LNKDIS   RECORD LENGTH OF K-DISTRIBUTION DEPENDENT BINARY FILES.
!       IMOL     LOOP INDEX FOR MOLECULAR SPECIES.
!       I_PROF   NUMERICAL VALUE OF C_PROF INPUT.
!       LOPEN    OPEN FILE FLAG
      CHARACTER*1 CSTR,C_PROF,CDTDIR,CDASTM
      CHARACTER*(NAMLEN) CKNAME,DFTSUN
      INTEGER NSPEED,L_USUN,IOS,LNKDIS,IMOL,I_PROF
      LOGICAL LOPEN

!     DATA:
!       CNAME    NAME OF UNIFORMLY MIXED MOLECULAR SPECIES.
!       YNAM16   NAME OF TRACE MOLECULAR SPECIES.
!       FRMT     FORMAT STRING.
!       CORKNM   NAME OF CORRELATED-K BINARY FILE.
!       CKNMSV   OLD CORRELATED-K DISTRIBUTIONS FILE NAME.
!       NSPDSV   SAVED VALUE OF VARIABLE NSPEED.
      CHARACTER*8 CNAME(4:NMOL+1),YNAM16(6:21)
      CHARACTER FRMT*9,CORKNM*10,CKNMSV*(NAMLEN)
      INTEGER NSPDSV,KNTRSV
      SAVE CNAME,YNAM16,FRMT,CORKNM,CKNMSV,NSPDSV,KNTRSV
      DATA CNAME,YNAM16/  '     N2O','      CO','     CH4','      O2',  &
     &         '      NO','     SO2','     NO2','     NH3','    HNO3',  &
     &         '      N2','      OH','      HF','     HCl','     HBr',  &
     &         '      HI','     ClO','     OCS','    H2CO','    HOCl',  &
     &         '      N2','     HCN','   CH3Cl','    H2O2','    C2H2',  &
     &         '    C2H6','     PH3'/,FRMT/'((A    ))'/,                &
     &  CORKNM/'CORK00.BIN'/,CKNMSV/' '/,NSPDSV,KNTRSV/2*0/

!     INPUT CARD1A:
      IF(LJMASS)THEN
          CALL INITCD('CARD1A')
          L_UMIX=.FALSE.
          L_XSEC=.FALSE.
          L_TRAC=.FALSE.
      ELSE
          USRSUN=' '
          READ(IRD,'(A3,I2,A1,F4.0,F10.0,2A10,2A1,4(1X,A1),F10.0,A1,    &
     &      F9.0,3F10.0,I10)',IOSTAT=IOS)DISSTR,NSTR,CSTR,SFWHM,CO2MX,  &
     &      H2OSTR,O3STR,C_PROF,USRSUN(1:1),BMROOT(1:1),FILTNM(1:1),    &
     &      H2OAER,CDTDIR,SOLCON,CDASTM,ASTMC,ASTMX,ASTMO,AERRH,NSSALB
          !WRITE(*,*) 'DISSTR = ',DISSTR
          !WRITE(*,*) 'NSTR = ',NSTR
          !WRITE(*,*) 'CSTR = ',CSTR
          !WRITE(*,*) 'SFWHM = ',SFWHM
          !WRITE(*,*) 'CO2MX  = ',CO2MX
          !WRITE(*,*) 'H2OSTR = ',H2OSTR
          !WRITE(*,*) 'O3STR = ',O3STR
          !WRITE(*,*) 'C_PROF = ',C_PROF
          !WRITE(*,*) 'USRSUN = ',USRSUN
          !WRITE(*,*) 'BMROOT = ',BMROOT
          !WRITE(*,*) 'FILTNM = ',FILTNM
          !WRITE(*,*) 'H2OAER = ',H2OAER
          !WRITE(*,*) 'CDTDIR = ',CDTDIR
          !WRITE(*,*) 'SOLCON = ',SOLCON
          !WRITE(*,*) 'CDASTM = ',CDASTM
          !WRITE(*,*) 'ASTMC =  ',ASTMC
          !WRITE(*,*) 'ASTMX = ',ASTMX
          !WRITE(*,*) 'ASTMO = ',ASTMO
          !WRITE(*,*) 'AERRH = ',AERRH
          !WRITE(*,*) 'NSSALB = ',NSSALB

          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/A)')'Error in CD1A:  Unable to'//            &
     &          ' successfully read CARD 1A from MODTRAN input file'
              STOP 'Error in CD1A:  Unable to read CARD 1A!'
          ENDIF

!         CHECK FOR EXPANDED NSTR INPUT:
          IF(INDEX('0123456789',CSTR).GT.0)                             &
     &      NSTR=10*NSTR+INDEX('0123456789',CSTR)-1
          LSUN=ABS(SFWHM).GT..0
          I_PROF=INDEX('1234567',C_PROF)
          L_UMIX=MOD(I_PROF,2).EQ.1
          L_XSEC=I_PROF/2.EQ.1 .OR. I_PROF/2.EQ.3
          L_TRAC=I_PROF.GE.4
      ENDIF
      IF(CO2MX.LE.0.)CO2MX=330.
      CO2RAT=CO2MX/330.

!     FOR BACKWARD COMPATIBILITY WITH OLDER MODTRAN INPUT FILES,
!     THE DISORT LOGICALS ARE READ INTO A CHARACTER STRING.
      CALL UPCASE(DISSTR)

!     DISORT COMPUTATION IS USED TO SCALE ISAACS IF LDSICL = .TRUE.
      IF(IMULT.EQ.0)THEN
          LDISCL=.FALSE.
          DIS=.FALSE.
      ELSE
          LDISCL=DISSTR(1:1).EQ.'S'
          DIS=LDISCL .OR. DISSTR(1:1).EQ.'T'
      ENDIF
      IF(DIS)THEN
          DISAZM=IEMSCT.EQ.2 .AND. DISSTR(2:2).EQ.'T'
          DISALB=DISSTR(3:3).EQ.'T'
          IF(DISALB)THEN
              DISALB=IEMSCT.EQ.2
              IF(.NOT.DISALB)WRITE(IPR,'(A,/19X,A)')                    &
     &          'Warning from CD1A:  The DISORT Spherical Albedo'//     &
     &          ' / Diffuse Transmittance',' Option was turned off'//   &
     &          ' because solar scattering is not being computed.'
          ENDIF
          NSTR=ABS(NSTR)
          IF(NSTR.GT.MXCMU)THEN
              WRITE(IPR,'(A,2(I4,A),/19X,A)')                           &
     &          'Warning from CD1A:  The requested number of DISORT'    &
     &          //' streams (NSTR =',NSTR,') exceeds the maximum'       &
     &          //' allowed (MXCMU =',MXCMU,').  Input NSTR',           &
     &          'is being reduced to MXCMU.  Increase parameter'//      &
     &          ' MXCMU in file PARAMS.h to use additional streams.'
              NSTR=MXCMU
          ENDIF
          NAZ=NSTR-1
          N2GAUS=NSTR/2
      ELSE
          DISAZM=.FALSE.
          DISALB=.FALSE.
          NSTR=1
          NAZ=-1
          N2GAUS=0
      ENDIF
      IF(MODEL.EQ.0)THEN

!         DO NOT SCALE PROFILE FOR CONSTANT PRESSURE PATH:
          LENN=0
          H2OSTR='          '
          O3STR='          '
          H2OAER=' '
      ELSEIF(H2OSTR.EQ.'          ')THEN
          H2OAER=' '
      ELSEIF(H2OAER.EQ.'t')THEN
          H2OAER='T'
      ENDIF
      IF(.NOT.LJMASS)THEN
          WRITE(IPR,'(/A,3L1,I3,F4.1,F10.5,2A10,2A1,4(1X,A1),           &
     &      F10.3,A1,F9.4,3F10.4,I10)')' CARD 1A *****',DIS,            &
     &      DISAZM,DISALB,NSTR,ABS(SFWHM),CO2MX,H2OSTR,O3STR,           &
     &      C_PROF,USRSUN(1:1),BMROOT(1:1),FILTNM(1:1),H2OAER,          &
     &      CDTDIR,SOLCON,CDASTM,ASTMC,ASTMX,ASTMO,AERRH,NSSALB
          IF(LDISCL)WRITE(IPR,'(/A)')                                   &
     &      ' A spectral grid of DISORT calculations will be'           &
     &      //' used to scale Isaacs'' scattering results.'
      ENDIF
      LASTM=ASTMX.NE.0.
      IF(LASTM)THEN

!         LDASTM(1) IS .TRUE. IF THE BOUNDARY LAYER AEROSOL OPTICAL
!         PROPERTIES ARE TO BE PERTURBED; LDASTM(2) IS .TRUE. IF
!         THE TROPOSPHERIC AEROSOL PROPERTIES ARE TO BE PERTURBED.
          LDASTM(2)=CDASTM.EQ.'D' .OR. CDASTM.EQ.'d'                    &
     &      .OR. CDASTM.EQ.'T' .OR. CDASTM.EQ.'t'
          LDASTM(1)=CDASTM.EQ.'B' .OR. CDASTM.EQ.'b' .OR. LDASTM(2)
          IF(.NOT.LDASTM(1))THEN

!             ANGSTROM LAW:  EXT=ASTMO+ASTMC*(550nm/WAVLEN)**ASTMX
              IF(ASTMO+ASTMC.LE.0.)THEN
                  WRITE(IPR,'(/A,/21X,A)')' Warning from CD1A:  '//     &
     &              'Angstrom Law extinction is not positive at 550nm.',&
     &              'Angstrom Law option is being turned off.'
                  LASTM=.FALSE.
                  ASTMC=0.
              ELSE

!                 NORMALIZE THE 550NM EXTINCTION:
                  ASTMC=ASTMC/(ASTMO+ASTMC)
                  ASTMO=1.-ASTMC
                  WRITE(IPR,'(3(/A),3(F10.5,A))')' The BOUNDARY LAYER'  &
     &              //' and TROPOSPHERIC aerosol spectral extinction',  &
     &              ' is being defined with an Angstrom Law'//          &
     &              ' dependence:',' EXT(lambda)/EXT(550nm) =',         &
     &              ASTMO,' +',ASTMC,' * (550nm/lambda) ^',ASTMX
              ENDIF
          ELSEIF(LDASTM(2))THEN
              WRITE(IPR,'(3(/A),F10.5)')' The BOUNDARY LAYER and'       &
     &          //' TROPOSPHERIC aerosol spectral extinction',          &
     &          ' is being perturbed with an Angstrom Law dependence:', &
     &          ' EXT(lambda)'' = EXT(lambda) * (550nm/lambda) ^',ASTMX
          ELSE
              WRITE(IPR,'(3(/A),F10.5)')                                &
     &          ' The BOUNDARY LAYER aerosol spectral extinction',      &
     &          ' is being perturbed with an Angstrom Law dependence:', &
     &          ' EXT(lambda)'' = EXT(lambda) * (550nm/lambda) ^',ASTMX
          ENDIF
      ELSE
          ASTMC=0.
          LDASTM(1)=.FALSE.
          LDASTM(2)=.FALSE.
      ENDIF
      IF(AERRH.GT.0.)THEN
         IF(AERRH.GT.99.)AERRH=99.
         IF(.NOT.LJMASS)WRITE(IPR,'(/A,F5.1,A)')' ***  a relative'      &
     &     //' humidity of',AERRH,'% will be used to define the'        &
     &     //' boundary layer and tropospheric aerosols  ***'
      ENDIF

      IF(.NOT.LJMASS)THEN

!         DEFINE INPUT FORMAT:
          WRITE(FRMT(4:7),'(I4.4)')NAMLEN

!         CARD 1A1:
          IF(USRSUN(1:1).EQ.'T' .OR. USRSUN(1:1).EQ.'t')THEN

!             USRSUN IS READ IN AND USED EVEN IF SFWHM=0.
              READ(IRD,FRMT)USRSUN
              LSUN=.TRUE.
          ELSEIF(USRSUN(1:1).EQ.'1' .OR. USRSUN(1:1).EQ.'2' .OR.        &
     &           USRSUN(1:1).EQ.'3' .OR. USRSUN(1:1).EQ.'4' .OR.        &
     &           USRSUN(1:1).EQ.'5' .OR. USRSUN(1:1).EQ.'6' .OR.        &
     &           USRSUN(1:1).EQ.'7')THEN

!             MODTRAN SUPPLIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE DATA
!             FILE READ IN AND USED EVEN IF SFWHM=0.
              LSUN=.TRUE.
          ELSE
              USRSUN=' '
          ENDIF

!         CARD 1A2:
          IF(BMROOT(1:1).EQ.'T' .OR. BMROOT(1:1).EQ.'t')THEN
              READ(IRD,FRMT)BMROOT
          ELSE
              BMROOT=' '
          ENDIF

!         CARD 1A3:
          IF(FILTNM(1:1).EQ.'T' .OR. FILTNM(1:1).EQ.'t')THEN
              READ(IRD,FRMT)FILTNM
          ELSE
              FILTNM=' '
          ENDIF

!         CARD 1A4:
          IF(CDTDIR.EQ.'T' .OR. CDTDIR.EQ.'t')THEN
              WRITE(FRMT(4:7),'(I4.4)')NAMLEN-LENSUN
              READ(IRD,FRMT)DATDIR
          ENDIF

!         CARD 1A5:
          IF(L_UMIX)THEN
              READ(IRD,'(10F5.0)',IOSTAT=IOS)S_UMIX
              IF(IOS.NE.0)THEN
                  WRITE(IPR,'(/A,/15X,A)')'Error in CD1A:  Unable to'   &
     &              //' successfully read profile scale factors for',   &
     &              ' uniformly mixed molecular species from CARD 1A5!'
                  STOP 'Error in CD1A:  Unable to read CARD 1A5!'
              ENDIF
              WRITE(IPR,'(/A)')'SCALE FACTORS FOR UNIFORMLY MIXED'      &
     &          //' MOLECULAR SPECIES DEFAULT PROFILES'
              DO IMOL=4,NMOL+1
                  IF(S_UMIX(IMOL).LT.0.)S_UMIX(IMOL)=1.
                  IF(ABS(S_UMIX(IMOL)-1).GT..00001)                     &
     &              WRITE(IPR,'(A,F9.5)')CNAME(IMOL),S_UMIX(IMOL)
              ENDDO
          ENDIF

!         CARD 1A6:
          IF(L_XSEC)THEN
              READ(IRD,'(15F5.0)',IOSTAT=IOS)S_XSEC
              IF(IOS.NE.0)THEN
                  WRITE(IPR,'(/A,/15X,A)')'Error in CD1A:  Unable to'   &
     &              //' successfully read profile scale factors for',   &
     &              ' cross-section molecular species from CARD 1A6!'
                  STOP 'Error in CD1A:  Unable to read CARD 1A6!'
              ENDIF
              WRITE(IPR,'(/A)')'SCALE FACTORS FOR CROSS-SECTION'        &
     &          //' MOLECULAR SPECIES DEFAULT PROFILES'
              DO IMOL=1,NMOLX
                  IF(S_XSEC(IMOL).LT.0.)S_XSEC(IMOL)=1.
                  IF(ABS(S_XSEC(IMOL)-1).GT..00001)                     &
     &              WRITE(IPR,'(A,F9.5)')CNAMEX(IMOL),S_XSEC(IMOL)
              ENDDO
          ENDIF

!         CARD 1A7:
          IF(L_TRAC)THEN
              READ(IRD,'(16F5.0)',IOSTAT=IOS)S_TRAC
              IF(IOS.NE.0)THEN
                  WRITE(IPR,'(/A,/15X,A)')'Error in CD1A:  Unable to'   &
     &              //' successfully read profile scale factors for',   &
     &              ' trace molecular species from CARD 1A7!'
                  STOP 'Error in CD1A:  Unable to read CARD 1A7!'
              ENDIF
              WRITE(IPR,'(/A)')'SCALE FACTORS FOR TRACE'                &
     &          //' MOLECULAR SPECIES DEFAULT PROFILES'
              DO IMOL=6,21
                  IF(S_TRAC(IMOL).LT.0.)S_TRAC(IMOL)=1.
                  IF(ABS(S_TRAC(IMOL)-1).GT..00001)                     &
     &              WRITE(IPR,'(A,F9.5)')YNAM16(IMOL),S_TRAC(IMOL)
              ENDDO
          ENDIF
      ENDIF
      IF(LJMASS)THEN
          LNFILT=0
      ELSE
          LNFILT=LENSTR(FILTNM)
      ENDIF
      LENDAT=LENSTR(DATDIR)
      CKNAME=' '

      CALL RSSALB(NSSALB,RHASYM)

!     OPEN MOLECULAR BAND MODEL PARAMETERS FILE AND READ HEADER:
      KNTRVL=1
      IF(MODTRN)THEN

!         OPEN BAND MODEL DATA FILES:
          CALL OPENBM(LENDAT,DATDIR(1:LENDAT),BMROOT,RESCHR,LABEL)

!         OPEN CORRELATED-K DISTRIBUTIONS FILE.
          IF((IEMSCT.EQ.1 .OR. IEMSCT.EQ.2) .AND. (CODE.EQ.'C' .OR.     &
     &      CODE.EQ.'c' .OR. CODE.EQ.'K' .OR. CODE.EQ.'k'))THEN

!             SET NSPEED
              IF(SPEED.EQ.'M' .OR. SPEED.EQ.'m' .OR. SPEED.EQ.'1')THEN
                  NSPEED=1
              ELSE
                  NSPEED=0
              ENDIF
              IF(LENSTR(CKNAME).EQ.0)THEN
                  CORKNM(5:6)=RESCHR
                  CKNAME=DATDIR(1:LENDAT)//CORKNM
              ENDIF
              IF(CKNAME.NE.CKNMSV .OR. NSPEED.NE.NSPDSV)THEN

!                 NEW CORRELATED-K DISTRIBUTIONS FILE.
                  CALL RDCORK(CKNAME,NSPEED,KNTRVL)
                  CKNMSV=CKNAME
                  NSPDSV=NSPEED
                  KNTRSV=KNTRVL
              ELSE

!                 RESTORE THE OLD "KNTRVL" NUMBER.
                  KNTRVL=KNTRSV
              ENDIF

!             SET UP FOR k-DISTRIBUTION OUTPUT:
              IF(KNTRVL.GT.1 .AND. KPRINT)THEN

!                 K-DISTRIBUTION DEPENDENT DATA REQUESTED:
                  IF(BINOUT)THEN

!                     BINARY FILES:
                      INQUIRE(JTKDIS,OPENED=LOPEN,RECL=LNKDIS)
                      IF(LOPEN)THEN

!                         CHECK RECORD LENGTH:
                          IF(IRECLN(KNTRVL+5).GT.LNKDIS)THEN
                              WRITE(IPR,'(/A,15X,A)')                   &
     &                          'Error in CD1A:  Number of k-intervals' &
     &                          //' is too large for (repeat run)',     &
     &                          ' binary k-distribution dependent'      &
     &                          //' output file.'
                              STOP 'Error in CD1A:  Too many k''s'
                          ENDIF
                      ELSE

!                         OPEN K-DEPENDENT TRANSMITTANCE FILE:
                          LNKDIS=IRECLN(KNTRVL+5)
                          IF(LNFLRT.GT.0)THEN
                              CALL OPNFL(JTKDIS,LNKDIS,FLRT(1:LNFLRT)// &
     &                          'b.t_k','UNKNOWN','UNFORMATTED','CD1A')
                          ELSE
                              CALL OPNFL(JTKDIS,LNKDIS,'t_kdis.bin',    &
     &                          'UNKNOWN','UNFORMATTED','CD1A')
                          ENDIF
                          NTKDIS=0
                      ENDIF
                      INQUIRE(JRKDIS,OPENED=LOPEN,RECL=LNKDIS)
                      IF(LOPEN)THEN

!                         CHECK RECORD LENGTH:
                          IF(IRECLN(KNTRVL+5).GT.LNKDIS)THEN
                              WRITE(IPR,'(/A,15X,A)')                   &
     &                          'Error in CD1A:  Number of k-intervals' &
     &                          //' is too large for (repeat run)',     &
     &                          ' binary k-distribution dependent'      &
     &                          //' output file.'
                              STOP 'Error in CD1A:  Too many k''s'
                          ENDIF
                      ELSE

!                         OPEN K-DEPENDENT RADIANCE FILE:
                          LNKDIS=IRECLN(KNTRVL+5)
                          IF(LNFLRT.GT.0)THEN
                              CALL OPNFL(JRKDIS,LNKDIS,FLRT(1:LNFLRT)// &
     &                          'b.r_k','UNKNOWN','UNFORMATTED','CD1A')
                          ELSE
                              CALL OPNFL(JRKDIS,LNKDIS,'r_kdis.bin',    &
     &                          'UNKNOWN','UNFORMATTED','CD1A')
                          ENDIF
                          NRKDIS=0
                      ENDIF
                  ELSE

!                     ASCII FILES:
                      INQUIRE(ITKDIS,OPENED=LOPEN)
                      IF(.NOT.LOPEN)THEN

!                         OPEN K-DISTRIBUTION TRANSMITTANCE FILE:
                          IF(LNFLRT.GT.0)THEN
                              CALL OPNFL(ITKDIS,0,FLRT(1:LNFLRT)//      &
     &                          '.t_k','UNKNOWN','FORMATTED','CD1A')
                          ELSE
                              CALL OPNFL(ITKDIS,0,'t_kdis.dat',         &
     &                          'UNKNOWN','FORMATTED','CD1A')
                          ENDIF
                      ENDIF
                      INQUIRE(IRKDIS,OPENED=LOPEN)
                      IF(.NOT.LOPEN)THEN

!                         OPEN K-DISTRIBUTION RADIANCE FILE:
                          IF(LNFLRT.GT.0)THEN
                              CALL OPNFL(IRKDIS,0,FLRT(1:LNFLRT)//      &
     &                          '.r_k','UNKNOWN','FORMATTED','CD1A')
                          ELSE
                              CALL OPNFL(IRKDIS,0,'r_kdis.dat',         &
     &                          'UNKNOWN','FORMATTED','CD1A')
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
      ELSE

!         SET BNDWID TO 5 CM-1 FOR LOWTRAN BAND MODEL.
          BNDWID=5.
          HBNDWD=2.5
          RESCHR='00'
      ENDIF
      IF(KNTRVL.EQ.1)KPRINT=.FALSE.

!     SOLAR IRRADIANCE DATA:
      IF(IEMSCT.LE.1)RETURN

!     CHECK FOR INPUT OF SOLAR CONSTANT:
      IF(SOLCON.NE.0. .AND. SFWHM.EQ.0.)THEN

!         SOLAR CONSTANT FROM THE DATABASES IS BEING SCALED.  THE
!         SOLAR FILE MUST BE READ EVEN IF LSUN IS FALSE IN TAPE5.
!         SET SFWHM TO THE DEFAULT VALUE (5 CM-1) IF LSUN IS FALSE.
          IF(.NOT.MODTRN)THEN
              SFWHM=20.
          ELSEIF(RESCHR.EQ.'01')THEN
              SFWHM=5.
          ELSE
              SFWHM=BNDWID
          ENDIF
          LSUN=.TRUE.
      ENDIF

!     CHECK IF SOLAR DATA MUST BE READ IN:
      IF(RESCHR.EQ.'p1')THEN

!         SOLAR DATA IS ALWAYS READ IN WITH THE 0.1 CM-1 BAND MODEL.
          L_USUN=LENSTR(USRSUN)
          IF(L_USUN.EQ.0)THEN

!             'USRSUN' IS BLANK, SO USE THE DEFAULT SOLAR FILE:
              L_USUN=1
              USRSUN(1:1)='7'
              DFTSUN=DATDIR(1:LENDAT)//'SUNp1kurucz1995.bin'
              CALL RSUNP1(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'1')THEN

!             THE KURUCZ 2005 SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUNp1kurucz2005.bin'
              CALL RSUNP1(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'5')THEN

!             THE FONTENLA SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUNp1fontenla.bin'
              CALL RSUNP1(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSE

!             IF 'USRSUN' IS 2, 3, 4 or 6 WITH 0.1 CM-1, RESET TO 7:
              IF(USRSUN(1:L_USUN).EQ.'2' .OR. USRSUN(1:L_USUN).EQ.'3'   &
     &          .OR. USRSUN(1:L_USUN).EQ.'4'                            &
     &          .OR. USRSUN(1:L_USUN).EQ.'6')USRSUN(1:1)='7'

!             'USRSUN' EITHER EQUALS '7' OR A USER FILE WAS SPECIFIED:
              DFTSUN=DATDIR(1:LENDAT)//'SUNp1kurucz1995.bin'
              CALL RSUNP1(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ENDIF

!         WRITE SOLAR IRRADIANCE FILE NAME:
          IF(.NOT.LJMASS)THEN
              IF(USRSUN.EQ.'1' .OR. USRSUN.EQ.'5'                       &
     &                         .OR. USRSUN.EQ.'7')THEN
                  WRITE(IPR,'(/A,30X,A)')                               &
     &              ' SOLAR IRRADIANCE FILE:',DFTSUN(1:LENSTR(DFTSUN))
              ELSE
                  WRITE(IPR,'(/A,30X,A)')                               &
     &              ' SOLAR IRRADIANCE FILE:',USRSUN(1:L_USUN)
              ENDIF
          ENDIF
      ELSEIF(LSUN)THEN

!         DETERMINE THE APPROPRIATE SOLAR FILE AND READ IT:
          L_USUN=LENSTR(USRSUN)
          IF(L_USUN.EQ.0)THEN

!             'USRSUN' IS BLANK, SO USE THE DEFAULT SOLAR FILE:
              L_USUN=1
              USRSUN(1:1)='6'
              !DFTSUN=DATDIR(1:LENDAT)//SUNFIL
              DFTSUN=DATDIR(1:LENDAT)//'SUN01kurucz1997.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'1')THEN

!             THE newkur SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUN01kurucz2005.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'2')THEN

!             THE CHANCE-KURUCZ SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUN01chkur.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'3')THEN

!             THE CEBULA-CHANCE-KURUCZ SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUN01cebchkur.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'4')THEN

!             THE THUILLIER-KURUCZ SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUN01thkur.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'5')THEN

!             THE FONTENLA SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUN01fontenla.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSEIF(USRSUN(1:L_USUN).EQ.'7')THEN

!             THE oldkur SOLAR FILE WAS SELECTED:
              DFTSUN=DATDIR(1:LENDAT)//'SUN01kurucz1995.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ELSE

!             'USRSUN' EITHER EQUALS '6' OR A USER FILE WAS SPECIFIED:
             !DFTSUN=DATDIR(1:LENDAT)//SUNFIL
              DFTSUN=DATDIR(1:LENDAT)//'SUN01kurucz1997.dat'
              CALL RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)
          ENDIF
          LRDSUN=ABS(ABS(SFWHM)-5).GT..01 .OR.                          &
     &      USRSUN(1:L_USUN).NE.'1' .OR. SOLCON.NE.0.

!         WRITE SOLAR IRRADIANCE FILE NAME:
          IF(.NOT.LJMASS)THEN
              IF(USRSUN.EQ.'1' .OR. USRSUN.EQ.'2' .OR. USRSUN.EQ.'3'    &
     &          .OR. USRSUN.EQ.'4' .OR. USRSUN.EQ.'5'                   &
     &          .OR. USRSUN.EQ.'6' .OR. USRSUN.EQ.'7')THEN
                  WRITE(IPR,'(/A,30X,A)')                               &
     &              ' SOLAR IRRADIANCE FILE:',DFTSUN(1:LENSTR(DFTSUN))
              ELSE
                  WRITE(IPR,'(/A,30X,A)')                               &
     &              ' SOLAR IRRADIANCE FILE:',USRSUN(1:L_USUN)
              ENDIF
          ENDIF
      ENDIF

!     RETURN TO DRIVER:
      RETURN
      END
