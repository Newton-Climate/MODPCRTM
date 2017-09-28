      SUBROUTINE GEO(ILOS,MSPATH,ICH1,CKRANG,REARTH,ZMAX,AERRH,LENN,    &
     &  H1ALT,H2ALT,OBSZEN,HRANGE,BETA,HMIN,HMAX,BCKZEN,IERROR,BEND)

!     GEO SERVES AS AN INTERFACE BETWEEN ROUTINE GEODRV AND
!     GEOMETRY SUBROUTINES INCLUDING GEOINP, REDUCE, FDBETA,
!     EXPINT, FINDMN, DPFISH, DPSCHT, DPANDX, DPRARF, RFPATH,
!     FILL AND RFGEOM.  THESE ROUTINES CALCULATE ABSORBER
!     AMOUNTS FOR A REFRACTED PATH THROUGH THE ATMOSPHERE.
      IMPLICIT NONE

!     PARAMETERS:
      DOUBLE PRECISION TLRNCE
      PARAMETER(TLRNCE=0.001D0)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       ILOS     PATH LINE-OF-SIGHT INDEX.
!       CKRANG   MAXIMUM PATH RANGE FOR K-DISTRIBUTION OUTPUT
!                (=0. FOR TOTAL PATH ONLY; <0. FOR ALL RANGES).
!       REARTH   RADIUS OF THE EARTH [KM].
!       ZMAX     MAXIMUM ATMOSPHERIC PROFILE ALTITUDE [KM].
!       AERRH    RELATIVE HUMIDITY USED FOR AEROSOL MODELS.
!       LENN     INPUT PATH LENGTH SWITCH (0=SHORT, 1=LONG).
!       H1ALT    INPUT OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    FINAL (TANGENT FOR LIMB PATH) ALTITUDE [KM].
!       OBSZEN   INPUT OBSERVER ZENITH ANGLE (H1ALT TO H2ALT) [DEG].
!       HRANGE   ACTUAL DISTANCE FROM H1ALT TO H2ALT [KM].
!       BETA     INPUT EARTH CENTER ANGLE BETWEEN H1ALT AND H2ALT [DEG].
!       HMIN     INPUT PATH MINIMUM ALTITUDE [KM].
!       HMAX     PATH MAXIMUM ALTITUDE (ITYPE=4 ONLY) [KM].
!       BCKZEN   INPUT ZENITH ANGLE FOR BACKWARD PATH (H2 TO H1) [DEG].
!       HRANGE   INPUT DISTANCE FROM H1ALT TO H2ALT [KM].
!       MSPATH   LOGICAL, TRUE FOR MULTIPLE SCATTERING VERTICAL PATH.
!       ICH1     NUMERIC LABEL FOR BOUNDARY LAYER AEROSOL MODEL.
      LOGICAL MSPATH
      INTEGER ILOS,ICH1
      DOUBLE PRECISION CKRANG,REARTH,ZMAX,HMAX

!     OUTPUT ARGUMENTS:
!       AERRH    UPDATED AER MODEL RELATIVE HUMIDITY IF INPUT WAS ZERO.
!       LENN     LENGTH SWITCH FOR ACTUAL PATH (0=SHORT, 1=LONG).
!       H1ALT    ACTUAL PATH OBSERVER (SENSOR) ALTITUDE (< ZMAX) [KM].
!       H2ALT    ACTUAL PATH FINAL / TANGENT ALTITUDE (< ZMAX) [KM].
!       OBSZEN   ACTUAL OBSERVER ZENITH ANGLE (H1ALT TO H2ALT) [DEG].
!       HRANGE   ACTUAL DISTANCE FROM H1ALT TO H2ALT [KM].
!       BETA     ACTUAL EARTH CENTER ANGLE BETWEEN H1ALT & H2ALT [DEG].
!       HMIN     ACTUAL PATH MINIMUM ALTITUDE [KM].
!       BCKZEN   ACTUAL ZENITH ANGLE FOR BACKWARD PATH (H2 TO H1) [DEG].
!       BEND     PATH BENDING FROM REFRACTION [DEG].
!       IERROR   INTEGER ERROR FLAG.
      REAL AERRH
      INTEGER LENN,IERROR
      DOUBLE PRECISION H1ALT,H2ALT,OBSZEN,HRANGE,BETA,HMIN,BCKZEN,BEND

!     COMMONS:
      INCLUDE 'YPROP.h'
      INCLUDE 'SEGDAT.h'
      INCLUDE 'IFIL.h'
      INCLUDE 'SOLS.h'

!     /FILLP/
!       ZP       ALTITUDE AT REFRACTED PATH LEVELS [KM].
!       RANGEP   SEGMENT RANGE [KM].
!       PP       PRESSURE AT REFRACTED PATH LEVELS [MBAR].
!       TP       TEMPERATURE AT REFRACTED PATH LEVELS [K].
!       RHP      RELATIVE HUMIDITY AT REFRACTED PATH LEVELS [%].
!       RFNDXP   REFRACTIVITIES AT REFRACTED PATH LEVELS.
!       DENPTH   ATMOSPHERIC LEVEL DENSITIES
!                [UNITS VARY WITH SPECIES BUT MOST OFTEN ATM CM / KM].
      DOUBLE PRECISION ZP,RANGEP,PP,RFNDXP,DENPTH
      REAL TP,RHP
      COMMON/FILLP/ZP(LAYDM1),RANGEP(LAYDM1),PP(LAYDM1),                &
     &  TP(LAYDM1),RHP(LAYDM1),RFNDXP(LAYDM1),DENPTH(0:MEXTXY,1:LAYDM1)

!     /PTHAMT/
!       PPSUM    DENSITY WEIGHTED PATH SEGMENT AVERAGE PRESSURE [MBAR].
!       TPSUM    DENSITY WEIGHTED PATH SEGMENT AVERAGE TEMPERATURE [K].
!       APSAPX   AEROSOL PATH SEGMENT SPECTRAL EXTINCTION OPTICAL DEPTH.
!       APSAPA   AEROSOL PATH SEGMENT SPECTRAL ABSORPTION OPTICAL DEPTH.
!       RHOPSM   REFRACTED PATH SEGMENT COLUMN DENSITY [G/CM3].
!       AMTPTH   PATH SEGMENT COLUMN DENSITIES
!                [UNITS VARY WITH SPECIES BUT OFTEN ATM CM].
      REAL PPSUM,TPSUM,APSAPX,APSAPA,RHOPSM,AMTPTH
      COMMON/PTHAMT/PPSUM(LAYDM1),TPSUM(LAYDM1),APSAPX(MWVSAP,LAYDM1),  &
     &  APSAPA(MWVSAP,LAYDM1),RHOPSM(LAYDM1),AMTPTH(0:MEXTXY,1:LAYDM1)

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

!     /SURFWV/
!       LAMBER   LOGICAL FLAG, .TRUE. FOR LAMBERTIAN SURFACE.
!       TPTEMP   TARGET-PIXEL SURFACE TEMPERATURES [K].
!       TPHDIR   TARGET-PIXEL HEMISPHERE DIRECTIONAL REFLECTANCE AT
!                VIEWING ANGLES.
!       TPBRDF   TARGET-PIXEL BIDIRECTIONAL REFLECTANCE DISTRIBUTION
!                FUNCTION AT VIEWING AND SUN ANGLES.
!       AATEMP   AREA-AVERAGED GROUND SURFACE TEMPERATURES [K].
!       AASALB   AREA-AVERAGED GROUND SURFACE ALBEDO.
!       AADREF   AREA-AVERAGED GROUND SURFACE DIRECTIONAL REFLECTIVITY
!                AT THE SOLAR ZENITH ANGLE.
!       EMU      GROUND DIRECTIONAL EMISSIVITY AT VIEWING ANGLES.
!       BEM      GROUND DIRECTIONAL EMISSIVITY AT QUADRATURE ANGLES.
!       RMU      GROUND BRDF AZIMUTH COMPONENTS AT VIEWING ANGLES
!                AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
!       BDR      GROUND BRDF AZIMUTH COMPONENTS AT QUADRATURE ANGLE
!                AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
      LOGICAL LAMBER
      REAL TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,EMU,BEM,RMU,BDR
      COMMON/SURFWV/LAMBER,TPTEMP(MLOS),TPHDIR(MLOS),TPBRDF(MLOS),      &
     &  AATEMP,AASALB,AADREF,EMU(MXUMU),BEM(MI),                        &
     &  RMU(1:MXUMU,0:MI,0:MAZ),BDR(1:MI,0:MI,0:MAZ)                    &

!     /ANGSRF/
!       CVWSRF   COSINE OF THE VIEW ZENITH ANGLE FROM THE SURFACE.
!       CSNSRF   COSINE OF THE SOLAR (LUNAR) ZENITH AT SURFACE.
!       AZMSRF   RELATIVE AZIMUTH ANGLE (SUN - SENSOR AT SURFACE) [RAD].
!       UMU1     COSINE OF THE PATH NADIR ANGLE.
!                (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       UMU0     COSINE OF THE SOLAR ZENITH ANGLE.
!                (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       PHI1     RELATIVE AZIMUTH ANGLE (SUN - LOS PATH @ SENSOR) [DEG].
!                (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       CMU      COSINE OF THE NADIR ANGLES USED IN DISORT.
      REAL CVWSRF,CSNSRF,AZMSRF,UMU1,UMU0,PHI1,CMU
      COMMON/ANGSRF/CVWSRF(MLOS),CSNSRF(MLOS),AZMSRF(MLOS),UMU1(MLOS),  &
     &  UMU0,PHI1(MLOS),CMU(MI)

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

      DOUBLE PRECISION DHALFR,DPRNG2
      COMMON/SMALL1/DHALFR,DPRNG2
      LOGICAL LSMALL
      COMMON/SMALL2/LSMALL
      LOGICAL LPRINT
      COMMON/CPRINT/LPRINT

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     /RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

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

!     /SAPDAT/
!       BSAPX    AEROSOL EXTINCTION COEFFICIENTS [KM-1].
!       BSAPA    AEROSOL ABSORPTION COEFFICIENTS [KM-1].
      REAL BSAPX,BSAPA
      COMMON/SAPDAT/BSAPX(MWVSAP,LAYDIM),BSAPA(MWVSAP,LAYDIM)

!     /SAPPTH/
!       BPSAPX   AEROSOL EXTINCTION COEFFICIENTS [KM-1].
!       BPSAPA   AEROSOL ABSORPTION COEFFICIENTS [KM-1].
      DOUBLE PRECISION BPSAPX,BPSAPA
      COMMON/SAPPTH/BPSAPX(MWVSAP,LAYDM1),BPSAPA(MWVSAP,LAYDM1)

!     /SAPDEP/
!       XSAP     LINE-OF-SIGHT PATH AEROSOL EXTINCTION OPTICAL DEPTHS.
!       XSAPM    VERTICAL MS PATH AEROSOL EXTINCTION OPTICAL DEPTHS.
!       ASAP     LINE-OF-SIGHT PATH AEROSOL ABSORPTION OPTICAL DEPTHS.
!       ASAPM    VERTICAL MS PATH AEROSOL ABSORPTION OPTICAL DEPTHS.
!       XSAPS    SOLAR PATH AEROSOL EXTINCTION OPTICAL DEPTHS FROM LOS.
!       XSAPSM   SOLAR PATH AEROSOL EXTINCTION OPTICAL DEPTHS FOR MS.
!       ASAPS    SOLAR PATH AEROSOL ABSORPTION OPTICAL DEPTHS FROM LOS.
!       ASAPSM   SOLAR PATH AEROSOL ABSORPTION OPTICAL DEPTHS FOR MS.
      REAL XSAP,XSAPM,ASAP,ASAPM,XSAPS,XSAPSM,ASAPS,ASAPSM
      COMMON/SAPDEP/                                                    &
     &  XSAP(1:MWVSAP,0:LAYTWO,0:MLOS),XSAPM(1:MWVSAP,0:LAYDIM),        &
     &  ASAP(1:MWVSAP,0:LAYTWO,0:MLOS),ASAPM(1:MWVSAP,0:LAYDIM),        &
     &  XSAPS(MWVSAP,LAYTWO,MLOS),XSAPSM(MWVSAP,LAYDIM),                &
     &  ASAPS(MWVSAP,LAYTWO,MLOS),ASAPSM(MWVSAP,LAYDIM)

!     /MOL_FM/
!       UMX_FM   ATM LEVEL MOLE FRACTION RATIO, RELATIVE TO CO2.
!       AUX_FM   ATM LEVEL AUX MOLECULAR GASES RELATIVE MOLE FRACTION.
      REAL UMX_FM,AUX_FM
      COMMON/MOL_FM/UMX_FM(3:NMOLXT,1:LAYDIM),AUX_FM(MMOLY,LAYDIM)

!     /MOL_FP/
!       UMX_FP   PATH LEVEL MOLE FRACTION RATIO, RELATIVE TO CO2.
!       AUX_FP   PATH LEVEL AUX MOLECULAR GASES RELATIVE MOLE FRACTION.
      REAL UMX_FP,AUX_FP
      COMMON/MOL_FP/UMX_FP(3:NMOLXT,1:LAYDM1),AUX_FP(MMOLY,LAYDM1)

!     /MOL_F/
!       UMX_F    LOS LEVEL MOLE FRACTION RATIO, RELATIVE TO CO2.
!       AUX_F    LOS LEVEL AUX MOLECULAR GASES RELATIVE MOLE FRACTION.
      REAL UMX_F,AUX_F
      COMMON/MOL_F/UMX_F(2:NMOLXT,0:LAYTWO),AUX_F(1:MMOLY,0:LAYTWO)

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER AMOD3D*1
      COMMON/CJM5/AMOD3D

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD,ATMCON,XMLATM

!     FUNCTIONS:
      INTEGER PSLCT
      REAL EXPINT

!     LOCAL VARIABLES:
!       NSEGX    NUMBER OF PATH SEGMENTS FOR CURRENT LINE-OF-SIGHT.
!       ISEG     PATH SEGMENT INDEX.
!       ISEGM1   ISEG MINUS 1.
!       IMOL     MOLECULAR INDEX.
!       KSPEC    SPECIES LOOP INDEX.
!       IWVSAP   SPECTRAL GRID INDEX FOR SPECTRAL AEROSOL PROFILES.
!       TLAYC    DENSITY WEIGHTED LAYER TEMPERATURE [C]
      INTEGER I,NSEGX,ISEG,ISEGM1,KSPEC,L,ILO,IHI,IH1ALT,IH2ALT,        &
     &  ISLCT,LENNSV,LDEL,IWVSAP,IMOL
      LOGICAL LSAVE,LNOGEO
      REAL FAC,TLAYC,PRCNT
      DOUBLE PRECISION H2SAV,RANGSV,BZENX

!     INITIALIZE CONSTANTS AND CLEAR CUMULATIVE VARIABLES:
      H2SAV=H2ALT
      LSMALL=.FALSE.
      LPRINT=.TRUE.
      LNOGEO=.FALSE.
      IF(.NOT.LJMASS .AND. GNDALT.GT.ZM(1) .AND. .NOT.LSSGEO)           &
     &  WRITE(IPR,'(/A,2F10.6,I5)')                                     &
     &  ' GNDALT is above first profile altitude:',GNDALT,ZM(1),MODEL
      IERROR=0

!     INITIALIZE CUMULATIVE VARIABLES:
      DO I=1,LAYDM1
          LAYSEG(I)=0
          RANGEP(I)=0.D0
          PPSUM(I)=0.
          TPSUM(I)=0.
          RHOPSM(I)=0.
          DO KSPEC=0,MEXTX+NMOLY
              AMTPTH(KSPEC,I)=0.
          ENDDO
          IF(LSAP)THEN
              DO IWVSAP=1,NWVSAP
                  APSAPX(IWVSAP,I)=0.
                  APSAPA(IWVSAP,I)=0.
              ENDDO
          ENDIF
      ENDDO
      IF(ITYPE.LT.2)THEN

!         HORIZONTAL PATH:  INTERPOLATE PROFILE TO H1ALT.
          ZP(1)=H1ALT
          IF(ML.EQ.1)THEN
              PP(1)=DBLE(PM(1))
              PATH_T(0,1)=TM(1)
              PATH_P(0,1)=PM(1)/PZERO
              TP(1)=TM(1)
              LAYSEG(1)=1
              DO KSPEC=0,MEXTX+NMOLY
                  DENPTH(KSPEC,1)=DBLE(DENSTY(KSPEC,1))
              ENDDO
              IF(LSAP)THEN
                  DO IWVSAP=1,NWVSAP
                      BPSAPX(IWVSAP,1)=DBLE(BSAPX(IWVSAP,1))
                      BPSAPA(IWVSAP,1)=DBLE(BSAPA(IWVSAP,1))
                  ENDDO
              ENDIF
          ELSEIF(H1ALT.LT.ZM(1))THEN
              WRITE(IPR,'(/3(A,F10.6))')' Error in GEO:  Horizontal'    &
     &          //' path H1ALT [=',H1ALT,'] < ZM(1) [=',ZM(1),'].'
              IERROR=1
              RETURN
          ELSE
              ILO=1
              DO IHI=2,ML
                  IF(H1ALT.LT.ZM(IHI))GOTO 10
                  ILO=IHI
              ENDDO
   10         FAC=SNGL((H1ALT-ZM(ILO))/(ZM(IHI)-ZM(ILO)))
              PP(1)=DBLE(EXPINT(PM(ILO),PM(IHI),FAC))
              TP(1)=TM(ILO)+(TM(IHI)-TM(ILO))*FAC
              PATH_T(0,1)=TP(1)
              PATH_P(0,1)=SNGL(PP(1))/PZERO
              LAYSEG(1)=ILO
              IF(FAC.GT..5)LAYSEG(1)=IHI
              DO KSPEC=0,MEXTX+NMOLY
                  DENPTH(KSPEC,1)=DBLE(EXPINT(DENSTY(KSPEC,ILO),        &
     &                                      DENSTY(KSPEC,IHI),FAC))
              ENDDO
              IF(LSAP)THEN
                  IF(IHI.LE.LEVSAP)THEN
                      DO IWVSAP=1,NWVSAP
                          BPSAPX(IWVSAP,1)                              &
     &                      =DBLE(EXPINT(BSAPX(IWVSAP,ILO),             &
     &                                   BSAPX(IWVSAP,IHI),FAC))
                          BPSAPA(IWVSAP,1)                              &
     &                      =DBLE(EXPINT(BSAPA(IWVSAP,ILO),             &
     &                                   BSAPA(IWVSAP,IHI),FAC))
                      ENDDO
                  ELSE
                      DO IWVSAP=1,NWVSAP
                          BPSAPX(IWVSAP,1)=0.D0
                          BPSAPA(IWVSAP,1)=0.D0
                      ENDDO
                  ENDIF
              ENDIF

!             USE LINEAR INTERPOLATION FOR CLOUDS
              DENPTH(16,1)=DBLE(DENSTY(16,ILO)                          &
     &          +FAC*(DENSTY(16,IHI)-DENSTY(16,ILO)))
              DENPTH(66,1)=DBLE(DENSTY(66,ILO)                          &
     &          +FAC*(DENSTY(66,IHI)-DENSTY(66,ILO)))
              DENPTH(67,1)=DBLE(DENSTY(67,ILO)                          &
     &          +FAC*(DENSTY(67,IHI)-DENSTY(67,ILO)))
          ENDIF

!         CALCULATE ABSORBER AMOUNTS FOR A HORIZONTAL PATH
          IF(.NOT.LJMASS)WRITE(IPR,'(/2(A,F11.3),A,I4)')                &
     &      ' HORIZONTAL PATH AT ALTITUDE =',H1ALT,                     &
     &      ' KM WITH HRANGE =',HRANGE,' KM, MODEL =',MODEL
          NSEG(1)=1
          IF(MODEL.EQ.0)THEN
              PP(1)=DBLE(PM(1))
              TP(1)=TM(1)
              PP(1)=DBLE(PM(1))
          ENDIF
          PATM(1,1)=SNGL(PP(1))/PZERO
          TSEG(1,1)=TP(1)
          DRNG(1,1)=SNGL(HRANGE)
          PTHRNG(1,1)=HRANGE
          DO KSPEC=0,MEXTX+NMOLY
              WTOTAL(KSPEC)=SNGL(DENPTH(KSPEC,1)*HRANGE)
              WPTH(KSPEC,1,1)=WTOTAL(KSPEC)
          ENDDO
          IF(LSAP)THEN
              DO IWVSAP=1,NWVSAP
                  XSAP(IWVSAP,0,1)=SNGL(BPSAPX(IWVSAP,1)*HRANGE)
                  XSAP(IWVSAP,1,1)=XSAP(IWVSAP,0,1)
                  ASAP(IWVSAP,0,1)=SNGL(BPSAPA(IWVSAP,1)*HRANGE)
                  ASAP(IWVSAP,1,1)=ASAP(IWVSAP,0,1)
              ENDDO
          ENDIF

!         H2O CONTINUUM BAND MODEL COLUMN DENSITY DATA:
          WTOTAL(9)=WTOTAL(5)*(296-TP(1))/36

!         OZONE CONTINUUM BAND MODEL COLUMN DENSITY DATA:
          TLAYC=TP(1)-TZERO
          WTOTAL(59)=WTOTAL(8)*.269*TLAYC
          WTOTAL(60)=WTOTAL(59)*TLAYC
          WPTH(9,1,1)=WTOTAL(9)
          WPTH(59,1,1)=WTOTAL(59)
          WPTH(60,1,1)=WTOTAL(60)
      ELSE

!         SLANT PATH SELECTED.  INTERPRET SLANT PATH PARAMETERS.
!         ISLCT IS POSITIVE ONLY FOR AN OPTICAL PATH WITH ITYPE=2.
          IF(ITYPE.NE.2 .OR. MSPATH)THEN
              ISLCT=0
          ELSE
              ISLCT=PSLCT(OBSZEN,HRANGE,BETA)

!             SPECIAL TREATMENT EXCEPT FOR CASE 2A (ISLCT=1):
              IF(ISLCT.GT.1)THEN

!                 IF HRANGE IS SMALL, CONVERT TO CASE 2C (ISLCT=3):
                  CALL SMPREP(REARTH,H1ALT,H2ALT,OBSZEN,                &
     &              HRANGE,BETA,ISLCT)
                  IF(HRANGE.GT.0.D0 .AND. HRANGE.LE.RSMALL)THEN
                      LSMALL=.TRUE.
                      RANGSV=HRANGE
                      ISLCT=3
                  ELSEIF(ISLCT.EQ.2)THEN

!                     CASE 2B:  (H1ALT,OBSZEN,HRANGE)
!                     IF PATH TYPE IS CASE 2B, CHECK THAT THE HRANGE
!                     USED IN THE CALCULATION EQUALS THE INPUT HRANGE.
!                     DETERMINE H2ALT AND LENN USING ROUTINE NEWH2.
                      LENNSV=LENN
                      RANGSV=HRANGE
                      CALL NEWH2(REARTH,ZMAX,H1ALT,H2ALT,               &
     &                  OBSZEN,HRANGE,BETA,LENN,HMIN,BCKZEN)
                      IF(LENN.EQ.0)HMIN=MIN(H2ALT,H1ALT)
                      LPRINT=.FALSE.
                      BZENX=BCKZEN
                      CALL RFPATH(REARTH,H1ALT,H2ALT,OBSZEN,LENN,       &
     &                  HMIN,.FALSE.,BZENX,BETA,HRANGE,BEND)
                      LPRINT=.TRUE.
                      PRCNT=SNGL(100*ABS(HRANGE-RANGSV)/RANGSV)
                      IF(.NOT.LJMASS)THEN
                          WRITE(IPR,'(/A)')' SOME INTERNAL DETAILS:'
                          WRITE(IPR,'(/(A))')' LOS IS INTERNAL CASE 2B' &
     &                      //' (H1ALT, OBSZEN, HRANGE).',' USING'      &
     &                      //' H2ALT OBTAINED FROM SUBROUTINE NEWH2:'
                          WRITE(IPR,'(/A,9(F12.5,/A),I12,/))')          &
     &                      ' H1ALT            =',H1ALT,                &
     &                      ' H2ALT            =',H2ALT,                &
     &                      ' OBSZEN           =',OBSZEN,               &
     &                      ' BCKZEN           =',BCKZEN,               &
     &                      ' BETA             =',BETA,                 &
     &                      ' HMIN (MIN ALT)   =',HMIN,                 &
     &                      ' HRANGE (OUTPUT)  =',HRANGE,               &
     &                      ' HRANGE (INPUT)   =',RANGSV,               &
     &                      ' PERCENT DIFF     =',PRCNT,                &
     &                      ' LENN             =',LENN
                      ENDIF
                      IF(OBSZEN.GT.90.D0 .AND.                          &
     &                  ABS(H2ALT-GNDALT).LT.TLRNCE)THEN
                          LNOGEO=.TRUE.
                          WRITE(IPR,'(/2A,/)')'***** WARNING ***** ',   &
     &                      ' PATH HITS THE EARTH.'
                      ELSEIF(OBSZEN.LE.90.D0 .AND.                      &
     &                  ABS(H2ALT-ZM(ML)).LT.TLRNCE)THEN
                          LNOGEO=.TRUE.
                          WRITE(IPR,'(/2A,/)')'***** WARNING ***** ',   &
     &                      ' PATH HITS THE UPPERMOST LAYER BOUNDARY.'
                      ELSEIF(PRCNT.LE.1. .OR. H2ALT.GE.ZM(ML))THEN
                          IF(.NOT.LJMASS)WRITE(IPR,'(/(2A))')           &
     &                      ' PERCENT DIFFERENCE IS THAN 1 OR THE',     &
     &                      ' PATH TERMINATES AT THE TOP OF THE',       &
     &                      ' ATMOSPHERE; THESE PATH PARAMETERS',       &
     &                      ' WILL BE USED WITHOUT CALLING GEOINP.'
                          LNOGEO=.TRUE.
                      ELSE
                          IF(.NOT.LJMASS)WRITE(IPR,'(3(/A))')           &
     &                      ' SINCE THE PERCENT DIFFERENCE EXCEEDS 1,', &
     &                      ' AN EQUIVALENT INTERNAL CASE 2C IS USED.'
                          LNOGEO=.FALSE.
                          HRANGE=RANGSV
                          BETA=0.D0
                          HMIN=0.D0
                          BCKZEN=0.D0
                          LENN=LENNSV
                          OBSZEN=0.D0
                          ISLCT=PSLCT(OBSZEN,HRANGE,BETA)
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
          LSAVE=LSMALL
          LSMALL=.FALSE.
          IF(.NOT.LNOGEO .AND. ITYPE.NE.4)                              &
     &      CALL GEOINP(REARTH,ZMAX,H1ALT,H2ALT,OBSZEN,HRANGE,          &
     &      BETA,ITYPE,LENN,HMIN,BCKZEN,IERROR,ISLCT)
          LSMALL=LSAVE

!         CHECK FOR NON-ZERO IERROR:
          IF(IERROR.NE.0)THEN
              IF(.NOT.LJMASS .AND. .NOT.LSSGEO)                         &
     &          WRITE(IPR,'(/2A)')'GEO:  IERROR NE 0:',                 &
     &          ' END THIS CALCULATION AND SKIP TO THE NEXT CASE'
              RETURN
          ELSEIF(LSMALL)THEN
              HRANGE=RANGSV
              CALL SMGEO(H1ALT,H2ALT,HRANGE,REARTH,OBSZEN,              &
     &          BETA,BCKZEN,DHALFR,DPRNG2,BEND,LENN,HMIN)
          ENDIF

!         INITIALIZE COLUMN DENSITY SUMS:
          DO KSPEC=1,MEXTX+NMOLY
              WTOTAL(KSPEC)=0.
          ENDDO
          IF(LSAP)THEN
              IF(MSPATH)THEN
                  DO IWVSAP=1,NWVSAP
                      XSAPM(IWVSAP,0)=0.
                      ASAPM(IWVSAP,0)=0.
                  ENDDO
              ELSE
                  DO IWVSAP=1,NWVSAP
                      XSAP(IWVSAP,0,ILOS)=0.
                      ASAP(IWVSAP,0,ILOS)=0.
                  ENDDO
              ENDIF
          ENDIF

!         USER-SPECIFIED OR CALCULATED REFRACTIVE PATH?
          IF(ITYPE.EQ.4)THEN

!             COMPUTE COLUMN DENSITIES FOR USER-DEFINED REFRACTED PATHS.
              NSEGX=NSEG(1)
              CALL GEOPTH(H1ALT,H2ALT,OBSZEN,HRANGE,BETA,REARTH,        &
     &          HMIN,HMAX,CKRANG,BCKZEN,BEND)
          ELSE

!             CALCULATE ATMOSPHERIC PATH.  FOR NEAR HALF TANGENT PATHS
!             [HRANGE EXCEEDING RSMALL, OBSZEN NEAR 90 DEG, AND H1ALT
!             NEAR HMIN], DO NOT LET LENN EQUAL 1 BECAUSE HMIN MAY BE
!             IDENTICAL TO H1ALT PRODUCING PROBLEMS IN ROUTINE FILL.
!             THE CHECK ON VARIABLE RSMALL SIMPLY POINTS OUT THE FACT
!             THAT SMALL PATHS ARE ALREADY DEALT WITH BY OTHER METHODS.
              IF(ABS(OBSZEN-90.D0).LT..001 .AND. HRANGE.GT.RSMALL       &
     &          .AND. ABS(HMIN-H1ALT).LT.TLRNCE)LENN=0
              CALL RFPATH(REARTH,H1ALT,H2ALT,OBSZEN,LENN,               &
     &          HMIN,.TRUE.,BCKZEN,BETA,HRANGE,BEND)
              IF(HMIN.LT.GNDALT-TLRNCE .AND.                            &
     &          (ISLCT.GT.0 .OR. ITYPE.EQ.3))THEN

!                 ISLCT>0 ONLY DEALS WITH ITYPE=2.
!                 THEREFORE, THE CHECK FOR ITYPE=3.
                  IF(.NOT.LJMASS)                                       &
     &              WRITE(IPR,'(2(/A,F12.6),/(A))')'GNDALT =',GNDALT,   &
     &              'HMIN =',HMIN,'HMIN is less than GNDALT.',          &
     &              'This run is being aborted; try next run.',' '
                  IERROR=1
                  RETURN
              ENDIF

!             LOAD AMTPTH LAYER AMOUNTS INTO WPTH FROM H1ALT TO H2ALT:
              DO I=1,JMAX
                  IF(H1ALT.EQ.ZP(I))IH1ALT=I
                  IF(H2ALT.EQ.ZP(I))IH2ALT=I
              ENDDO
              NSEGX=(JMAX-1)+LENN*(MIN0(IH1ALT,IH2ALT)-1)
              IF(NSEGX.GT.LAYTWO)THEN
                  WRITE(IPR,'(/A)')'Error in GEO:  Number of path'//    &
     &              ' segments, NSEGX is greater than parameter LAYTWO.'
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP 'NSEGX is greater than parameter LAYTWO'
              ENDIF

!             DETERMINE LAYSEG(ISEG), THE LAYER NUMBER L=LAYSEG(ISEG) IN
!             AMTPTH(K,L), STARTING FROM HMIN, WHICH CORRESPONDS TO THE
!             LAYER ISEG IN WPTH(K,ISEG+MS,ILOS), STARTING FROM H1ALT.
              IF(LENN.EQ.0 .AND. H1ALT.LE.H2ALT)THEN

!                 INITIAL DIRECTION OF PATH IS UP.
                  L=0
                  LDEL=1
              ELSE

!                 INITIAL DIRECTION OF PATH IS DOWN.
                  L=IH1ALT
                  LDEL=-1
              ENDIF
              JTURN=0
              DO ISEG=1,NSEGX+1

!                 TEST FOR REVERSING DIRECTION OF PATH FROM DOWN TO UP
                  IF(L.EQ.1 .AND. LDEL.EQ.-1)THEN
                      JTURN=ISEG
                      L=0
                      LDEL=1
                  ENDIF
                  L=L+LDEL
                  LAYSEG(ISEG)=L
              ENDDO

!             CHECK FOR SPECIAL K-DATA OUTPUT INDEX CONDITIONS:
              NSEG(ILOS)=NSEGX
              IF(.NOT.LSSGEO .AND. .NOT.MSPATH)THEN
                  IF(CKRANG.LT.0.D0 .OR. CKRANG.GE.HRANGE)THEN

!                     WRITE K DATA AT ALL LEVELS.
                      IKOUT(ILOS)=NSEGX
                  ELSEIF(CKRANG.GT.0.D0)THEN

!                     SET IKOUT=NSEGX+1 TO FLAG A POSITIVE RANGE:
                      IKOUT(ILOS)=NSEGX+1
                  ELSE

!                     CKRANG IS ZERO:  WRITE K DATA FOR TOTAL PATH ONLY.
                      IKOUT(ILOS)=0
                  ENDIF

!                 DETERMINE OPTICAL PATH IKHMIN:
                  IF(JTURN.GT.0)THEN

!                     TANGENT POINT OCCURRED:
                      IKHMIN(ILOS)=JTURN-1
                  ELSEIF(LDEL.EQ.-1)THEN

!                     HMIN IS AT H2ALT:
                      IKHMIN(ILOS)=NSEGX
                  ELSE

!                     HMIN IS AT H1ALT:
                      IKHMIN(ILOS)=0
                  ENDIF
              ENDIF

!             LOAD TSEG, PATM, DRNG, PTHRNG AND WPTH:
              ISEGM1=0
              DO ISEG=1,NSEGX
                  L=LAYSEG(ISEG)
                  IF(L.GE.ML)L=ML

!                 H2O CONTINUUM BAND MODEL COLUMN DENSITY DATA:
                  AMTPTH(9,L)=AMTPTH(5,L)*(296-TPSUM(L))/36

!                 OZONE CONTINUUM BAND MODEL COLUMN DENSITY DATA:
                  TLAYC=TPSUM(L)-TZERO
                  AMTPTH(59,L)=.269*AMTPTH(8,L)*TLAYC
                  AMTPTH(60,L)=AMTPTH(59,L)*TLAYC
                  IF(MSPATH)THEN
                      TSEGMS(ISEG)=TPSUM(L)
                      PATMMS(ISEG)=PPSUM(L)/PZERO
                      DO KSPEC=0,MEXTX+NMOLY
                          WPTHMS(KSPEC,ISEG)=AMTPTH(KSPEC,L)
                          WTOTAL(KSPEC)=WTOTAL(KSPEC)+AMTPTH(KSPEC,L)
                      ENDDO
                      IF(LSAP)THEN
                          DO IWVSAP=1,NWVSAP
                              XSAPM(IWVSAP,ISEG)=APSAPX(IWVSAP,L)
                              XSAPM(IWVSAP,0)                           &
     &                          =XSAPM(IWVSAP,0)+APSAPX(IWVSAP,L)
                              ASAPM(IWVSAP,ISEG)=APSAPA(IWVSAP,L)
                              ASAPM(IWVSAP,0)                           &
     &                          =ASAPM(IWVSAP,0)+APSAPA(IWVSAP,L)
                          ENDDO
                      ENDIF
                  ELSE
                      TSEG(ISEG,ILOS)=TPSUM(L)
                      PATM(ISEG,ILOS)=PPSUM(L)/PZERO
                      DRNG(ISEG,ILOS)=SNGL(RANGEP(L))
                      IF(.NOT.LSSGEO)THEN

!                         STORE PATH RANGE; COMPUTE K-DATA OUTPUT INDEX:
                          PTHRNG(ISEG,ILOS)                             &
     &                      =PTHRNG(ISEGM1,ILOS)+RANGEP(L)
                          IF(PTHRNG(ISEG,ILOS).GE.CKRANG .AND.          &
     &                      IKOUT(ILOS).GT.NSEGX)IKOUT(ILOS)=ISEG
                      ENDIF
                      DO KSPEC=0,MEXTX+NMOLY
                          WPTH(KSPEC,ISEG,ILOS)=AMTPTH(KSPEC,L)
                          WTOTAL(KSPEC)=WTOTAL(KSPEC)+AMTPTH(KSPEC,L)
                      ENDDO
                      IF(LSAP)THEN
                          DO IWVSAP=1,NWVSAP
                              XSAP(IWVSAP,ISEG,ILOS)=APSAPX(IWVSAP,L)
                              XSAP(IWVSAP,0,ILOS)                       &
     &                          =XSAP(IWVSAP,0,ILOS)+APSAPX(IWVSAP,L)
                              ASAP(IWVSAP,ISEG,ILOS)=APSAPA(IWVSAP,L)
                              ASAP(IWVSAP,0,ILOS)                       &
     &                          =ASAP(IWVSAP,0,ILOS)+APSAPA(IWVSAP,L)
                          ENDDO
                      ENDIF
                  ENDIF
                  ISEGM1=ISEG
              ENDDO
          ENDIF

!NEXT
!         PRINT LAYER PATH SEGMENT ABSORBER AMOUNTS
          IF(.NOT.LSSGEO)THEN

!             IF INPUT SURFACE TEMPERATURE IS 0K & SLANT PATH INTERSECTS
!             EARTH, SET SURFACE TEMPERATURE TO GROUND AIR TEMPERATURE.
              IF(ILOS.LE.MLOS .AND. H2ALT.LE.ZM(1))THEN
                  IF(TPTEMP(ILOS).LE.0.)TPTEMP(ILOS)=TM(1)
              ENDIF

!             UNLESS ITYPE IS 4, IPTHHT AND PATH_T MUST BE DETERMINED:
              IF(ITYPE.NE.4)THEN

                  IF(MSPATH)THEN

!                     MULTIPLE SCATTERING PATH:
                      DO ISEG=1,NSEGX
                          PTH_MS(ISEG)=ZP(LAYSEG(ISEG)+1)
                      ENDDO
                  ELSE

!                     LINEARLY INTERPOLATE TEMPERATURE AT H1ALT:
                      ILO=1
                      DO IHI=2,ML-1
                          IF(ZM(IHI).GE.H1ALT)GOTO 20
                          ILO=IHI
                      ENDDO
                      IHI=ML
   20                 CONTINUE
                      FAC=SNGL((H1ALT-ZM(ILO))/(ZM(IHI)-ZM(ILO)))
                      PATH_T(0,ILOS)=TM(ILO)+FAC*(TM(IHI)-TM(ILO))
                      PATH_P(0,ILOS)                                    &
     &                  =(PM(ILO)/PZERO)*(PM(IHI)/PM(ILO))**FAC
                      PTHALT(0,ILOS)=H1ALT
                      IPTHHT(0)=NINT(1000*H1ALT)

!                     LDEL=0 => GOING DOWN;    LDEL=1 => GOING UP.
                      IF(LENN.EQ.1 .OR. H1ALT.GT.H2ALT)THEN
                          LDEL=0
                      ELSE
                          LDEL=1
                      ENDIF
                      DO ISEG=1,NSEGX
                          IF(ISEG.EQ.JTURN)LDEL=1
                          L=LAYSEG(ISEG)+LDEL
                          PTHALT(ISEG,ILOS)=ZP(L)
                          IPTHHT(ISEG)=NINT(1000*ZP(L))
                          PATH_T(ISEG,ILOS)=TP(L)
                          PATH_P(ISEG,ILOS)=SNGL(PP(L))/PZERO
                      ENDDO

!                     CHECK FOR 'C' (AMOD3D ONLY RESET TO 'M' IN CARD5):
                      IF(AMOD3D.EQ.'C')THEN
                          UMX_F(2,0)=1.
                          DO IMOL=3,NMOLXT
                              UMX_F(IMOL,0)=EXPINT                      &
     &                          (UMX_FM(IMOL,ILO),UMX_FM(IMOL,IHI),FAC)
                          ENDDO
                          DO IMOL=1,NMOLY
                              AUX_F(IMOL,0)=EXPINT                      &
     &                          (AUX_FM(IMOL,ILO),AUX_FM(IMOL,IHI),FAC)
                          ENDDO
                          DO ISEG=1,NSEGX
                              IF(ISEG.EQ.JTURN)LDEL=1
                              L=LAYSEG(ISEG)+LDEL
                              UMX_F(2,ISEG)=1.
                              DO IMOL=3,NMOLXT
                                  UMX_F(IMOL,ISEG)=UMX_FP(IMOL,L)
                              ENDDO
                              DO IMOL=1,NMOLY
                                  AUX_F(IMOL,ISEG)=AUX_FP(IMOL,L)
                              ENDDO
                          ENDDO
                      ENDIF
                  ENDIF
              ENDIF

!             OUTPUTS AVAILABLE IN COMMON BLOCKS
              IF(.NOT.LJMASS)THEN
                  IF(MSPATH)THEN
                      IF(NPR.LT.1)CALL GEOTBL(NSEGX,                    &
     &                  PTH_MS,TSEGMS,WPTHMS,XSAPM,LSAP,MODTRN)
                      WRITE(IPR,'(//A)')' SUMMARY OF MULTIPLE'//        &
     &                  ' SCATTERING VERTICAL PATH GEOMETRY CALCULATION'
                  ELSE
                      IF(NPR.LT.1)                                      &
     &                  CALL GEOTBL(NSEGX,PTHALT(0,ILOS),TSEG(1,ILOS),  &
     &                  WPTH(0,1,ILOS),XSAP(1,0,ILOS),LSAP,MODTRN)
                      WRITE(IPR,'(//A,I3,A)')' SUMMARY OF LINE-OF-SIGHT'&
     &                  //' No.',ILOS,' GEOMETRY CALCULATION'
                  ENDIF
                  WRITE(IPR,'(9(/10X,A,F11.5,1X,A),/10X,A,I11)')        &
     &              'H1ALT   =', H1ALT, 'KM','H2ALT   =', H2ALT, 'KM',  &
     &              'OBSZEN  =',OBSZEN,'DEG','HRANGE  =',HRANGE, 'KM',  &
     &              'BETA    =',  BETA,'DEG','BCKZEN  =',BCKZEN,'DEG',  &
     &              'HMIN    =',  HMIN, 'KM','BENDING =',  BEND,'DEG',  &
     &              'CKRANG  =',CKRANG, 'KM','LENN    =',  LENN
              ENDIF

!             SAVE ZENITH ANGLE FROM H2ALT TO H1ALT AND INITIALIZE THE
!             COSINE SOLAR ZENITH AND RELATIVE AZIMUTH ANGLES AT H2ALT.
              IF(ILOS.LE.MLOS)THEN
                  CVWSRF(ILOS)=COS(SNGL(BCKZEN)/DEG)
                  AZMSRF(ILOS)=PI
                  IF(IMULT.EQ.1)THEN
                      UMU1(ILOS)=-COS(SNGL(OBSZEN)/DEG)
                  ELSE
                      UMU1(ILOS)=CVWSRF(ILOS)
                  ENDIF

!                 AVOID UMU1(ILOS) NEAR ZERO FOR PLANE-PARALLEL
!                 MULTIPLE SCATTERING:

!*************VINCENT ROSS CHANGED FOR BETTER HORIZONTAL CALC***************
#ifdef DISORT_HORIZON
                  IF(ABS(UMU1(ILOS)).LT..002 .OR.                       &
     &              (ITYPE.EQ.3 .AND. H2SAV.GT.0.))THEN
                      IF(UMU1(ILOS).GE.0.)THEN
                          UMU1(ILOS)=.002
                      ELSE
                          UMU1(ILOS)=-.002
                      ENDIF
                  ENDIF
#else
                  IF(ABS(UMU1(ILOS)).LT..05 .OR.                        &
     &              (ITYPE.EQ.3 .AND. H2SAV.GT.0.))THEN
                      IF(UMU1(ILOS).GE.0.)THEN
                          UMU1(ILOS)=.05
                      ELSE
                          UMU1(ILOS)=-.05
                      ENDIF
                  ENDIF
#endif
!****************************END VINCENT ROSS********************************

                  IF(ILOS.EQ.1)THEN
                      CSNSRF(1)=0.
                      UMU0=0.
                  ENDIF
              ENDIF
          ENDIF
      ENDIF

!     CALCULATE THE AEROSOL WEIGHTED MEAN RELATIVE HUMIDITY (RH):
      IF(WTOTAL(7).GT.0. .AND. ICH1.LE.7)THEN
          WTOTAL(15)=100-EXP(WTOTAL(15)/WTOTAL(7))
      ELSEIF(WTOTAL(12).GT.0. .AND. ICH1.GT.7)THEN
          WTOTAL(15)=100-EXP(WTOTAL(15)/WTOTAL(12))
      ELSE
          WTOTAL(15)=0.
      ENDIF

!     USE THE MS PATH AVERAGE RH FOR SELECTING AEROSOLS:
      IF(MSPATH .AND. AERRH.LE.0.)AERRH=WTOTAL(15)

!     CONVERT MOLECULAR BAND WPTH TO CUMULATIVE PATH AMOUNTS IF LOWTRAN:
      IF(.NOT.MODTRN)THEN
          ISEGM1=1
          DO ISEG=2,NSEGX
              DO KSPEC=17,57
                  WPTH(KSPEC,ISEG,ILOS)                                 &
     &              =WPTH(KSPEC,ISEGM1,ILOS)+WPTH(KSPEC,ISEG,ILOS)
              ENDDO
              ISEGM1=ISEG
          ENDDO
      ENDIF

!     WRITE TOTAL PATH AMOUNTS:
      CALL GEOTOT(H1ALT,H2ALT,ILOS,LSSGEO,MSPATH,LSAP,IERROR)
      RETURN
      END

      SUBROUTINE GEOTBL(NSEGX,PTHALT,TSEG,WPTH,XSAP,LSAP,MODTRN)

!     GEOTBL WRITES OUT PATH COLUMN DENSITY TABLES
      IMPLICIT NONE

!     PARAMETERS:
!       RAY550   550 NM RAYLEIGH (C6DTA) SCATTERING COEFFICIENT [KM-1].
      REAL RAY550
      PARAMETER(RAY550=.0121181)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       NSEGX    NUMBER OF PATH SEGMENTS FOR CURRENT LINE-OF-SIGHT.
!       PTHALT   ALTITUDES AT PATH BOUNDARIES [KM].
!       TSEG     PATH SEGMENT DENSITY WEIGHTED TEMPERATURE [K].
!       WPTH     PATH INCREMENTAL COLUMN DENSITIES.
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
!       MODTRN   MODTRAN BAND MODEL FLAG.
      INTEGER NSEGX
      REAL TSEG(*),WPTH(0:MEXTXY,1:*),XSAP(1:MWVSAP,0:*)
      DOUBLE PRECISION PTHALT(0:*)
      LOGICAL LSAP,MODTRN

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /AER/ THERE ARE "MAER=17" PARTICULATE COMPONENTS:
!           1       AEROSOL 1 (NOMINALLY, BOUNDARY LAYER AEROSOL).
!           2       AEROSOL 2 (NOMINALLY, TROPOSPHERIC AEROSOL).
!           3       AEROSOL 3 (NOMINALLY, STRATOSPHERIC AEROSOL).
!           4       AEROSOL 4 (NOMINALLY, VOLCANIC AEROSOL).
!           5       CIRRUS CLOUD.
!           6       CLOUD 1 (NOMINALLY, WATER CLOUD).
!           7       CLOUD 2 (NOMINALLY, ICE CLOUD).
!           8-17    NOVAM (NAVY OCEANIC VERTICAL AEROSOL MODEL) LAYERS.
!       NAER     NUMBER OF ACTIVE AEROSOLS.
!       EXTV     SPECTRAL EXTINCTION (NORMALIZED TO 1 AT 550 NM).
!       ABSV     SPECTRAL ABSORPTION (1-ABSV/EXTV=SCATTERING ALBEDO).
!       ASYV     SPECTRAL LEGENDRE MOMENT (DIVIDED BY 2N+1).
!       FRAC5    5 CM-1 GRID SPECTRAL INTERPOLATION FRACTION.
!       ASYVLO   ASYMMETRY FACTOR FROM PREVIOUS SPECTRAL FREQUENCY.
      INTEGER NAER
      REAL EXTV,ABSV,ASYV,FRAC5,ASYVLO
      COMMON/AER/NAER,EXTV(MAER),ABSV(MAER),ASYV(MXCMU,MAER),           &
     &  FRAC5,ASYVLO(MAER)

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

!     /NAMEX/
!       CNAMEX   NAME OF CROSS-SECTION (X) SPECIES.
      CHARACTER CNAMEX*8
      COMMON/NAMEX/CNAMEX(NMOLX)
      SAVE /NAMEX/

!     LOCAL VARIABLES:
!       ISEG     PATH SEGMENT INDEX.
!       IW9SAP   INDEX USED TO WRITE 9 SPECTRAL POINTS AT A TIME.
!       IWVLIM   INDEX USED TO WRITE 9 SPECTRAL POINTS AT A TIME.
!       IWVSAP   SPECTRAL GRID INDEX FOR SPECTRAL AEROSOL PROFILES.
!       KSPEC    SPECIES LOOP INDEX.
      INTEGER ISEG,IW9SAP,IWVLIM,IWVSAP,KSPEC

!     WRITE TABLES:
      WRITE(IPR,'(////A,//(2A))')' LAYER ABSORBER AMOUNTS FOR THE '     &
     &  //'PATH SEGMENT ENDING AT Z',' ISEG    Z        TBAR       ',   &
     &  'HNO3         O3_UV     CNTMSLF1     CNTMSLF2      CNTMFRN  '   &
     &  //'      O2         WAT_DROP     ICE_PART','       (KM)    ',   &
     &  '    (K)     (ATM_CM)     (ATM_CM)    (MOL/CM2)    (MOL/CM2)'   &
     &  //'    (MOL/CM2)    (MOL/CM2)   (KM GM/M3)   (KM GM/M3)'
      WRITE(IPR,'(/(I3,0P,F10.4,F9.2,1P,6E13.3,0P,2F13.7))')            &
     &  (ISEG,PTHALT(ISEG),TSEG(ISEG),WPTH(11,ISEG),WPTH(8,ISEG),       &
     &  WPTH(5,ISEG),WPTH(9,ISEG),WPTH(10,ISEG),WPTH(58,ISEG),          &
     &  WPTH(66,ISEG),WPTH(67,ISEG),ISEG=1,NSEGX)

!     FETCH CLOUD 550 NM EXTINCTION COEFFICIENTS:
      CALL AEREXT(18181.82,6)
      WRITE(IPR,'(///2(A,5X,A,8X),4(A,8X),2(A,5X),A,/8X,A,16X,A,69X,A)')&
     &  ' ISEG','Z','N2_CONT','MOL_SCAT','AER_1','AER_2','AER_3',       &
     &  'AER_4','CIRRUS','WAT_DROP','ICE_PART','(KM)','(550nm EXT)',    &
     &  '(550 NM OPTICAL DEPTH)'
      WRITE(IPR,'(/(I3,0P,F10.4,9F13.7))')(ISEG,                        &
     &  PTHALT(ISEG),WPTH(4,ISEG),RAY550*WPTH(6,ISEG),WPTH(7,ISEG),     &
     &  WPTH(12,ISEG),WPTH(13,ISEG),WPTH(14,ISEG),WPTH(16,ISEG),        &
     &  WPTH(66,ISEG)*EXTV(6),WPTH(67,ISEG)*EXTV(7),ISEG=1,NSEGX)
      IF(LSAP)THEN
          WRITE(IPR,'(///A)')'  INCREMENTAL "SAP" AEROSOL EXTINCTION'   &
     &      //' OPTICAL DEPTHS FOR PATH SEGMENTS ENDING AT Z'
          DO IW9SAP=1,NWVSAP,9
              IWVLIM=MIN(IW9SAP+8,NWVSAP)
              WRITE(IPR,'(//A,9(F11.5,A))')' ISEG   Z(KM)',             &
     &          (WAVSAP(IWVSAP),'um',IWVSAP=IW9SAP,IWVLIM)
              DO ISEG=1,NSEGX
                  WRITE(IPR,'(I3,F10.4,9F13.7)')ISEG,PTHALT(ISEG),      &
     &              (XSAP(IWVSAP,ISEG),IWVSAP=IW9SAP,IWVLIM)
              ENDDO
          ENDDO
      ENDIF

!     PRINT PATH SUMMARY:
      WRITE(IPR,'(///A)')                                               &
     &  '    ISEG  Z       H2O       O3        CO2       CO        CH4' &
     &  //'       N2O       O2        NH3       NO        NO2       SO2'
      IF(MODTRN)THEN
          WRITE(IPR,'(8X,A,54X,A,45X,A)')'(KM)    (','ATM CM',')'
      ELSE
          WRITE(IPR,'(8X,A,44X,A,45X,A)')                               &
     &      '(KM)    (G/CM2)   (','ATM CM',')'
      ENDIF
      WRITE(IPR,'(/(I4,0P,F10.4,1P,11E10.2))')(ISEG,PTHALT(ISEG),       &
     &  WPTH(17,ISEG),WPTH(31,ISEG),WPTH(36,ISEG),WPTH(44,ISEG),        &
     &  WPTH(46,ISEG),WPTH(47,ISEG),WPTH(50,ISEG),WPTH(52,ISEG),        &
     &  WPTH(54,ISEG),WPTH(55,ISEG),WPTH(56,ISEG),ISEG=1,NSEGX)
      WRITE(IPR,'(///A14,13(1X,A8:),/(14X,13(1X,A8)))')                 &
     &  '  ISEG    Z   ',(CNAMEX(KSPEC),KSPEC=1,NMOLX)
      WRITE(IPR,'(8X,3(A:,53X))')'(KM)    (','ATM CM',')'
      DO ISEG=1,NSEGX
          WRITE(IPR,'(I4,F10.4,1P,13E9.2:,/(14X,13E9.2:))')             &
     &      ISEG,PTHALT(ISEG),(WPTH(MEXT+KSPEC,ISEG),KSPEC=1,NMOLX)
      ENDDO

!     RETURN TO GEO:
      RETURN
      END

      SUBROUTINE GEOTOT(H1ALT,H2ALT,ILOS,LSSGEO,MSPATH,LSAP,IERROR)

!     WRITES TOTAL PATH AMOUNTS:
      IMPLICIT NONE

!     PARAMETERS:
!       RAY550   550 NM RAYLEIGH (C6DTA) SCATTERING COEFFICIENT [KM-1].
      REAL RAY550
      PARAMETER(RAY550=.0121181)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       H1ALT    OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    FINAL (TANGENT FOR LIMB PATH) ALTITUDE [KM].
!       ILOS     PATH LINE-OF-SIGHT INDEX.
!       LSSGEO   LOGICAL FLAG, .TRUE. FOR SOLAR PATHS.
!       MSPATH   LOGICAL, TRUE FOR MULTIPLE SCATTERING VERTICAL PATH.
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      DOUBLE PRECISION H1ALT,H2ALT
      INTEGER ILOS
      LOGICAL LSSGEO,MSPATH,LSAP

!     OUTPUT ARGUMENTS:
!       IERROR   INTEGER ERROR FLAG.
      INTEGER IERROR

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'SEGDAT.h'
      INCLUDE 'YPROP.h'
      INCLUDE 'YPROPC.h'

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

!     /SAPDEP/
!       XSAP     LINE-OF-SIGHT PATH AEROSOL EXTINCTION OPTICAL DEPTHS.
!       XSAPM    VERTICAL MS PATH AEROSOL EXTINCTION OPTICAL DEPTHS.
!       ASAP     LINE-OF-SIGHT PATH AEROSOL ABSORPTION OPTICAL DEPTHS.
!       ASAPM    VERTICAL MS PATH AEROSOL ABSORPTION OPTICAL DEPTHS.
!       XSAPS    SOLAR PATH AEROSOL EXTINCTION OPTICAL DEPTHS FROM LOS.
!       XSAPSM   SOLAR PATH AEROSOL EXTINCTION OPTICAL DEPTHS FOR MS.
!       ASAPS    SOLAR PATH AEROSOL ABSORPTION OPTICAL DEPTHS FROM LOS.
!       ASAPSM   SOLAR PATH AEROSOL ABSORPTION OPTICAL DEPTHS FOR MS.
      REAL XSAP,XSAPM,ASAP,ASAPM,XSAPS,XSAPSM,ASAPS,ASAPSM
      COMMON/SAPDEP/                                                    &
     &  XSAP(1:MWVSAP,0:LAYTWO,0:MLOS),XSAPM(1:MWVSAP,0:LAYDIM),        &
     &  ASAP(1:MWVSAP,0:LAYTWO,0:MLOS),ASAPM(1:MWVSAP,0:LAYDIM),        &
     &  XSAPS(MWVSAP,LAYTWO,MLOS),XSAPSM(MWVSAP,LAYDIM),                &
     &  ASAPS(MWVSAP,LAYTWO,MLOS),ASAPSM(MWVSAP,LAYDIM)

!     /NAMEX/
!       CNAMEX   NAME OF CROSS-SECTION (X) SPECIES.
      CHARACTER CNAMEX*8
      COMMON/NAMEX/CNAMEX(NMOLX)
      SAVE /NAMEX/

!     LOCAL VARIABLES:
!       IW9SAP   INDEX USED TO WRITE 9 SPECTRAL POINTS AT A TIME.
!       IWVLIM   INDEX USED TO WRITE 9 SPECTRAL POINTS AT A TIME.
!       IWVSAP   SPECTRAL GRID INDEX FOR SPECTRAL AEROSOL PROFILES.
!       INAMLM   LIMIT ON INAM LOOP INDEX.
!       INAM     LOOP INDEX FOR X SPECIES NAME LIST.
      INTEGER IW9SAP,IWVLIM,IWVSAP,INAMLM,INAM

!     BEGIN WRITING:
      IF(H1ALT.LT.ZM(1))THEN
          WRITE(IPR,'(/A,F12.8,A,/15X,A,F12.8,A)')' Error in GEO: '//   &
     &      ' INITIAL ALTITUDE [H1ALT =',H1ALT,' KM] IS BELOW',         &
     &      ' THE BOTTOM OF THE ATMOSPHERE [ZM(1) =',ZM(1),' KM].'
          IERROR=1
          RETURN
      ELSEIF(H2ALT.LT.ZM(1) .AND. ITYPE.NE.1)THEN
          WRITE(IPR,'(/A,F12.8,A,/16X,A,F12.8,A)')' Error in GEO:  '//  &
     &      'FINAL OR TANGENT ALTITUDE [H2ALT =',H2ALT,' KM] IS',       &
     &      'BELOW THE BOTTOM OF THE ATMOSPHERE [ZM(1) =',ZM(1),' KM].'
          IERROR=1
          RETURN
      ENDIF
      IF(LSSGEO .OR. LJMASS)RETURN
      IF(.NOT.MODTRN)THEN
          WRITE(IPR,'(//A)')                                            &
     &      '   EQUIVALENT SEA LEVEL TOTAL ABSORBER AMOUNTS:'
      ELSEIF(MSPATH)THEN
          WRITE(IPR,'(//A)')'   TOTAL COLUMN ABSORBER'                  &
     &      //' AMOUNTS FOR A VERTICAL PATH FROM GROUND TO SPACE:'
      ELSE
          WRITE(IPR,'(//A,I3,A)')'   TOTAL COLUMN ABSORBER'             &
     &      //' AMOUNTS FOR LINE-OF-SIGHT PATH No.',ILOS,':'
      ENDIF
      WRITE(IPR,'(2(/15X,2A),/10X,1P,6E12.4,0P,2(1X,F11.8),             &
     &  /(/2(/15X,2A),/11X,4(1X,F11.8),3(F11.6,1X),F11.4))')            &
     &  '  HNO3      O3 UV      CNTMSLF1    CNTMSLF2    CNTMFRN  ',     &
     &  '   N2 CONT     MOL SCAT   TOTAL AER',                          &
     &  '(ATM CM)   (ATM CM)   (MOL CM-2)  (MOL CM-2)  (MOL CM-2)',     &
     &  '              ( 550 NM EXTINCTION )',WTOTAL(11),WTOTAL(8),     &
     &  WTOTAL(5),WTOTAL(9),WTOTAL(10),WTOTAL(4),RAY550*WTOTAL(6),      &
     &  WTOTAL(7)+WTOTAL(12)+WTOTAL(13)+WTOTAL(14),                     &
     &  ' AER 1       AER 2       AER 3       AER 4      CIRRUS   ',    &
     &  '  WAT DROP    ICE PART   MEAN AER RH',                         &
     &  ' (                 550 NM EXTINCTION                 )   ',    &
     &  ' (KM GM/M3)  (KM GM/M3)   (PERCENT)',                          &
     &  WTOTAL(7),WTOTAL(12),WTOTAL(13),WTOTAL(14),WTOTAL(16),          &
     &  WTOTAL(66),WTOTAL(67),WTOTAL(15),' H2O   ',                     &
     &  '      O3          CO2         CO          CH4         N2O'
      IF(MODTRN)THEN
          WRITE(IPR,'(13X,A,2(35X,A))')'(','ATM CM',')'
      ELSE
          WRITE(IPR,'(13X,A,2(22X,A))')'(G/CM**2)     (','ATM CM',')'
      ENDIF
      WRITE(IPR,'((10X,1P,6E12.4,//2(/15X,A),2(22X,A)))')WTOTAL(17),    &
     &  WTOTAL(31),WTOTAL(36),WTOTAL(44),WTOTAL(46),WTOTAL(47),         &
     &  ' O2          NH3         NO          NO2         SO2',         &
     &  '(','ATM CM',')',WTOTAL(50),WTOTAL(52),WTOTAL(54),WTOTAL(55),   &
     &  WTOTAL(56)
      IF(LSAP)THEN
          IF(MSPATH)THEN
              DO IW9SAP=1,NWVSAP,9
                  IWVLIM=MIN(IW9SAP+8,NWVSAP)
                  WRITE(IPR,'(//11X,A,2(24X,A),/10X,9(1X,F9.5,A))')     &
     &              '(','SPECTRAL AEROSOL PROFILES (SAP)'               &
     &              //' EXTINCTION OPTICAL DEPTHS',')',                 &
     &              (WAVSAP(IWVSAP),'um',IWVSAP=IW9SAP,IWVLIM)
                  WRITE(IPR,'(10X,9F12.6)')                             &
     &              (XSAPM(IWVSAP,0),IWVSAP=IW9SAP,IWVLIM)
              ENDDO
          ELSE
              DO IW9SAP=1,NWVSAP,9
                  IWVLIM=MIN(IW9SAP+8,NWVSAP)
                  WRITE(IPR,'(//11X,A,2(24X,A),/10X,9(1X,F9.5,A))')     &
     &              '(','SPECTRAL AEROSOL PROFILES (SAP)'               &
     &              //' EXTINCTION OPTICAL DEPTHS',')',                 &
     &              (WAVSAP(IWVSAP),'um',IWVSAP=IW9SAP,IWVLIM)
                  WRITE(IPR,'(10X,9F12.6)')                             &
     &              (XSAP(IWVSAP,0,ILOS),IWVSAP=IW9SAP,IWVLIM)
              ENDDO
          ENDIF
      ENDIF
      INAMLM=7
   10 CONTINUE
      IF(NMOLX.GT.INAMLM)THEN
          WRITE(IPR,'(//8X,7(4X,A),/13X,A,2(36X,A),                     &
     &      /10X,1P,7E12.4)')(CNAMEX(INAM),INAM=INAMLM-6,INAMLM),       &
     &      '(','ATM CM',')',(WTOTAL(MEXT+INAM),INAM=INAMLM-6,INAMLM)
          INAMLM=INAMLM+7
          GOTO 10
      ENDIF
      IF(NMOLX.GE.INAMLM-6)THEN
          WRITE(IPR,'(//8X,7(4X,A))')(CNAMEX(INAM),INAM=INAMLM-6,NMOLX)
          WRITE(IPR,'(13X,A,2(36X,A),/10X,1P,7E12.4)')                  &
     &      '(','ATM CM',')',(WTOTAL(MEXT+INAM),INAM=INAMLM-6,NMOLX)
      ENDIF
      INAMLM=7
   20 CONTINUE
      IF(NMOLY.GT.INAMLM)THEN
          WRITE(IPR,'(//8X,7(4X,A),/13X,A,2(36X,A),                     &
     &      /10X,1P,7E12.4)')(CNAMEY(INAM),INAM=INAMLM-6,INAMLM),       &
     &      '(','ATM CM',')',(WTOTAL(MEXTX+INAM),INAM=INAMLM-6,INAMLM)
          INAMLM=INAMLM+7
          GOTO 20
      ENDIF
      IF(NMOLY.GE.INAMLM-6)THEN
        WRITE(IPR,'(//8X,7(4X,A))')(CNAMEY(INAM),INAM=INAMLM-6,NMOLY)
        WRITE(IPR,'(13X,A,2(36X,A),/10X,1P,7E12.4)')                    &
     &    '(','ATM CM',')',(WTOTAL(MEXTX+INAM),INAM=INAMLM-6,NMOLY)
      ENDIF

!     RETURN TO GEO:
      RETURN
      END
