      SUBROUTINE LOOP(VBAND5,NO_WRT,IBINPT,VCEN,IPH,SUMSSS,LVBND5,UNIF, &
     &  TRACE,THMCUM,S0,TSNOBS,TSNREF,KNTRVL,GROUND,LABEL)

!     LOOP LOOPS OVER RADIANCE PATH SEGMENTS (NO MS) AT EACH FREQUENCY.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       IBINPT   BIN NUMBER OF CURRENT SPECTRAL POINT.
!                (CENTER FREQUENCY = IBINPT * BNDWID + OSHIFT).
!       VBAND5   5 CM-1 RESOLUTION DATA SPECTRAL FREQUENCY [CM-1].
!       NO_WRT   OUTPUT FLAG FOR SPECTRAL DATA AT CURRENT FREQUENCY.
!       VCEN     SPECTRAL-BIN CENTRAL FREQUENCY [CM-1].
!       LVBND5   LOGICAL FLAG, TRUE WHEN 5 CM-1 DATA IS TO BE OBTAINED.
!       KNTRVL   NUMBER OF CORRELATED-K SUB-INTERVALS.
!       GROUND   LOGICAL FLAG, TRUE IF LINE OF SIGHT INTERSECTS GROUND.

!     OUTPUT ARGUMENTS:
!       THMCUM   PATH EMIT & SCAT THERMAL RADIANCE [W CM-2 SR-1 / CM-1].
      INTEGER IBINPT,IPH,KNTRVL
      REAL VBAND5,VCEN,SUMSSS,UNIF,TRACE,THMCUM,S0,TSNOBS,TSNREF
      LOGICAL NO_WRT,LVBND5,GROUND
      CHARACTER LABEL*6

!     COMMONS:
      INCLUDE 'YPROP.h'
      INCLUDE 'BASE.h'
      INCLUDE 'SEGDAT.h'
      INCLUDE 'IFIL.h'
      INCLUDE 'WSOL.h'

!     /RCNSTN/
!       PI       THE CONSTANT PI.
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

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

!     /SURFWV/
!       LAMBER  LOGICAL FLAG, .TRUE. FOR LAMBERTIAN SURFACE.
!       TPTEMP  TARGET-PIXEL SURFACE TEMPERATURES [K].
!       TPHDIR  TARGET-PIXEL HEMISPHERE DIRECTIONAL REFLECTANCE AT
!               VIEWING ANGLES.
!       TPBRDF  TARGET-PIXEL BIDIRECTIONAL REFLECTANCE DISTRIBUTION
!               FUNCTION AT VIEWING AND SUN ANGLES.
!       AATEMP  AREA-AVERAGED GROUND SURFACE TEMPERATURES [K].
!       AASALB  AREA-AVERAGED GROUND SURFACE ALBEDO.
!       AADREF  AREA-AVERAGED GROUND SURFACE DIRECTIONAL REFLECTIVITY
!               AT THE SOLAR ZENITH ANGLE.
!       EMU     GROUND DIRECTIONAL EMISSIVITY AT VIEWING ANGLES.
!       BEM     GROUND DIRECTIONAL EMISSIVITY AT QUADRATURE ANGLES.
!       RMU     GROUND BRDF AZIMUTH COMPONENTS AT VIEWING ANGLES
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
!       BDR     GROUND BRDF AZIMUTH COMPONENTS AT QUADRATURE ANGLE
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
      LOGICAL LAMBER
      REAL TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,EMU,BEM,RMU,BDR
      COMMON/SURFWV/LAMBER,TPTEMP(MLOS),TPHDIR(MLOS),TPBRDF(MLOS),      &
     &  AATEMP,AASALB,AADREF,EMU(MXUMU),BEM(MI),                        &
     &  RMU(1:MXUMU,0:MI,0:MAZ),BDR(1:MI,0:MI,0:MAZ)                    &

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

!     /MSRD/
!       CSSCAT   COSINE OF THE SCATTERING ANGLE.
!                (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       SLEGEN   Nth LEGENDRE POLYNOMIAL EVALUATED AT THE COSINE OF THE
!                SCATTERING ANGLE TIMES (2N+1)/4pi (N=0 TO NSTR-1).
!       CSZEN0   LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.
!       CSZEN    LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.
!       CSZENX   AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING
!                (AWAY FROM EARTH) THE CURRENT LAYER.
!       BBGRND   THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].
!       BBNDRY   LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1].
!       TCONT    LAYER CONTINUUM OPTICAL DEPTH.
!       TAUT     LAYER TOTAL VERTICAL EXTINCTION OPTICAL DEPTH.
!       GTSCAT   SUM OVER SCATTERING SOURCES OF SCATTERING OPTICAL DEPTH
!                AND PHASE FUNCTION LEGENDRE COEFFICIENT PRODUCTS.
!       COSBAR   LAYER EFFECTIVE SCATTERING ASYMMETRY FACTOR.
!       DEPRAT   FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.
!       S0DEP    OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.
!       S0TRN    TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]
!       UPF      LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].
!       DNF      LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].
!       UPFS     LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].
!       DNFS     LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].
!       CO_LIN   TRUE IF LOS AND SOLAR PATHS ARE NEARLY IDENTICAL.
      REAL CSSCAT,SLEGEN,CSZEN0,CSZEN,CSZENX,TCONT,TAUT,GTSCAT,COSBAR,  &
     &  BBGRND,BBNDRY,S0DEP,S0TRN,DEPRAT,UPF,DNF,UPFS,DNFS
      LOGICAL CO_LIN
      COMMON/MSRD/CSSCAT(MLOS),SLEGEN(0:MAZ,MLOS),CSZEN0(LAYDIM),       &
     &  CSZEN(LAYDIM),CSZENX(LAYDIM),TCONT(LAYDIM),TAUT(MXKSUB,LAYDIM), &
     &  GTSCAT(0:MXCMU,1:LAYDIM),COSBAR(LAYDIM),BBGRND,BBNDRY(LAYDIM),  &
     &  S0DEP(MXKSUB,LAYTWO),S0TRN(MXKSUB,LAYTWO),DEPRAT(MXKSUB,LAYDIM),&
     &  UPF(MXKSUB,LAYDIM),DNF(MXKSUB,LAYDIM),UPFS(MXKSUB,LAYDIM),      &
     &  DNFS(MXKSUB,LAYDIM),CO_LIN(MLOS)

!     /SEG5DT/
!       SPCHI    LINE-OF-SIGHT SPECTRAL DATA @ HIGHER 5 CM-1 GRID POINT.
!       SPCMHI   VERTICAL PATH SPECTRAL DATA @ HIGHER 5 CM-1 GRID POINT.
!       SPCLO    LINE-OF-SIGHT SPECTRAL DATA @ LOWER 5 CM-1 GRID POINT.
!       SPCMLO   VERTICAL PATH SPECTRAL DATA @ LOWER 5 CM-1 GRID POINT.
!            1   SEGMENT SCATTERING OD WEIGHTED ASYMMETRY FACTOR
!            2   SEGMENT AEROSOL SCATTERING OPTICAL DEPTH
!            3   TOTAL O2 CONTINUUM TRANSMITTANCE.
!            4   N2 CONTINUUM TRANSMITTANCE.
!            5   TOTAL H2O CONTINUUM TRANSMITTANCE.
!            6   RAYLEIGH MOLECULAR SCATTERED TRANSMITTANCE.
!            7   TRANSMITTANCE FROM AEROSOL EXTINCTION.
!            8   TOTAL OZONE CONTINUUM TRANSMITTANCE.
!            9   ALL CONTINUUM TRANSMITTANCES EXCEPT O2 AND HNO3.
!           10   TRANSMITTANCE FROM AEROSOL ABSORPTION.
!           11   HNO3 TRANSMITTANCE.
!           12   MOLECULAR CONTINUUM OPTICAL DEPTH
!           13   SEGMENT AEROSOL EXTINCTION OPTICAL DEPTH
!           14   COMBINED SPECIES INCREMENTAL CONTINUUM OPTICAL DEPTH
!           15   RAYLEIGH MOLECULAR SCATTERING OPTICAL DEPTH
!           16   CIRRUS CLOUD TRANSMITTANCE (ICLD=20 ONLY)
!           17   UV/VIS NO2 TRANSMISSION
!           18   UV/VIS SO2 TRANSMISSION
!           19   SEGMENT WATER DROPLET SCATTERING OPTICAL DEPTH
!           20   SEGMENT ICE PARTICLE SCATTERING OPTICAL DEPTH
!           21   SEGMENT BOUNDARY LAYER AEROSOL SCATTERING OPTICAL DEPTH
!           22   SEGMENT TROPOSPHERIC AEROSOL SCATTERING OPTICAL DEPTH
!           23   SEGMENT STRATOSPHERIC AEROSOL SCATTERING OPTICAL DEPTH
!           24   SEGMENT MESOSPHERIC+ AEROSOL SCATTERING OPTICAL DEPTH
!           25   SEGMENT NOVAM AEROSOL SCATTERING OPTICAL DEPTH
!           26   SEGMENT NOVAM AEROSOL ASYMMETRY FACTOR
!           27   SEGMENT STD OR SUB-VIS CIRRUS SCATTERING OPTICAL DEPTH
!           28   SEGMENT CLOUD EXTINCTION OPTICAL DEPTH
!           29   SEGMENT CLOUD SCATTERING OD WEIGHTED ASYMMETRY FACTOR
!           30   SEGMENT RAIN EXTINCTION OPTICAL DEPTH
!           31   SEGMENT RAIN SCATTERING OPTICAL DEPTH
!           32   SEGMENT RAIN ASYMMETRY FACTOR
!           33   "14" WITHOUT RAIN, CLOUD AND AEROSOL.
!      NSEG5-2   VIS/NIR CH4 OPTICAL DEPTH
!      NSEG5-1   H2-H2 DIMER OPTICAL DEPTH
!        NSEG5   H2-HE DIMER OPTICAL DEPTH
!       SSAPH    LOS SAP SCATTERING OPTICAL DEPTH @ HIGHER 5 CM-1 POINT.
!       SSAPHM   MS SAP SCATTERING OPTICAL DEPTH @ HIGHER 5 CM-1 POINT.
!       SSAPL    LOS SAP SCATTERING OPTICAL DEPTH @ LOWER 5 CM-1 POINT.
!       SSAPLM   MS SAP SCATTERING OPTICAL DEPTH @ LOWER 5 CM-1 POINT.
!       TXLEG    SPECTRALLY INTERPOLATED LEGENDRE COEFFICIENTS FOR THE
!                SCATTERING PHASE FUNCTION (OPTICAL DEPTH WEIGHTED SUM).
      REAL SPCHI,SPCMHI,SPCLO,SPCMLO,SSAPH,SSAPHM,SSAPL,SSAPLM,TXLEG
      COMMON/SEG5DT/SPCHI(NSEG5,LAYTWO,MLOS,3),SPCMHI(NSEG5,LAYDIM,3),  &
     &  SPCLO(NSEG5,LAYTWO,MLOS,3),SPCMLO(NSEG5,LAYDIM,3),              &
     &  SSAPH(LAYTWO,MLOS,3),SSAPHM(LAYDIM,3),SSAPL(LAYTWO,MLOS,3),     &
     &  SSAPLM(LAYDIM,3),TXLEG(2:MXCMU,1:LAYDIM,0:1)
      SAVE/SEG5DT/

!     /CORKDT/
!       WTKSUB   SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.
!       WTKSAV   SAVED SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.
!       DEPLAY   INCREMENTAL OPTICAL DEPTHS.
!       TRNLAY   INCREMENTAL TRANSMITTANCES.
!       TRNCUM   CUMULATIVE TRANSMITTANCES.
!       K2TAIL   POINTER FROM K BIN TO LINE TAIL SUB-BIN
!                (=0 IF MULTIPLE LINE TAIL SUB-BINS CONTRIBUTE).
!       CONTWT   WEIGHTS FOR PARTITIONING LINE TAILS INTO K'S
!                (ONLY USED IF K2TAIL IS 0).
      REAL WTKSUB(MXKSUB),WTKSAV(NTLSUB),DEPLAY(MXKSUB),                &
     &  TRNLAY(MXKSUB),TRNCUM(MXKSUB),CONTWT(NTLSUB,MXKSUB)
      INTEGER K2TAIL(MXKSUB)
      COMMON/CORKDT/K2TAIL,WTKSUB,WTKSAV,DEPLAY,TRNLAY,TRNCUM,CONTWT
      SAVE /CORKDT/

!       SUBINT   SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.
!       UPFLX    LAYER BOUNDARY UPWARD THERMAL SPECTRAL
!                FLUX (INCLUDES SCATTERED SOLAR IF DISORT
!                & NO AZIMUTH DEPENDENCE) [W CM-2 / CM-1].
!       DNFLX    LAYER BOUNDARY DOWNWARD THERMAL SPECTRAL
!                FLUX (INCLUDES SCATTERED SOLAR IF DISORT
!                & NO AZIMUTH DEPENDENCE) [W CM-2 / CM-1].
!       UPFLXS   LAYER BOUNDARY UPWARD SCATTERED SOLAR
!                SPECTRAL FLUX (USED WITH DISORT ONLY
!                IF AZIMUTH DEPENDENT) [W CM-2 / CM-1].
!       DNFLXS   LAYER BOUNDARY DOWNWARD SCATTERED SOLAR
!                SPECTRAL FLUX (USED WITH DISORT ONLY
!                IF AZIMUTH DEPENDENT) [W CM-2 / CM-1].
!       NTFLX    LAYER BOUNDARY NET (THERMAL PLUS SCATTERED SOLAR PLUS
!                DIRECT SOLAR) UPWARD SPECTRAL FLUX [W CM-2 / CM-1].
      REAL SUBINT,UPFLX,DNFLX,UPFLXS,DNFLXS,NTFLX
      COMMON/NETFLX/SUBINT(MXKSUB),UPFLX(LAYDIM),DNFLX(LAYDIM),         &
     &  UPFLXS(LAYDIM),DNFLXS(LAYDIM),NTFLX(LAYDIM)

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER AMOD3D*1
      COMMON/CJM5/AMOD3D

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

!     /LSTART/
!       TOLABS   OUTPUT WARNING THRESHOLD FOR EXCESS ABSORPTION.
!       LLORMN   OUTPUT WARNING FLAG FOR SMALL LORENTZ HALF WIDTH.
!       LLORMX   OUTPUT WARNING FLAG FOR LARGE LORENTZ HALF WIDTH.
!       LDOPMN   OUTPUT WARNING FLAG FOR SMALL DOPPLER HALF WIDTH.
!       LDOPMX   OUTPUT WARNING FLAG FOR LARGE DOPPLER HALF WIDTH.
!       LLINMN   OUTPUT WARNING FLAG FOR SMALL EFFECTIVE LINE NUMBER.
!       LSSA     OUTPUT WARNING FLAG FOR SCATTERING ALBEDO NEAR 1.
!       LAEREX   OUTPUT WARNING FLAG FOR ROUTINE AEREXT.
!       LBMCRK   OUTPUT WARNING FLAG FOR ROUTINE BMCRKS.
!       LFLUXS   OUTPUT WARNING FLAG FOR ROUTINE FLUXES.
!       LSSRAD   OUTPUT WARNING FLAG FOR ROUTINE SSRAD.
!       LTRLAY   OUTPUT WARNING FLAG FOR ROUTINE TRLAY.
!       LNVERS   OUTPUT WARNING FLAG FOR ROUTINE DENFAC.
!       LGEOM    OUTPUT WARNING FLAG FOR INDEX OF REFRACTION GRADIENT.
!       LO3TRN   OUTPUT WARNING FLAG FOR OZONE CURTIS-GODSON PROBLEM.
      REAL TOLABS
      LOGICAL LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,                  &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      COMMON/LSTART/TOLABS,LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,     &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      SAVE /LSTART/

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      REAL BBFN

!     /K_FRMT/
!       TFRMT    FORMAT FOR OUTPUT OF K-DEPENDENT TRANSMITTANCES.
!       RFRMT    FORMAT FOR OUTPUT OF K-DEPENDENT RADIANCES.
!       HTFRMT   FORMAT FOR OUTPUT OF K-DEPENDENT TRANSMITTANCES HEADER.
!       HRFRMT   FORMAT FOR OUTPUT OF K-DEPENDENT RADIANCES HEADER.
      CHARACTER TFRMT*33,RFRMT*38,HTFRMT*36,HRFRMT*34
      COMMON/K_FRMT/TFRMT,RFRMT,HTFRMT,HRFRMT
      SAVE /K_FRMT/

!     /LOOP_0/
!       O3MAX    TOLERANCE FOR 3-PARAMETER CURTIS-GODSON MONOTONICITY.
!       WDKDIS   BAND MODEL WIDTH [CM-1].
!       V0KDIS   CALCULATION INITIAL FREQUENCY.
!       NMWAVE   WAVELENGTH BIN OF PREVIOUS DATA [NM].
!       NTHEAD   HEADER RECORD FOR K-DEPENDENT TRANSMITTANCE FILE.
!       NRHEAD   HEADER RECORD FOR K-DEPENDENT RADIANCE FILE.
!       N4BYTE   NUMBER OF 4 BYTE VARIABLES PER RECORD.
!       NTREC    NUMBER OF RECORDS PER SPECTRAL POINT IN TRANSM FILE.
!       NRREC    NUMBER OF RECORDS PER SPECTRAL POINT IN RADIANCE FILE.
!       NRANGE   NUMBER OF OUTPUT PATH RANGES.
!       KPRINT   LOGICAL FLAG, DICTATING K-DEPENDENT DATA OUTPUT.
      REAL O3MAX,WDKDIS,V0KDIS
      INTEGER NMWAVE,NTHEAD,NRHEAD,N4BYTE,NTREC,NRREC,NRANGE
      LOGICAL KPRINT
      COMMON/LOOP_0/O3MAX,WDKDIS,V0KDIS,                                &
     &  NMWAVE,NTHEAD,NRHEAD,N4BYTE,NTREC,NRREC,NRANGE,KPRINT
      SAVE /LOOP_0/

!     /ANGSRF/
!       CVWSRF  COSINE OF THE VIEW ZENITH ANGLE FROM THE SURFACE.
!       CSNSRF  COSINE OF THE SOLAR (LUNAR) ZENITH AT SURFACE.
!       AZMSRF  RELATIVE AZIMUTH ANGLE (SUN - SENSOR AT SURFACE) [RAD].
!       UMU1    COSINE OF THE PATH NADIR ANGLE.
!               (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       UMU0    COSINE OF THE SOLAR ZENITH ANGLE.
!               (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       PHI1    RELATIVE AZIMUTH ANGLE (SUN - LOS PATH AT SENSOR) [DEG].
!               (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       CMU     COSINE OF THE NADIR ANGLES USED IN DISORT.
      REAL CVWSRF,CSNSRF,AZMSRF,UMU1,UMU0,PHI1,CMU
      COMMON/ANGSRF/CVWSRF(MLOS),CSNSRF(MLOS),AZMSRF(MLOS),UMU1(MLOS),  &
     &  UMU0,PHI1(MLOS),CMU(MI)

!     LOCAL VARIABLES:
!       IPATH    PATH TYPE LABEL (1 FOR SENSOR TO SUN, 2 FOR SCATTERING
!                POINT TO SUN, AND 3 FOR SENSOR TO SCATTERING POINT).
!       PFMOL0   MOLECULAR PHASE FUNCTION ISOTROPIC SCATTERING TERM.
!       PFMOL2   MOLECULAR PHASE FUNCTION ANISOTROPIC SCATTERING TERM.
!       P2RAY    SECOND LEGENDRE EXPANSION COEFFICIENT OVER 5 (= 2N+1).
!       INTRVL   LOOP INDEX FOR CORRELATED-K SUB-INTERVALS.
!       JNTRVL   NUMBER OF K-DISTRIBUTION INTERVALS FOR CURRENT FREQ.
!       ITLSUB   LOOP INDEX FOR LINE TAIL SUB-INTERVALS.
!       TL_H2O   H2O TRANSMITTANCE FOR L-SHAPE PATH.
!       TL_UMX   UNIFORMLY MIXED GASES TRANSMITTANCE FOR L-SHAPE PATH.
!       TL_O3    O3 TRANSMITTANCE FOR L-SHAPE PATH.
!       TL_AUX   AUXILIARY SPECIES (Y) TRANSMITTANCE FOR L-SHAPE PATH.
!       TRAN_X   COMBINED TRANSMITTANCE OF CROSS-SECTION (X) MOLECULES.
!       TRAN_Y   COMBINED TRANSMITTANCE OF AUXILIARY (Y) MOLECULES.
!       H2OCUM   CUMULATIVE DIRECT PATH H2O CONTINUUM TRANSMITTANCE.
!       UMXCUM   CUMULATIVE DIRECT PATH CONTINUUM TRANSMITTANCE OF
!                UNIFORMLY MIXED GASES.
!       O3CUM    CUMULATIVE DIRECT PATH O3 CONTINUUM TRANSMITTANCE.
!       T_CH4    METHANE SHORTER WAVE TRANSMITTANCE.
!       GREY     GREY-BODY EMISSION [W SR-1 CM-2 / CM-1].
!       SURFAC   DIRECT SOLAR REFLECTION SURFACE CONTRIBUTION [SR-1].
!       SSAP     SCATTERING OPTICAL DEPTH OF SPECTRAL AEROSOL PROFILE.
!       RADCUM   CUMULATIVE PATH RADIANCE K-DATA [W SR-1 CM-2 / CM-1].
!       LBMCEN   FLAG, TRUE IF S/d & 1/d DATA FOR CURRENT FREQUENCY BIN.
!       LBMFLG   FLAG, TRUE IF CENTER OR TAIL DATA FOR CURRENT FREQ BIN.
      INTEGER K,IPATH,ISEG,ISEGP1,INTRVL,JNTRVL,ITLSUB,IEXT
      REAL STORE,BLAYER,BBOUND,THMLAY,TRNOLD,TRNNEW,TX9LAY,TX9CUM,TMOL, &
     &  OPTDEP,EMSLAY,PFMOL0,PFMOL2,P2RAY,TRAN_X,TRAN_Y,TX14SV,TL_H2O,  &
     &  TL_UMX,TL_O3,TL_AUX,H2OCUM,UMXCUM,O3CUM,T_CH4,GREY,SURFAC,SSAP, &
     &  O3TRAN,RADCUM(MXKSUB)
      DOUBLE PRECISION SMTRNL,SMTHML
      LOGICAL LBMCEN,LBMFLG

!     DEFINE THE LAYER INDEPENDENT 5 CM-1 DATA
      IF(LVBND5)CALL FRQ5DT(.FALSE.,VBAND5)
      IPATH=1

!     LINE-OF-SIGHT DEPENDENT INITIALIZATIONS:
      THMCUM=0.

!     FOR EACH WAVENUMBER, BMOD0 PERFORMS INITIALIZATIONS FOR BMOD:
      IF(MODTRN)CALL BMOD0(IBINPT,VCEN,.TRUE.,LBMCEN,LBMFLG)

!     INITIALIZE PARAMETERS:
      IPATH=1
      TRNOLD=1.
      O3TRAN=1.
      TX9CUM=1.

!     SKIP CORRELATED-K IF NO BAND MODEL LINE CENTER ABSORPTION:
      IF(KNTRVL.EQ.1)THEN

!         NO BAND MODEL DATA:  NO K-DISTRIBUTION USED.
          JNTRVL=1
          SUBINT(1)=1.
          WTKSUB(1)=1.
      ELSEIF(LBMCEN)THEN

!         BAND MODEL LINE CENTER DATA:  USE FULL K-DISTRIBUTION.
          JNTRVL=KNTRVL
          DO ITLSUB=1,NTLSUB
              SUBINT(ITLSUB)=WTKSAV(ITLSUB)
              WTKSUB(ITLSUB)=WTKSAV(ITLSUB)
          ENDDO
      ELSEIF(LBMFLG)THEN

!         BAND MODEL LINE TAIL DATA ONLY:  USE TAIL K-DISTRIBUTION.
          JNTRVL=NTLSUB
          DO ITLSUB=1,NTLSUB
              SUBINT(ITLSUB)=RTLSUB
              WTKSUB(ITLSUB)=RTLSUB
          ENDDO
      ELSE

!         NO BAND MODEL DATA:  NO K-DISTRIBUTION USED.
          JNTRVL=1
          SUBINT(1)=1.
          WTKSUB(1)=1.
      ENDIF
      DO INTRVL=1,JNTRVL
          TRNCUM(INTRVL)=1.
          RADCUM(INTRVL)=0.
      ENDDO

!     RAYLEIGH SCATTERING PHASE FUNCTION SPECTRAL COEFFICIENTS:
      CALL MOLSCT(VCEN,PFMOL0,PFMOL2,P2RAY)
      IF(AMOD3D.EQ.'M')THEN
          H2OCUM=1.
          UMXCUM=1.
          O3CUM=1.
      ENDIF

!     K-DISTRIBUTION DEPENDENT OUTPUT HEADER:
      IF(KPRINT .AND. .NOT.NO_WRT)THEN
          IF(BINOUT)THEN
              NTKDIS=NTKDIS+1
              WRITE(JTKDIS,REC=NTKDIS)VCEN,JNTRVL,PTHALT(0,1),          &
     &          (SUBINT(INTRVL),INTRVL=1,JNTRVL)
              NRKDIS=NRKDIS+1
              WRITE(JRKDIS,REC=NRKDIS)VCEN,JNTRVL,PTHALT(0,1),          &
     &          (SUBINT(INTRVL),INTRVL=1,JNTRVL)
          ELSE
              WRITE(HTFRMT(28:29),'(I2)')JNTRVL
              WRITE(ITKDIS,HTFRMT)VCEN,JNTRVL,PTHALT(0,1),              &
     &          'WEIGHTS:',(SUBINT(INTRVL),INTRVL=1,JNTRVL)
              WRITE(HRFRMT(25:26),'(I2)')JNTRVL
              WRITE(IRKDIS,HRFRMT)VCEN,JNTRVL,PTHALT(0,1),              &
     &          'WEIGHTS:',(SUBINT(INTRVL),INTRVL=1,JNTRVL)
          ENDIF
      ENDIF

!     BEGINNING OF LAYER LOOP:
      DO ISEG=1,NSEG(1)

!         PASS APPROPRIATE ABSORBER AMOUNTS:
          IF(IEMSCT.EQ.1)GOTO 20

!         FOR FIRST LAYER (ISEG=1) SET UP FOR SUN-TO-OBSERVER PATH:
          IF(IPATH.EQ.1)THEN
              IF(WSPTH(36,1,1).GE.0.)THEN

!                 DEFINE ALTITUDE DEPENDENT 5 CM-1 DATA:
                  IF(LVBND5)CALL SEG5(VBAND5,ISEG,1,NSEG(1),1,          &
     &              WSPTH(0,1,1),XSAPS(1,1,1),ASAPS(1,1,1),.FALSE.)
                  GOTO 30
              ENDIF

!             SCATTERING POINT IS IN SHADE IF WSPTH(36,1,1) IS NEGATIVE:
              CALL SHADE(IPH,1,ISEG,1,JNTRVL,VCEN,PFMOL0,PFMOL2,        &
     &          TSNOBS,TSNREF,SUMSSS,RADCUM)
          ENDIF
   10     CONTINUE

!         SUN-TO-SCATTERING_POINT PATHS:
          IPATH=2
          ISEGP1=ISEG+1
          IF(WSPTH(36,ISEGP1,1).GE.0.)THEN

!             DEFINE ALTITUDE DEPENDENT 5 CM-1 DATA [WHEN CORRELATED-K
!             APPROACH IS NOT INVOKED (JNTRVL=1), MOST ABSORBER AMOUNTS
!             INCLUDE OBSERVER-TO-SCATTERING_POINT CONTRIBUTIONS]:
              IF(LVBND5)CALL SEG5(VBAND5,ISEG,2,                        &
     &          NSEG(1),1,WSPTH(0,ISEGP1,1),                            &
     &          XSAPS(1,ISEGP1,1),ASAPS(1,ISEGP1,1),.FALSE.)
              GOTO 30
          ENDIF

!         SCATTERING POINT IS IN SHADE IF WSPTH(36,ISEGP1,1) < 0:
          CALL SHADE(IPH,1,ISEG,2,JNTRVL,VCEN,PFMOL0,PFMOL2,            &
     &      TSNOBS,TSNREF,SUMSSS,RADCUM)

!         DEFINE OPTICAL PATH ALTITUDE DEPENDENT 5 CM-1 DATA:
   20     CONTINUE
          IPATH=3
          IF(LVBND5)CALL SEG5(VBAND5,ISEG,3,NSEG(1),1,                  &
     &      WPTH(0,ISEG,1),XSAP(1,ISEG,1),ASAP(1,ISEG,1),.FALSE.)

!         DEFINE TX ARRAY
   30     CONTINUE

!         INTERPOLATE COARSE SPECTRAL RESOLUTION DATA:
          IF(ABS(FRAC5).LT..001)THEN

!             GRID POINT: NO SPECTRAL INTERPOLATION:
              DO K=1,16
                  TX(K)=SPCHI(K,ISEG,1,IPATH)
              ENDDO
              DO K=17,NSEG5
                  TX(47+K)=SPCHI(K,ISEG,1,IPATH)
              ENDDO

!             VIS/NIR CH4:
              T_CH4=EXP(-SPCHI(NSEG5-2,ISEG,1,IPATH))

!             H2-H2 AND H2-HE DIMERS:
              TX(MEXTX-1)=EXP(-SPCHI(NSEG5-1,ISEG,1,IPATH))
              TX(MEXTX)=EXP(-SPCHI(NSEG5,ISEG,1,IPATH))

!             SAP AEROSOLS:
              SSAP=SSAPH(ISEG,1,IPATH)
          ELSE

!             PERFORM SPECTRAL INTERPOLATION:
              DO K=1,16
                  STORE=SPCHI(K,ISEG,1,IPATH)
                  TX(K)=STORE+FRAC5*(SPCLO(K,ISEG,1,IPATH)-STORE)
              ENDDO
              DO K=17,NSEG5
                  STORE=SPCHI(K,ISEG,1,IPATH)
                  TX(47+K)=STORE+FRAC5*(SPCLO(K,ISEG,1,IPATH)-STORE)
              ENDDO

!             VIS/NIR CH4:
              STORE=SPCHI(NSEG5-2,ISEG,1,IPATH)
              T_CH4                                                     &
     &          =EXP(-(STORE+FRAC5*(SPCLO(NSEG5-2,ISEG,1,IPATH)-STORE)))

!             H2-H2 AND H2-HE DIMERS:
              STORE=SPCHI(NSEG5-1,ISEG,1,IPATH)
              TX(MEXTX-1)                                               &
     &          =EXP(-(STORE+FRAC5*(SPCLO(NSEG5-1,ISEG,1,IPATH)-STORE)))
              STORE=SPCHI(NSEG5,ISEG,1,IPATH)
              TX(MEXTX)                                                 &
     &          =EXP(-(STORE+FRAC5*(SPCLO(NSEG5,ISEG,1,IPATH)-STORE)))

!             SAP AEROSOLS:
              SSAP=SSAPH(ISEG,1,IPATH)
              SSAP=SSAP+FRAC5*(SSAPL(ISEG,1,IPATH)-SSAP)
          ENDIF
          IF(KNTRVL.LE.1)THEN

!             NO CORRELATED-K:
              IF(MODTRN)                                                &
     &          CALL BMOD(1,ISEG,NSEG(1),IPATH,IBINPT,LBMCEN,LBMFLG)

!             COMBINE TRANSMISSIONS IF NOT USING CORRELATED-K METHOD
!               UNIF    UNIFORMLY MIXED GASES TRANSMITTANCE
!               TRACE   TRACE GASES TRANSMITTANCE
              UNIF=TX(36)*TX(44)*TX(46)*TX(47)*TX(50)
              TRACE=TX(52)*TX(54)*TX(55)*TX(56)*TX(11)
              IF(IPATH.EQ.3)THEN

!                 EXCLUDE H2-H2 & H2-HE DIMERS FROM X TRANSMITTANCE:
                  TRAN_X=TX(MEXT+1)
                  DO IEXT=MEXT+2,MEXTX-2
                      TRAN_X=TRAN_X*TX(IEXT)
                  ENDDO

!                 AUXILIARY SPECIES (Y) TRANSMITTANCE:
                  TRAN_Y=1.
                  DO IEXT=MEXTX+1,MEXTX+NMOLY
                      TRAN_Y=TRAN_Y*TX(IEXT)
                  ENDDO
                  IF(TRNOLD.LE.0.)THEN

!                     TRANSMITTANCE RATIO ILL-DEFINED; SET IT TO ZERO.
                      TRNLAY(1)=0.
                      DEPLAY(1)=BIGEXP
                  ELSE
                      DEPLAY(1)=0.
                      TRNLAY(1)=1.
                      TRNNEW=TX(17)*UNIF*TX(31)*TRACE*TRAN_X*TRAN_Y

!                     IF A PATH TRANSMITTANCE THROUGH "ISEG" SEGMENTS IS
!                     EXCEEDINGLY SMALL OR THE GROUND-TO-LEVEL VERTICAL
!                     TRANSMITTANCE IS RELATIVELY SMALL AND BEGUN TO
!                     INCREASE, SET THE TRANSMITTANCE TO EXACTLY ZERO
!                     SO THAT THE MOLECULAR OPTICAL DEPTH OF THE CURRENT
!                     AND SUBSEQUENT LAYERS IS ASSIGNED "BIGEXP".
                      IF(TRNNEW.LT.1.E-30)TRNNEW=0.
                      IF(TRNNEW.LT.TRNOLD)THEN

!                         DETERMINE DECREASE IN TRANSMITTANCE.
                          TRNLAY(1)=TRNNEW/TRNOLD
                          DEPLAY(1)=BIGEXP
                          IF(TRNLAY(1).GT.1.E-30)                       &
     &                      DEPLAY(1)=-LOG(TRNLAY(1))
                          TRNOLD=TRNNEW
                      ENDIF
                  ENDIF

!                 CHECK OZONE CURVE-OF-GROWTH:
                  IF(TX(31).LT.O3TRAN)THEN

!                     ONLY UPDATE O3TRAN IF NEW VALUE IS LOWER:
                      O3TRAN=TX(31)
                  ELSEIF(TX(31)-O3TRAN.GT.O3MAX .AND. LO3TRN)THEN

!                     3-PARAMETER CURTIS-GODSON HAS FAILED FOR OZONE:
                      O3MAX=TX(31)-O3TRAN
                      WRITE(IPR,'(A,F9.2,2(A,/21X),3(A,F10.7),/21X,A)') &
     &                  ' Warning from LOOP:  At spectral frequency',   &
     &                  VCEN,' CM-1, the 3-parameter Curtis-Godson',    &
     &                  'formulation is not producing a monotonic '//   &
     &                  'curve-of-growth.','The O3 transmittance '//    &
     &                  'increased',O3MAX,' from',O3TRAN,' to',TX(31),  &
     &                  'Run the Correlated-k algorithm to avoid this'//&
     &                  ' problem.  This warning will not be repeated.'
                      LO3TRN=.FALSE.
                  ENDIF
                  DEPLAY(1)=DEPLAY(1)+TX(14)
                  TX9LAY=EXP(-TX(14))
                  TX9CUM=TX9CUM*TX9LAY
                  TX(9)=TX9CUM*TRNNEW
                  TRNLAY(1)=TRNLAY(1)*TX9LAY
                  IF(IEMSCT.EQ.2)THEN
                      CALL SSRAD(IPH,1,ISEG,IPATH,VCEN,S0,PFMOL0,       &
     &                  PFMOL2,TRNNEW,SSAP,TX(1),TSNOBS,TSNREF,SUMSSS)
                      IF(AMOD3D.EQ.'M')THEN
                          H2OCUM=H2OCUM*TX(5)
                          UMXCUM=UMXCUM*TX(3)*TX(64)*TX(4)*TX(65)*T_CH4
                          O3CUM=O3CUM*TX(8)
                          TRAN_X=TRAN_X*TX(MEXTX-1)*TX(MEXTX)
                          CALL MCMOL(VCEN,ISEG,NSEG(1),NMWAVE,TX(17)    &
     &                      *H2OCUM,UNIF*TRACE*TRAN_X*UMXCUM,TX(31)     &
     &                      *O3CUM,TRAN_Y,TL_H2O,TL_UMX,TL_O3,TL_AUX)
                      ENDIF
                  ELSEIF(AMOD3D.EQ.'C')THEN
                      CALL MCCONT(VCEN,ISEG,NSEG(1),NMWAVE)
                  ENDIF
              ELSEIF(IEMSCT.EQ.2)THEN

!                 INCLUDE H2-H2 & H2-HE DIMERS IN XY TRANSMITTANCE:
                  TRAN_X=TX(MEXT+1)
                  DO IEXT=MEXT+2,MEXTX
                      TRAN_X=TRAN_X*TX(IEXT)
                  ENDDO

!                 AUXILIARY SPECIES (Y) TRANSMITTANCE:
                  TL_AUX=1.
                  DO IEXT=MEXTX+1,MEXTX+NMOLY
                      TL_AUX=TL_AUX*TX(IEXT)
                  ENDDO

!                 CALCULATE TMOL, THE MOLECULAR TRANSMITTANCES
!                 TO THE SUN.  VALUES OF TMOL LESS THAN 1.E-8
!                 ARE NOT ACCURATE AND ARE SET TO ZERO.
                  TMOL=TX(17)*UNIF*TX(31)*TRACE*TRAN_X*TL_AUX
                  IF(TMOL.LT.1.E-8)TMOL=0.
                  CALL SSRAD(IPH,1,ISEG,IPATH,VCEN,S0,PFMOL0,PFMOL2,    &
     &              TMOL,SSAP,TX(1),TSNOBS,TSNREF,SUMSSS)
                  IF(AMOD3D.NE.'M')GOTO(10,20),IPATH

!                 COMPUTE L-PATH TRANSMITTANCES FOR MCSCENE DATABASES.
                  TL_H2O=TX(17)*TX(5)
                  TL_UMX=UNIF*TRACE*TRAN_X                              &
     &              *TX(3)*TX(64)*TX(4)*TX(65)*T_CH4
                  TL_O3=TX(31)*TX(8)
                  IF(IPATH.EQ.2)GOTO 20

!                 PROCESS MCSCENE DATA FOR DIRECT SENSOR TO SUN PATH:
                  CALL MCMOL(VCEN,0,NSEG(1),NMWAVE,                     &
     &              1.,1.,1.,1.,TL_H2O,TL_UMX,TL_O3,TL_AUX)
                  GOTO 10
              ENDIF
          ELSEIF(IPATH.LT.3)THEN

!             CK SOLAR PATH:
              CALL CKLOSS(JNTRVL,ISEG+IPATH-1,1)
              CALL SSCORK(.FALSE.,IPH,1,ISEG,IPATH,JNTRVL,VCEN,S0,      &
     &          PFMOL0,PFMOL2,SSAP,TX(1),TSNOBS,SUMSSS,RADCUM)
              GOTO(10,20),IPATH
          ELSE

!             CK OPTICAL (LINE-OF-SIGHT) PATH:
              IF(AMOD3D.EQ.'T')THEN

!                 EXCLUDE AEROSOL, CLOUD & RAIN FROM TX(14) FOR MOD3D:
                  TX14SV=TX(14)
                  TX(14)=TX(80)
                  CALL CKLOS(JNTRVL,ISEG,1)
                  TX(14)=TX14SV
              CALL M3D(JNTRVL,ISEG,NSEG(1),SUBINT,WPTH(0,ISEG,1),LABEL)   
              ELSE

!                 LINE-OF-SIGHT PATH WITHOUT MOD3D:
                  CALL CKLOS(JNTRVL,ISEG,1)
              ENDIF
              IF(IEMSCT.EQ.2)                                           &
     &          CALL SSCORK(.FALSE.,IPH,1,ISEG,IPATH,JNTRVL,VCEN,       &
     &          S0,PFMOL0,PFMOL2,SSAP,TX(1),TSNOBS,SUMSSS,RADCUM)
          ENDIF

!         INITIALIZE TOTAL TRANSMITTANCE:
          TX(9)=0.

!         DEFINE BLACKBODY FUNCTIONS.
          BLAYER=BBFN(TSEG(ISEG,1),VCEN)
          BBOUND=BBFN(PATH_T(ISEG-1,1),VCEN)

!         LOOP OVER CORRELATED-K METHOD SUB-INTERVALS:
          SMTRNL=0.D0
          SMTHML=0.D0
          DO INTRVL=1,JNTRVL
              OPTDEP=DEPLAY(INTRVL)
              IF(OPTDEP.LT..02)THEN
                  THMLAY=OPTDEP*(BLAYER+OPTDEP*(BBOUND-4*BLAYER)/6)
              ELSE
                  EMSLAY=1.-TRNLAY(INTRVL)
                  THMLAY=EMSLAY*BBOUND+2*(BLAYER-BBOUND)                &
     &              *(EMSLAY/OPTDEP-TRNLAY(INTRVL))
              ENDIF
              THMLAY=TRNCUM(INTRVL)*THMLAY
              SMTHML=SMTHML+DBLE(SUBINT(INTRVL)*THMLAY)
              RADCUM(INTRVL)=RADCUM(INTRVL)+THMLAY
              TRNCUM(INTRVL)=TRNCUM(INTRVL)*TRNLAY(INTRVL)
              TX(9)=TX(9)+SUBINT(INTRVL)*TRNCUM(INTRVL)
              SMTRNL=SMTRNL+DBLE(SUBINT(INTRVL)*TRNLAY(INTRVL))
          ENDDO
          THMCUM=THMCUM+REAL(SMTHML)

!         IF NOPRNT=-1, PRINT WEIGHTING FUNCTION DATA.
          IF(NOPRNT.LE.-1 .AND. .NOT.NO_WRT)THEN
              IF(BINOUT)THEN
                  CALL BNWT9(VCEN,PTHALT(ISEG-1,1),PTHALT(ISEG,1),      &
     &              BLAYER,BBOUND,TX(9),SMTRNL,SMTHML,THMCUM)
              ELSE
                  WRITE(IPR1,'(F8.2,2F10.5,3X,2(1P,2E11.3,0P,2F11.7))') &
     &              VCEN,PTHALT(ISEG-1,1),PTHALT(ISEG,1),BLAYER,        &
     &              BBOUND,TX(9),SMTRNL,SMTHML,THMCUM
              ENDIF
          ENDIF

!         K-DISTRIBUTION DEPENDENT OUTPUT:
          IF(KPRINT .AND. (ISEG.LE.IKOUT(1) .OR. ISEG.EQ.NSEG(1)))THEN
              IF(BINOUT)THEN
                  NTKDIS=NTKDIS+1
                  WRITE(JTKDIS,REC=NTKDIS)PTHALT(ISEG,1),PTHRNG(ISEG,1),&
     &              TX(9),(TRNCUM(INTRVL),INTRVL=1,JNTRVL)
                  NRKDIS=NRKDIS+1
                  WRITE(JRKDIS,REC=NRKDIS)PTHALT(ISEG,1),PTHRNG(ISEG,1),&
     &              THMCUM+SUMSSS,(RADCUM(INTRVL),INTRVL=1,JNTRVL)
              ELSE
                  WRITE(ITKDIS,TFRMT)PTHALT(ISEG,1),PTHRNG(ISEG,1),     &
     &              TX(9),(TRNCUM(INTRVL),INTRVL=1,JNTRVL)
                  WRITE(IRKDIS,RFRMT)PTHALT(ISEG,1),PTHRNG(ISEG,1),     &
     &              THMCUM+SUMSSS,(RADCUM(INTRVL),INTRVL=1,JNTRVL)
              ENDIF
          ENDIF

!         IF TOTAL TRANSMISSION HAS DROPPED TO ZERO AND LVBND5 IS
!         FALSE, EXIT LAYER LOOP UNLESS THE CORRELATED-K APPROACH
!         IS BEING USED AND WEIGHTING FUNCTIONS ARE BEING PRINTED.
          IF(TX(9).LE.0. .AND. .NOT.LVBND5 .AND. AMOD3D.EQ.' '          &
     &      .AND. (KNTRVL.EQ.1 .OR. NOPRNT.LE.-1))GOTO 40
      ENDDO

!     LAYER LOOP EXIT:
   40 CONTINUE

!     LINE-OF-SIGHT PATH:
      IF(IEMSCT.EQ.2 .AND. KNTRVL.GT.1)THEN
          ISEG=NSEG(1)+1
          S0TRN(1,ISEG)=S0TRN(1,ISEG)*TRNCUM(1)
          TSNREF=SUBINT(1)*S0TRN(1,ISEG)

!         IF THE GROUND SCATTERING POINT IS IN SHADOW, THE GROUND
!         REFLECTED SOLAR IRRADIANCE, TSNREF, IS ZERO.
          IF(TSNREF.GT.0.)THEN
              DO INTRVL=2,JNTRVL
                  S0TRN(INTRVL,ISEG)=S0TRN(INTRVL,ISEG)*TRNCUM(INTRVL)
                  TSNREF=TSNREF+SUBINT(INTRVL)*S0TRN(INTRVL,ISEG)
              ENDDO
          ELSE
              DO INTRVL=2,JNTRVL
                  S0TRN(INTRVL,ISEG)=0.
              ENDDO
          ENDIF
      ENDIF

!     NO FLUXES SINCE MULTIPLE SCATTERING NOT DONE.
      IF(.NOT.KPRINT .OR. NO_WRT)RETURN

!     OUTPUT SURFACE TERMS TO *.r_k FILE:
      IF(TPTEMP(1).GT.0. .AND. TPHDIR(1).LT.1.)THEN

!         THERMAL EMISSION AT BOUNDARY:
          GREY=BBFN(TPTEMP(1),VCEN)*(1.-TPHDIR(1))
      ELSE

!         NO THERMAL EMISSION AT BOUNDARY:
          GREY=0.
      ENDIF
      IF(BINOUT)THEN

!         UPDATE K-DEPENDENT DATA BINARY FILE HEADER:
          NRKDIS=NRKDIS+1
          WRITE(JTKDIS,REC=NTHEAD)                                      &
     &      N4BYTE,NTREC,NRANGE,NTKDIS,WDKDIS,V0KDIS,.TRUE.
          WRITE(JRKDIS,REC=NRHEAD)                                      &
     &      N4BYTE,NRREC,NRANGE,NRKDIS,WDKDIS,V0KDIS,.TRUE.
          IF(TSNREF.GT.0. .AND. GROUND)THEN

!             SURFACE EMISSION WITH DIRECT SOLAR REFLECTION.
              SURFAC=TPBRDF(1)*CSNSRF(1)/PI
              WRITE(JRKDIS,REC=NRKDIS)PTHALT(NSEG(1),1),                &
     &          PTHRNG(NSEG(1),1),TX(9)*GREY+SURFAC*TSNREF,             &
     &          (GREY*TRNCUM(INTRVL)+SURFAC*S0TRN(INTRVL,ISEG),         &
     &          INTRVL=1,JNTRVL)
          ELSE

!             SURFACE EMISSION WITH NO REFLECTION.
              WRITE(JRKDIS,REC=NRKDIS)PTHALT(NSEG(1),1),                &
     &          PTHRNG(NSEG(1),1),TX(9)*GREY,                           &
     &          (GREY*TRNCUM(INTRVL),INTRVL=1,JNTRVL)
          ENDIF
      ELSE
          RFRMT(2:2)='/'
          IF(TSNREF.GT.0. .AND. GROUND)THEN

!             SURFACE EMISSION WITH DIRECT SOLAR REFLECTION.
              SURFAC=TPBRDF(1)*CSNSRF(1)/PI
              WRITE(IRKDIS,RFRMT)PTHALT(NSEG(1),1),PTHRNG(NSEG(1),1),   &
     &          GREY*TX(9)+SURFAC*TSNREF,(GREY*TRNCUM(INTRVL)           &
     &          +SURFAC*S0TRN(INTRVL,ISEG),INTRVL=1,JNTRVL)
          ELSE

!             SURFACE EMISSION WITH NO REFLECTION.
              WRITE(IRKDIS,RFRMT)PTHALT(NSEG(1),1),PTHRNG(NSEG(1),1),   &
     &          GREY*TX(9),(GREY*TRNCUM(INTRVL),INTRVL=1,JNTRVL)
          ENDIF
          RFRMT(2:2)=' '
      ENDIF
      RETURN
      END
