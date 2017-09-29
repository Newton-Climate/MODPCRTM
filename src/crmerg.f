      SUBROUTINE CRMERG

!     THIS ROUTINE MERGES TOGETHER CLOUD/RAIN AND OLD ATMOSPHERIC
!     PROFILES.  WATER DROPLET DENSITIES [GM/M3], ICE PARTICLE
!     DENSITIES [GM/M3] AND RAIN RATES [MM/HR] ARE STORED IN
!     DENSTY(66,.), DENSTY(67,.) AND DENSTY(3,.), RESPECTIVELY.
!     LAYER BOUNDARIES ARE MERGED TOGETHER IF THEY DIFFER BY LESS
!     THAN "ZTOL" KM (HALF A METER).
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'YPROP.h'
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'

!     /M_PTWO/
!       PPROF    PRESSURE PROFILE [MB].
!       TPROF    TEMPERATURE PROFILE [K].
!       WH2O     H2O VOLUME MIXING RATIO PROFILE [PPMV].
!       WO3      O3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL PPROF,TPROF,WH2O,WO3
      COMMON/M_PTWO/PPROF(LAYDIM),TPROF(LAYDIM),WH2O(LAYDIM),WO3(LAYDIM)

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

!     /MDATXY/
!       WMOLXT   CROSS-SECTION MOLECULE DENSITY PROFILE [PPMV].
!       WMOLYT   AUXILIARY (Y) MOLECULE DENSITY PROFILE [PPMV].
      REAL WMOLXT,WMOLYT
      COMMON/MDATXY/WMOLXT(NMOLX,LAYDIM),WMOLYT(MMOLY,LAYDIM)

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

!     /CLDRN/
      DOUBLE PRECISION ZCLDRN(NZCLD)
      INTEGER NBND
      REAL DRPWAT(NZCLD),PRTICE(NZCLD),RNPROF(NZCLD)
      COMMON/CLDRN/ZCLDRN,NBND,DRPWAT,PRTICE,RNPROF

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      DOUBLE PRECISION CTHIK,CALT
      REAL CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,        &
     &  CCOLIP,CHUMID,ASYMWD,ASYMIP

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

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       MSAP    MAXIMUN SAP LEVEL AFTER OFFSETING.
      INTEGER NEWBND,MBND,IBND,I,MLOLD,IWVSAP,ILGSAP,IANSAP,MSAP
      LOGICAL LSET,WARNWD,WARNIP
      REAL FAC,TRATIO,CLDHUM
      DOUBLE PRECISION DPFAC

!     FUNCTIONS:
      REAL EXPINT

!     DATA:
!       TFREEZ   TEMPERATURE BELOW WHICH A WARNING IS GENERATED IF
!                LIQUID WATER DROPLET DENSITY IS POSITIVE [K]
!       TMELT    TEMPERATURE ABOVE WHICH A WARNING IS GENERATED IF
!                ICE PARTICLE DENSITY IS POSITIVE [K]
      REAL TFREEZ,TMELT
      DATA TFREEZ/260./,TMELT/278./

!     CHECK THAT CLOUD/RAIN ALTITUDES ARE BOUNDED BY ZM ALTITUDES
      IF(ZCLDRN(1).LT.ZM(1) .OR. ZCLDRN(NBND).GT.ZM(ML))THEN
          WRITE(IPR,'(/A,2(F10.5,A),/18X,A,2(F10.5,A))')                &
     &      ' Error in CRMERG:  Cloud/Rain model bounding altitudes (', &
     &      ZCLDRN(1),' and',ZCLDRN(NBND),' km) are',                   &
     &      ' not bracketed by atmosphere bounding altitudes (',        &
     &      ZM(1),' AND', ZM(ML),' km).'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP
      ENDIF

!     CHECK THAT RAIN RATE AND CLOUD WATER DROPLET AND ICE PARTICLE
!     DENSITIES ARE INDEED ZERO AT THE TOP OF THE PROFILE.
      IF(DRPWAT(NBND).NE.0. .OR. PRTICE(NBND).NE.0.                     &
     &                      .OR. RNPROF(NBND).NE.0.)THEN
          WRITE(IPR,'(/A,/18X,A)')' Error in CRMERG:  Rain rate and'    &
     &      //' cloud water droplet and ice particle densities and',    &
     &      ' are not zero at the top of the cloud/rain profile.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP
      ENDIF

!     INITIALIZE WARNING MESSAGE VARIABLES
!     WARNWD   A LOGICAL VARIABLE THAT IS SET TO FALSE AFTER THE USER
!              HAS BEEN WARNED THAT WATER DROPLETS EXIST BELOW TFREEZ.
!     WARNIP   A LOGICAL VARIABLE THAT IS SET TO FALSE AFTER THE USER
!              HAS BEEN WARNED THAT ICE PARTICLES EXIST ABOVE TMELT.
      WARNWD=.TRUE.
      WARNIP=.TRUE.

!     DETERMINE THE RELATIVE HUMIDITY [%] WITHIN THE CLOUD.
!     THE RELATIVE HUMIDITY MUST BE POSITIVE (A DRY ATMOSPHERE
!     IS NOT ALLOWED), AND SUPER-SATURATION IS LIMITED TO 5%.
      IF(CHUMID.GT.0.)THEN
          CLDHUM=CHUMID
      ELSE
          CLDHUM=100.
      ENDIF

!     TEMPORARILY TRANSLATE BOUNDARY LAYER DATA TO MAXIMUM
!     LAYER INDICES TO MAKE SPACE FOR ADDITIONAL LAYERS.
      NEWBND=LAYDM1
      DO MBND=ML,2,-1
          NEWBND=NEWBND-1
          ZM(NEWBND)=ZM(MBND)
          PPROF(NEWBND)=PPROF(MBND)
          TPROF(NEWBND)=TPROF(MBND)
          RELHUM(NEWBND)=RELHUM(MBND)
          WH2O(NEWBND)=WH2O(MBND)
          LRHSET(NEWBND)=LRHSET(MBND)
          WCO2(NEWBND)=WCO2(MBND)
          WO3(NEWBND)=WO3(MBND)
          WN2O(NEWBND)=WN2O(MBND)
          WCO(NEWBND)=WCO(MBND)
          WCH4(NEWBND)=WCH4(MBND)
          WO2(NEWBND)=WO2(MBND)
          WHNO3(NEWBND)=WHNO3(MBND)
          WNO(NEWBND)=WNO(MBND)
          WSO2(NEWBND)=WSO2(MBND)
          WNO2(NEWBND)=WNO2(MBND)
          WNH3(NEWBND)=WNH3(MBND)
          DO I=1,NMOLX
              WMOLXT(I,NEWBND)=WMOLXT(I,MBND)
          ENDDO
          DO I=1,NMOLY
              WMOLYT(I,NEWBND)=WMOLYT(I,MBND)
          ENDDO
          DENSTY(7,NEWBND)=DENSTY(7,MBND)
          DENSTY(12,NEWBND)=DENSTY(12,MBND)
          DENSTY(13,NEWBND)=DENSTY(13,MBND)
          DENSTY(14,NEWBND)=DENSTY(14,MBND)
          DENSTY(15,NEWBND)=DENSTY(15,MBND)
          DENSTY(16,NEWBND)=DENSTY(16,MBND)

!         SAP OPTION:
          IF(LSAP)THEN
              IF(MBND.LE.LEVSAP)THEN
                  DO IWVSAP=1,NWVSAP
                      BSAPX(IWVSAP,NEWBND)=BSAPX(IWVSAP,MBND)
                      BSAPA(IWVSAP,NEWBND)=BSAPA(IWVSAP,MBND)
                      DO ILGSAP=1,NLGSAP
                          LEGSAP(ILGSAP,IWVSAP,NEWBND)                  &
     &                      =LEGSAP(ILGSAP,IWVSAP,MBND)
                      ENDDO
                      DO IANSAP=1,NANSAP
                          PFSAP(IANSAP,IWVSAP,NEWBND)                   &
     &                      =PFSAP(IANSAP,IWVSAP,MBND)
                      ENDDO
                  ENDDO
              ENDIF
          ENDIF
      ENDDO

!     INITIALIZE RAIN RATE AND WATER DROPLET AND
!     ICE PARTICLE DENSITIES AT THE GROUND.
      ML=1
      DENSTY(66,1)=0.
      DENSTY(67,1)=0.
      DENSTY(3,1)=0.

!     LOOP OVER OLD ATMOSPHERIC LAYERS
      IBND=1
      IF(LSAP)MSAP=LEVSAP+NEWBND-1
      DO MBND=NEWBND,LAYDIM
          LSET=.FALSE.
   10     CONTINUE
          IF(ZCLDRN(IBND).LE.ZM(MBND))THEN

!             CLOUD LAYER BOUNDARY WITHIN OLD ATMOSPHERIC LAYER
              IF(ZCLDRN(IBND).LT.ZM(ML)+ZTOL)THEN

!                 MERGE INTO LOWER LAYER BOUNDARY:
                  DENSTY(66,ML)=DRPWAT(IBND)
                  DENSTY(67,ML)=PRTICE(IBND)
                  DENSTY(3,ML)=RNPROF(IBND)
                  IF(CLDHUM.LE.105. .AND. (DRPWAT(IBND).GT.0. .OR.      &
     &              PRTICE(IBND).GT.0. .OR. RNPROF(IBND).GT.0.))THEN
                      RELHUM(ML)=CLDHUM
                      TRATIO=273.15/TPROF(ML)
                      WH2O(ML)=.01*CLDHUM*TRATIO*                       &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.TRUE.
                      IF(WARNWD .AND. DRPWAT(IBND).GT.0.                &
     &                          .AND. TPROF(ML).LT.TFREEZ)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE LIQUID WATER DROPLET DENSITY',    &
     &                      ' IS POSITIVE [',DRPWAT(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      TPROF(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNWD=.FALSE.
                      ELSEIF(WARNIP .AND. PRTICE(IBND).GT.0.            &
     &                              .AND. TPROF(ML).GT.TMELT)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE ICE PARITCLE DENSITY',            &
     &                      ' IS POSITIVE [',PRTICE(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      TPROF(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNIP=.FALSE.
                      ENDIF
                  ENDIF
              ELSEIF(ZCLDRN(IBND).GT.ZM(MBND)-ZTOL)THEN

!                 MERGE INTO UPPER LAYER BOUNDARY
                  LSET=.TRUE.
                  DENSTY(66,MBND)=DRPWAT(IBND)
                  DENSTY(67,MBND)=PRTICE(IBND)
                  DENSTY(3,MBND)=RNPROF(IBND)
                  IF(CLDHUM.LE.105. .AND. (DRPWAT(IBND).GT.0. .OR.      &
     &              PRTICE(IBND).GT.0. .OR. RNPROF(IBND).GT.0.))THEN
                      RELHUM(MBND)=CLDHUM
                      TRATIO=273.15/TPROF(MBND)
                      WH2O(MBND)=.01*CLDHUM*TRATIO*                     &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.TRUE.
                  ENDIF
              ELSEIF(ML+1.LT.MBND)THEN

!                 ADD NEW ZM LAYER BOUNDARY
                  MLOLD=ML
                  ML=ML+1
                  ZM(ML)=ZCLDRN(IBND)
                  DENSTY(66,ML)=DRPWAT(IBND)
                  DENSTY(67,ML)=PRTICE(IBND)
                  DENSTY(3,ML)=RNPROF(IBND)
                  DPFAC=(ZM(ML)-ZM(MLOLD))/(ZM(MBND)-ZM(MLOLD))
                  FAC=SNGL(DPFAC)
                  PPROF(ML)=EXPINT(PPROF(MLOLD),PPROF(MBND),FAC)
                  TPROF(ML)=EXPINT(TPROF(MLOLD),TPROF(MBND),FAC)
                  TRATIO=273.15/TPROF(ML)
                  IF(CLDHUM.LE.105. .AND. (DRPWAT(IBND).GT.0. .OR.      &
     &              PRTICE(IBND).GT.0. .OR. RNPROF(IBND).GT.0.))THEN
                      RELHUM(ML)=CLDHUM
                      WH2O(ML)=.01*CLDHUM*TRATIO*                       &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.TRUE.
                      IF(WARNWD .AND. DRPWAT(IBND).GT.0.                &
     &                  .AND. TPROF(ML).LT.TFREEZ)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE LIQUID WATER DROPLET DENSITY',    &
     &                      ' IS POSITIVE [',DRPWAT(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      TPROF(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNWD=.FALSE.
                      ELSEIF(WARNIP .AND. PRTICE(IBND).GT.0.            &
     &                              .AND. TPROF(ML).GT.TMELT)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE ICE PARITCLE DENSITY',            &
     &                      ' IS POSITIVE [',PRTICE(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      TPROF(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNIP=.FALSE.
                      ENDIF
                  ELSE
                      RELHUM(ML)=EXPINT(RELHUM(MLOLD),RELHUM(MBND),FAC)
                      WH2O(ML)=.01*RELHUM(ML)*TRATIO*                   &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.FALSE.
                  ENDIF
                  WCO2(ML)=EXPINT(WCO2(MLOLD),WCO2(MBND),FAC)
                  WO3(ML)=EXPINT(WO3(MLOLD),WO3(MBND),FAC)
                  WN2O(ML)=EXPINT(WN2O(MLOLD),WN2O(MBND),FAC)
                  WCO(ML)=EXPINT(WCO(MLOLD),WCO(MBND),FAC)
                  WCH4(ML)=EXPINT(WCH4(MLOLD),WCH4(MBND),FAC)
                  WO2(ML)=EXPINT(WO2(MLOLD),WO2(MBND),FAC)
                  WHNO3(ML)=EXPINT(WHNO3(MLOLD),WHNO3(MBND),FAC)
                  WNO(ML)=EXPINT(WNO(MLOLD),WNO(MBND),FAC)
                  WSO2(ML)=EXPINT(WSO2(MLOLD),WSO2(MBND),FAC)
                  WNO2(ML)=EXPINT(WNO2(MLOLD),WNO2(MBND),FAC)
                  WNH3(ML)=EXPINT(WNH3(MLOLD),WNH3(MBND),FAC)
                  DO I=1,NMOLX
                      WMOLXT(I,ML)=                                     &
     &                  EXPINT(WMOLXT(I,MLOLD),WMOLXT(I,MBND),FAC)
                  ENDDO
                  DO I=1,NMOLY
                      WMOLYT(I,ML)=                                     &
     &                  EXPINT(WMOLYT(I,MLOLD),WMOLYT(I,MBND),FAC)
                  ENDDO
                  DENSTY(7,ML)=                                         &
     &              EXPINT(DENSTY(7,MLOLD),DENSTY(7,MBND),FAC)
                  DENSTY(12,ML)=                                        &
     &              EXPINT(DENSTY(12,MLOLD),DENSTY(12,MBND),FAC)
                  DENSTY(13,ML)=                                        &
     &              EXPINT(DENSTY(13,MLOLD),DENSTY(13,MBND),FAC)
                  DENSTY(14,ML)=                                        &
     &              EXPINT(DENSTY(14,MLOLD),DENSTY(14,MBND),FAC)
                  DENSTY(15,ML)=                                        &
     &              EXPINT(DENSTY(15,MLOLD),DENSTY(15,MBND),FAC)
                  DENSTY(16,ML)=                                        &
     &              EXPINT(DENSTY(16,MLOLD),DENSTY(16,MBND),FAC)

!                 SAP OPTION:
                  IF(LSAP)THEN
                      IF(MBND.LE.MSAP)THEN
                          LEVSAP=LEVSAP+1
                          DO IWVSAP=1,NWVSAP
                              BSAPX(IWVSAP,ML)                          &
     &                          =EXPINT(BSAPX(IWVSAP,MLOLD),            &
     &                                  BSAPX(IWVSAP,MBND),FAC)
                              BSAPA(IWVSAP,ML)                          &
     &                          =EXPINT(BSAPA(IWVSAP,MLOLD),            &
     &                                  BSAPA(IWVSAP,MBND),FAC)
                              DO ILGSAP=1,NLGSAP
                                  LEGSAP(ILGSAP,IWVSAP,ML)              &
     &                              =LEGSAP(ILGSAP,IWVSAP,MLOLD)        &
     &                              +FAC*(LEGSAP(ILGSAP,IWVSAP,MBND)    &
     &                                   -LEGSAP(ILGSAP,IWVSAP,MLOLD))
                              ENDDO
                              DO IANSAP=1,NANSAP
                                  PFSAP(IANSAP,IWVSAP,ML)               &
     &                              =PFSAP(IANSAP,IWVSAP,MLOLD)         &
     &                              +DPFAC*(PFSAP(IANSAP,IWVSAP,MBND)   &
     &                                     -PFSAP(IANSAP,IWVSAP,MLOLD))
                              ENDDO
                          ENDDO
                      ENDIF
                  ENDIF
              ELSE

!                 NO MORE SPACE IN ARRAYS FOR AN ADDITIONAL LAYER
                  WRITE(IPR,'(/3A,/14X,A,2(A,I4))')' FATAL ERROR: ',    &
     &              ' FILE "PARAMS.h" PARAMETER "LAYDIM" MUST',         &
     &              ' BE INCREASED.',' IT SUFFICES TO INCREASE',        &
     &              ' LAYDIM FROM',LAYDIM,' TO',LAYDIM+NBND-IBND+1
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP
              ENDIF

!             EXIT LOOP IF ALL CLOUD BOUNDARIES HAVE BEEN INTEGRATED.
              IF(IBND.GE.NBND)GOTO 20

!             INCREMENT CLOUD BOUNDARY INDEX AND START AGAIN
              IBND=IBND+1
              GOTO 10
          ENDIF

!         LINEARLY INTERPOLATE WATER PARTICLE DENSITIES AND RAIN RATE.
          IF(.NOT.LSET)THEN
              FAC=SNGL((ZM(MBND)-ZM(ML))/(ZCLDRN(IBND)-ZM(ML)))
              !FAC=1.
              DENSTY(66,MBND)=DENSTY(66,ML)                             &
     &          +FAC*(DRPWAT(IBND)-DENSTY(66,ML))
              DENSTY(67,MBND)=DENSTY(67,ML)                             &
     &          +FAC*(PRTICE(IBND)-DENSTY(67,ML))
              DENSTY(3,MBND)=DENSTY(3,ML)                               &
     &          +FAC*(RNPROF(IBND)-DENSTY(3,ML))
              IF(CLDHUM.LE.105. .AND. (DENSTY(66,MBND).GT.0. .OR.       &
     &          DENSTY(67,MBND).GT.0. .OR. DENSTY( 3,MBND).GT.0.))THEN
                  RELHUM(MBND)=CLDHUM
                  TRATIO=273.15/TPROF(MBND)
                  WH2O(MBND)=.01*CLDHUM*TRATIO*                         &
     &              EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                  LRHSET(ML)=.TRUE.
              ENDIF
          ENDIF

!         TRANSLATE LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX.
          ML=ML+1
          ZM(ML)=ZM(MBND)
          DENSTY(66,ML)=DENSTY(66,MBND)
          DENSTY(67,ML)=DENSTY(67,MBND)
          DENSTY(3,ML)=DENSTY(3,MBND)
          PPROF(ML)=PPROF(MBND)
          TPROF(ML)=TPROF(MBND)
          IF(WARNWD .AND. DENSTY(66,ML).GT.0.                           &
     &              .AND. TPROF(ML).LT.TFREEZ)THEN
              WRITE(IPR,'(//A,F9.4,A,/(11X,A,F9.4,A))')                 &
     &          ' WARNING:  AT ALTITUDE',ZM(ML),                        &
     &          ' KM, THE LIQUID WATER DROPLET DENSITY',                &
     &          'IS POSITIVE [',DENSTY(66,ML),' GM/M3] EVEN THOUGH THE',&
     &          'TEMPERATURE IS',TPROF(ML)-273.15,' DEGREES CELSIUS.'
              WARNWD=.FALSE.
          ENDIF
          IF(WARNIP .AND. DENSTY(67,ML).GT.0.                           &
     &              .AND. TPROF(ML).GT.TMELT)THEN
              WRITE(IPR,'(//A,F9.4,A,/(11X,A,F9.4,A))')' WARNING: '//   &
     &          ' AT ALTITUDE',ZM(ML),' KM, THE ICE PARITCLE DENSITY',  &
     &          'IS POSITIVE [',DENSTY(67,ML),' GM/M3] EVEN THOUGH THE',&
     &          'TEMPERATURE IS',TPROF(ML)-273.15,' DEGREES CELSIUS.'
              WARNIP=.FALSE.
          ENDIF
          RELHUM(ML)=RELHUM(MBND)
          WH2O(ML)=WH2O(MBND)
          LRHSET(ML)=LRHSET(MBND)
          WCO2(ML)=WCO2(MBND)
          WO3(ML)=WO3(MBND)
          WN2O(ML)=WN2O(MBND)
          WCO(ML)=WCO(MBND)
          WCH4(ML)=WCH4(MBND)
          WO2(ML)=WO2(MBND)
          WHNO3(ML)=WHNO3(MBND)
          WNO(ML)=WNO(MBND)
          WSO2(ML)=WSO2(MBND)
          WNO2(ML)=WNO2(MBND)
          WNH3(ML)=WNH3(MBND)
          DO I=1,NMOLX
              WMOLXT(I,ML)=WMOLXT(I,MBND)
          ENDDO
          DO I=1,NMOLY
             WMOLYT(I,ML)=WMOLYT(I,MBND)
          ENDDO
          DENSTY(7,ML)=DENSTY(7,MBND)
          DENSTY(12,ML)=DENSTY(12,MBND)
          DENSTY(13,ML)=DENSTY(13,MBND)
          DENSTY(14,ML)=DENSTY(14,MBND)
          DENSTY(15,ML)=DENSTY(15,MBND)
          DENSTY(16,ML)=DENSTY(16,MBND)

!         SAP OPTION:
          IF(LSAP)THEN
              IF(MBND.LE.MSAP)THEN
                  DO IWVSAP=1,NWVSAP
                      BSAPX(IWVSAP,ML)=BSAPX(IWVSAP,MBND)
                      BSAPA(IWVSAP,ML)=BSAPA(IWVSAP,MBND)
                      DO ILGSAP=1,NLGSAP
                          LEGSAP(ILGSAP,IWVSAP,ML)                      &
     &                      =LEGSAP(ILGSAP,IWVSAP,MBND)
                      ENDDO
                      DO IANSAP=1,NANSAP
                          PFSAP(IANSAP,IWVSAP,ML)                       &
     &                      =PFSAP(IANSAP,IWVSAP,MBND)
                      ENDDO
                  ENDDO
              ENDIF
          ENDIF
      ENDDO
   20 CONTINUE

!     TRANSLATE REMAINING LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX
!     SETTING CLOUD PARTICLE DENSITIES AND RAIN RATES TO ZERO.
      NEWBND=MBND
      DO MBND=NEWBND,LAYDIM
          ML=ML+1
          ZM(ML)=ZM(MBND)
          DENSTY(66,ML)=0.
          DENSTY(67,ML)=0.
          DENSTY(3,ML)=0.
          PPROF(ML)=PPROF(MBND)
          TPROF(ML)=TPROF(MBND)
          RELHUM(ML)=RELHUM(MBND)
          WH2O(ML)=WH2O(MBND)
          LRHSET(ML)=LRHSET(MBND)
          WCO2(ML)=WCO2(MBND)
          WO3(ML)=WO3(MBND)
          WN2O(ML)=WN2O(MBND)
          WCO(ML)=WCO(MBND)
          WCH4(ML)=WCH4(MBND)
          WO2(ML)=WO2(MBND)
          WHNO3(ML)=WHNO3(MBND)
          WNO(ML)=WNO(MBND)
          WSO2(ML)=WSO2(MBND)
          WNO2(ML)=WNO2(MBND)
          WNH3(ML)=WNH3(MBND)
          DO I=1,NMOLX
              WMOLXT(I,ML)=WMOLXT(I,MBND)
          ENDDO
          DO I=1,NMOLY
              WMOLYT(I,ML)=WMOLYT(I,MBND)
          ENDDO
          DENSTY(7,ML)=DENSTY(7,MBND)
          DENSTY(12,ML)=DENSTY(12,MBND)
          DENSTY(13,ML)=DENSTY(13,MBND)
          DENSTY(14,ML)=DENSTY(14,MBND)
          DENSTY(15,ML)=DENSTY(15,MBND)
          DENSTY(16,ML)=DENSTY(16,MBND)

!         SAP OPTION:
          IF(LSAP)THEN
              IF(MBND.LE.MSAP)THEN
                  DO IWVSAP=1,NWVSAP
                      BSAPX(IWVSAP,ML)=BSAPX(IWVSAP,MBND)
                      BSAPA(IWVSAP,ML)=BSAPA(IWVSAP,MBND)
                      DO ILGSAP=1,NLGSAP
                          LEGSAP(ILGSAP,IWVSAP,ML)                      &
     &                      =LEGSAP(ILGSAP,IWVSAP,MBND)
                      ENDDO
                      DO IANSAP=1,NANSAP
                          PFSAP(IANSAP,IWVSAP,ML)                       &
     &                      =PFSAP(IANSAP,IWVSAP,MBND)
                      ENDDO
                  ENDDO
              ENDIF
          ENDIF
      ENDDO
      RETURN
      END
