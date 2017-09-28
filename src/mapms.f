      SUBROUTINE MAPMS(NLOS)

!     ROUTINE MAPMS DETERMINES WHICH ATMOSPHERIC LAYER
!     CONTAINS EACH LINE-OF-SIGHT PATH SEGMENT.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       NLOS    NUMBER OF LINE-OF-SIGHT PATHS.
      INTEGER NLOS

!     COMMONS:

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

!     LOCAL VARIABLES:
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
      INTEGER I,J,ILO,IHI,ILOS
      DOUBLE PRECISION AVGALT

!     LOOP OVER LINES-OF-SIGHT:
      DO ILOS=1,NLOS

!         FIND INTERPOLATION INDEX AND FRACTION FOR SENSOR ALTITUDE:
          ILO=1
          DO IHI=2,ML-1
              IF(ZM(IHI).GT.PTHALT(0,ILOS))GOTO 10
              ILO=IHI
          ENDDO
          IHI=ML
   10     CONTINUE
          LOWALT(0,ILOS)=ILO
          FACALT(0,ILOS)                                                &
     &      =SNGL((PTHALT(0,ILOS)-ZM(ILO))/(ZM(IHI)-ZM(ILO)))
          IF(ABS(FACALT(0,ILOS)).LT..00001)FACALT(0,ILOS)=0.

!         FIND INTERPOLATION INDEX AND FRACTION FOR SEGMENT ENDPOINTS:
          DO J=1,NSEG(ILOS)
              ILO=1
              DO IHI=2,ML-1
                  IF(ZM(IHI).GT.PTHALT(J,ILOS))GOTO 20
                  ILO=IHI
              ENDDO
              IHI=ML
   20         CONTINUE
              LOWALT(J,ILOS)=ILO
              FACALT(J,ILOS)                                            &
     &          =SNGL((PTHALT(J,ILOS)-ZM(ILO))/(ZM(IHI)-ZM(ILO)))
              IF(ABS(FACALT(J,ILOS)).LT..00001)FACALT(J,ILOS)=0.

!             DETERMINE ATMOSPHERIC LAYER CONTAINING SEGMENT MIDPOINTS:
              AVGALT=(PTHALT(J-1,ILOS)+PTHALT(J,ILOS))/2
              MAPPTH(J,ILOS)=1
              DO I=2,ML-1
                  IF(ZM(I).GE.AVGALT)GOTO 30
                  MAPPTH(J,ILOS)=I
              ENDDO
   30         CONTINUE
          ENDDO

!     END LINES-OF-SIGHT LOOP:
      ENDDO
      RETURN
      END