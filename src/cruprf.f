      SUBROUTINE CRUPRF(NCRALT)

!     THIS ROUTINE READS IN USER-DEFINED CLOUD/RAIN MODEL PROFILES.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       NCRALT    NUMBER OF CLOUD/RAIN PROFILES BOUNDARY ALTITUDES
      INTEGER NCRALT

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /CLDRR/
!       ZCLD     INPUT (*,0) & MODEL (*,1) CLOUD ALTITUDE PROFILES [KM].
!       PCLD     CLOUD PROFILE PRESSURE [MB].
!       CLD      INPUT (*,0) & MODEL (*,>0) WATER CLOUD PROFILES [G/M3].
!       CLDICE   INPUT (*,0) & MODEL (*,>0) ICE CLOUD PROFILES [G/M3].
!       RR       INPUT (*,0) & MODEL (*,>0) RAIN RATE PROFILES [MM/HR].
      DOUBLE PRECISION ZCLD
      REAL PCLD,CLD,CLDICE,RR
      COMMON/CLDRR/ZCLD(1:NZCLD,0:1),PCLD,CLD(1:NZCLD,0:5),             &
     &  CLDICE(1:NZCLD,0:1),RR(1:NZCLD,0:5)

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

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /CLDRR/
      EXTERNAL DEVCBD,MDTA

!     LOCAL VARIABLES:
!       TCLD     CLOUD PROFILE TEMPERATURE [K].
!       PCLDSV   PREVIOUS CLOUD PROFILE PRESSURE [MB].
!       ZCLOUD   CLOUD PROFILE ABOVE SEA-LEVEL ALTITUDE [KM].
      INTEGER IALTM1,ICRALT,K,KM1,KSAV,KSAVM1
      REAL TCLD,PCLDSV
      DOUBLE PRECISION ZCLOUD

!     FUNCTIONS:
!       HYSTAT   INTEGRATES HYDROSTATIC EQUATION TO FIND ALTITUDE.
      DOUBLE PRECISION HYSTAT

!     CHECK NUMBER OF BOUNDARY ALTITUDES (NCRALT)
      IF(NCRALT.GT.NZCLD)THEN
          WRITE(IPR,'(/A,I3,A,/(18X,A,I3,A))')' Error in CRUPRF:  The'//&
     &      ' input number of cloud/rain level altitudes (',NCRALT,')', &
     &      ' exceeds parameter NZCLD (=',NZCLD,').  Increase NZCLD in',&
     &      ' "PARAMS.h" & make necessary changes in block data "MDTA".'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error in CRUPRF:  Too many cloud/rain level altitudes!'
      ENDIF

!     WRITE OUT USER-DEFINED CLOUD/RAIN MODEL PROFILES MESSAGE
      IF(.NOT.LJMASS)WRITE(IPR,'(/A,I3,A)')' User-defined cloud/rain'   &
     &  //' model profiles with',NCRALT,' level altitudes.'

!     READ IN DATA AT FIRST ALTITUDE
      ICRALT=1
      IF(MODEL.EQ.8)THEN
          IF(LJMASS)THEN
              CALL INITCD('CARD2E1')
          ELSE
              READ(IRD,'(4F10.0)',ERR=30)                               &
     &          PCLD,CLD(1,0),CLDICE(1,0),RR(1,0)
          ENDIF

!         INTEGRATE HYDROSTATIC EQUATION TO FIND ALTITUDE:
          KM1=1
          IF(PCLD.GT.PPROF(1))PCLD=PPROF(1)
          DO K=2,ML
              IF(PCLD.GE.PPROF(K))THEN
                  TCLD=TPROF(KM1)+(TPROF(K)-TPROF(KM1))                 &
     &              *LOG(PCLD/PPROF(KM1))/LOG(PPROF(K)/PPROF(KM1))

!                 FIND CLOSER PRESSURE LEVEL:
                  IF(PCLD-PPROF(K).LT.PPROF(KM1)-PCLD)THEN
                      ZCLOUD=HYSTAT(ZM(K),TPROF(K),PPROF(K),TCLD,PCLD)
                  ELSE
                      ZCLOUD=HYSTAT(ZM(KM1),                            &
     &                  TPROF(KM1),PPROF(KM1),TCLD,PCLD)
                  ENDIF
                  ZCLD(1,0)=ZCLOUD-ZM(1)
                  IF(ABS(ZCLD(1,0)).LT..0001D0)ZCLD(1,0)=0.D0
                  KSAVM1=KM1
                  PCLDSV=PCLD
                  GOTO 10
              ENDIF
              KM1=K
          ENDDO
          GOTO 30
      ELSEIF(LJMASS)THEN
          CALL INITCD('CARD2E1')
      ELSE
          READ(IRD,'(4F10.0)',ERR=30)                                   &
     &      ZCLD(1,0),CLD(1,0),CLDICE(1,0),RR(1,0)
      ENDIF
   10 CONTINUE

!     THE BASE ALTITUDE (DEFINED RELATIVE TO THE GROUND), THE
!     DENSITIES AND THE RAIN RATE MUST BE ALL BE NON-NEGATIVE
      IF(ZCLD(1,0).LT.0.D0)THEN
          WRITE(IPR,'(/A,F10.5,A)')' Error in CRUPRF:  The cloud/rain'  &
     &      //' profile base altitude is',ZCLD(1,0),' km below ground.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' CLOUD/RAIN PROFILE BASE ALTITUDE IS BELOW GROUND!'
      ENDIF
      IF(CLD(1,0).LT.0. .OR. CLDICE(1,0).LT.0. .OR. RR(1,0).LT.0.)THEN
          WRITE(IPR,'(/A,F10.5,A,/(18X,A,F10.5,A))')' Error in CRUPRF: '&
     &      //' Negative cloud/rain density/rate input at',ZCLD(1,0),   &
     &      ' KM.',' Water droplet density',CLD(1,0),' gm/m3',          &
     &             ' Ice particle density ',CLDICE(1,0),' gm/m3',       &
     &             ' Rain rate            ',RR(1,0),' mm/hr'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error in CRUPRF: Negative cloud/rain density/rate!'
      ENDIF

!     LOOP OVER PROFILE ALTITUDE
      IALTM1=1
      DO ICRALT=2,NCRALT

!         READ IN DATA:
          IF(MODEL.EQ.8)THEN
              IF(LJMASS)THEN
                  CALL INITCD('CARD2E1')
              ELSE
                  READ(IRD,'(4F10.0)',ERR=30)                           &
     &              PCLD,CLD(ICRALT,0),CLDICE(ICRALT,0),RR(ICRALT,0)
              ENDIF

!             INTEGRATE HYDROSTATIC EQUATION TO FIND ALTITUDE:
              KM1=KSAVM1
              KSAV=KM1+1
              DO K=KSAV,ML
                  IF(PCLDSV.LE.PCLD)THEN
                      WRITE(IPR,'(/A,/(18X,A,F10.5,A))')                &
     &                  ' Error in CRUPRF: MODTRAN does'//              &
     &                  ' not allow pressure inversions.',              &
     &                  ' Pressure at previous level',PCLDSV,' mbar.',  &
     &                  ' Pressure at  current level',PCLD,' mbar.'
                      IF(LJMASS)CALL WRTBUF(FATAL)
                      STOP 'Error in CRUPRF:  Pressure inversion!'
                  ENDIF
                  IF(PCLD.GE.PPROF(K))THEN
                      TCLD=TPROF(KM1)+(TPROF(K)-TPROF(KM1))             &
     &                  *LOG(PCLD/PPROF(KM1))/LOG(PPROF(K)/PPROF(KM1))
                      IF(PCLD-PPROF(K).LT.PPROF(KM1)-PCLD)THEN
                          ZCLOUD=HYSTAT(ZM(K),                          &
     &                      TPROF(K),PPROF(K),TCLD,PCLD)
                      ELSE
                          ZCLOUD=HYSTAT(ZM(KM1),                        &
     &                      TPROF(KM1),PPROF(KM1),TCLD,PCLD)
                      ENDIF
                      ZCLD(ICRALT,0)=ZCLOUD-ZM(1)
                      KSAVM1=KM1
                      PCLDSV=PCLD
                      GOTO 20
                  ENDIF
                  KM1=K
              ENDDO
              GOTO 30
          ELSEIF(LJMASS)THEN
              CALL INITCD('CARD2E1')
          ELSE
              READ(IRD,'(4F10.0)',ERR=30)ZCLD(ICRALT,0),CLD(ICRALT,0),  &
     &          CLDICE(ICRALT,0),RR(ICRALT,0)
          ENDIF
   20     CONTINUE

!         CHECK FOR INCREASING ALTITUDES.
          IF(ZCLD(ICRALT,0).LE.ZCLD(IALTM1,0))THEN
              WRITE(IPR,'(/A,/18X,A,/(18X,I5,F10.5))')' Error in'//     &
     &          ' CRUPRF:  The cloud/rain profile altitudes must be',   &
     &          ' read in increasing order.  Values thus far are:',     &
     &          (IALTM1,ZCLD(IALTM1,0),IALTM1=1,ICRALT)
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' Cloud/rain altitudes not monotonically increasing!'
          ENDIF

!         CHECK FOR NEGATIVE CLOUD DENSITIES OR RAIN RATES.
          IF(CLD(ICRALT,0).LT.0. .OR. CLDICE(ICRALT,0).LT.0.            &
     &      .OR. RR(ICRALT,0).LT.0.)THEN
              WRITE(IPR,'(/A,F10.5,A,/(18X,A,F10.5,A))')                &
     &          ' Error in CRUPRF:  Negative cloud/rain'//              &
     &          ' density/rate input at',ZCLD(ICRALT,0),' KM.',         &
     &          ' WATER DROPLET DENSITY',CLD(ICRALT,0),' GM/M3',        &
     &          ' ICE PARTICLE DENSITY ',CLDICE(ICRALT,0),' GM/M3',     &
     &          ' RAIN RATE            ',RR(ICRALT,0),' MM/HR'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'Error in CRUPRF:  Negative cloud/rain density/rate!'
          ENDIF
          IALTM1=ICRALT
      ENDDO

!     RETURN TO CRPROF
      RETURN

!     ERROR READING CLOUD/RAIN DATA
   30 CONTINUE
      WRITE(IPR,'(/A,I3)')' Error in CRUPRF:  Unable to read'           &
     &  //' CLOUD/RAIN PROFILE DATA for layer',ICRALT
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP 'ERROR reading CLOUD/RAIN PROFILE DATA in routine CRUPRF!'
      END
