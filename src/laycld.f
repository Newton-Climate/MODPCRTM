      SUBROUTINE LAYCLD(K,EQLWCZ,RRATZ,ICLD1,GNDALT)

!     ROUTINE TO DETERMINE CLOUD DENSITY AND RAIN RATE FOR LAYER K
!       ZM      COMMON /MPROF/ FINAL ALTITUDES FOR LOWTRAN
!       ZK      EFFECTIVE CLOUD ALTITUDES
!       ZCLD    CLOUD ALTITUDE ARRAY
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       K        LAYER INDEX.
!       GNDALT   GROUND ALTITUDE [KM].
      INTEGER K,ICLD1
      REAL EQLWCZ,RRATZ
      DOUBLE PRECISION GNDALT

!     COMMONS:

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

!     LOCAL VARIABLES:
      INTEGER MC,MR,MKM1,MK
      REAL FAC
      DOUBLE PRECISION ZK

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /CLDRR/
      EXTERNAL MDTA
      IF(ICLD1.LE.0 .OR. ICLD1.GT.11)RETURN
      ZK=ZM(K)-GNDALT
      IF(ZK.LT.0.D0)ZK=0.D0
      IF(ICLD1.LE.5)THEN

!         ICLD1 IS 1- 5 ONE OF 5 SPECIFIC CLOUD MODELS IS CHOSEN
          MC=ICLD1
          MR=6
      ELSE

!         ICLD1 IS 6-11 ONE OF 5 SPECIFIC CLOUD/RAIN MODELS CHOSEN
          MC=1
          IF(ICLD1.EQ.6)MC=3
          IF(ICLD1.EQ.7 .OR. ICLD1.EQ.8)MC=5
          MR=ICLD1-5
      ENDIF
      EQLWCZ=0.
      RRATZ=0.
      IF(ZK.GE.ZCLD(1,1))THEN
          MKM1=1
          DO MK=2,16
              IF(ZK.LE.ZCLD(MK,1))THEN
                  FAC=SNGL((ZCLD(MK,1)-ZK)/(ZCLD(MK,1)-ZCLD(MKM1,1)))
                  EQLWCZ=CLD(MK,MC)+FAC*(CLD(MKM1,MC)-CLD(MK,MC))
                  IF(MR.LE.5)RRATZ=RR(MK,MR)+FAC*(RR(MKM1,MR)-RR(MK,MR))
                  RETURN
              ENDIF
              MKM1=MK
          ENDDO
      ENDIF
      RETURN
      END
