      SUBROUTINE CD2B

!     PROCESS CARD2B INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

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

!     /CARD2B/
      REAL ZCVSA,ZTVSA,ZINVSA
      COMMON/CARD2B/ZCVSA,ZTVSA,ZINVSA

!     /ZVSALY/
      DOUBLE PRECISION ZVSA
      REAL RHVSA,AHVSA
      INTEGER IHVSA
      COMMON/ZVSALY/ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)

!     CHECK FOR VERTICAL STRUCTURE ALGORITHM (VSA) OPTION:
      IF(IVSA.EQ.1)THEN
          IF(LJMASS)THEN
              CALL INITCD('CARD2B')
          ELSE
              READ(IRD,'(3F10.0)')ZCVSA,ZTVSA,ZINVSA
              WRITE(IPR,'(/A,3F10.5)')                                  &
     &          ' CARD 2B *****',ZCVSA,ZTVSA,ZINVSA
          ENDIF
          CALL VSA(IHAZE,VIS,ZCVSA,ZTVSA,ZINVSA,ZVSA,RHVSA,AHVSA,IHVSA)
      ELSE
          ZCVSA=-99.
          ZTVSA=-99.
          ZINVSA=-99.
      ENDIF

!     RETURN TO DRIVER:
      RETURN
      END
