      SUBROUTINE IRFXN(X,RX,RXRAT)

!     THIS ROUTINE FINDS INDEX OF REFRACTION AND ITS DERIVATIVE AT X.
      IMPLICIT NONE

      INCLUDE 'PARAMS.h'
!     X = HEIGHT FROM THE SURFACE OF THE EARTH
!     RX = INDEX OF REFRACTION
!     DRX = DERIVATIVE
!     RXRAT = DRX/RX

      DOUBLE PRECISION  X, RX, RXRAT, XN, H
      INTEGER I, J

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

      DO J=2,ML-1
         IF(X.LE.ZM(J))GOTO 10
      ENDDO
  10  CONTINUE
      I=J-1

      H=DBLE(LOG((RFNDX(J)+1.E-30)/(RFNDX(I)+1.E-30)))/(ZM(J)-ZM(I))
      XN=DBLE(RFNDX(I))*EXP((X-ZM(I))*H)
      RX=1+XN
      RXRAT=H*XN/RX
      RETURN
      END
