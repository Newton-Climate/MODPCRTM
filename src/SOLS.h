
!     /SOLS/
!       ECALAY   PATH SEGMENT EARTH CENTER ANGLE IN ALTITUDE LAYER.
!       ZENLAY   PATH ZENITH (NADIR) ANGLE IN ALTITUDE LAYER.
!       JTURN    PATH SEGMENT INDEX OF TANGENT POINT (IF > 0).
!       LAYSEG   MAP FROM PATH SEGMENT INDEX TO ALTITUDE LAYER INDEX.
!       NANGLS   NUMBER OF ANGLES IN USER-DEFINED PHASE FUNCTIONS.
!       ANGF     GRID OF ANGLES FOR USER-DEFINED PHASE FUNCTIONS [DEG].
      DOUBLE PRECISION ECALAY,ZENLAY
      INTEGER JTURN,LAYSEG,NANGLS
      REAL ANGF
      COMMON/SOLS/ECALAY(LAYDM1),ZENLAY(LAYDM1),JTURN,LAYSEG(LAYTWO+1), &
     &  NANGLS,ANGF(MANGLS)
