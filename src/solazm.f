      DOUBLE PRECISION FUNCTION SOLAZM(REL_AZ,DELO,BETA,IARB,IARBO)

!     FUNCTION SOLAZM RETURNS THE VALUE OF SOLAR AZIMUTH RELATIVE TO
!     THE LINE OF SIGHT, AT THE CURRENT SCATTERING LOCATION
      IMPLICIT NONE

!     PARAMETERS:
!       EPSILN   SMALL ANGLE CUT-OFF [DEG].
      DOUBLE PRECISION EPSILN
      PARAMETER(EPSILN=1.D-5)

!     INPUT ARGUMENTS:
!       REL_AZ   RELATIVE SOLAR AZIMUTH ANGLE AT SENSOR [DEG].
!       DELO     OBSERVER-TO-SUN EARTH CENTERED ANGLE [DEG].
!       BETA     OBSERVER-TO-SCATTER POINT EARTH CENTERED ANGLE [DEG].
!       IARBO    INITIAL GUESS AT CASE INDEX.
      DOUBLE PRECISION REL_AZ,DELO,BETA
      INTEGER IARBO

!     OUTPUT ARGUMENTS:
!       IARBO    ACTUAL CASE INDEX.
      INTEGER IARB

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     LOCAL VARIABLES:
!       DELOR    OBSERVER-TO-SUN EARTH CENTERED ANGLE [RAD].
!       SNDELO   SIN OF OBSERVER-TO-SUN EARTH CENTERED ANGLE.
!       RELAZR   RELATIVE SOLAR AZIMUTH ANGLE AT SENSOR [RAD].
!       BETAR    OBSERVER-TO-SCATTER POINT EARTH CENTERED ANGLE [RAD].
!       DENOM    DENOMINATOR OF ARC-TANGENT ARGUMENT.
      DOUBLE PRECISION DELOR,SNDELO,RELAZR,BETAR,DENOM

!     CASE BRANCHING:
      IF(IARBO.EQ.0)THEN

!         INITIALIZE SOLAZM AND CASE (IARB) FOR GENERAL CASE:
          SOLAZM=REL_AZ
          IARB=0

!         SPECIAL CASES:  NUMERATOR NEAR ZERO IN THE FOLLOWING 3 CASES.
          IF(DELO.LE.EPSILN)THEN

!     1)      DELO=0.
              IF(BETA.LE.EPSILN)THEN
                  IARB=2
              ELSE
                  SOLAZM=180.D0
              ENDIF
          ELSEIF(ABS(REL_AZ).LE.EPSILN)THEN

!     2)      REL_AZ=0.
              IF(ABS(BETA-DELO).LT.EPSILN)THEN

!                 SCATTERING POINT IS DIRECTLY UNDER THE SUN
                  IARB=2
              ELSEIF(BETA.LT.DELO)THEN
                  SOLAZM=0.D0
              ELSE
                  SOLAZM=180.D0
              ENDIF
          ELSEIF(ABS(REL_AZ).GE.180-EPSILN)THEN

!     3)      REL_AZ=180.
              SOLAZM=180.D0
          ELSE
              DELOR=DELO/DPDEG
              SNDELO=SIN(DELOR)
              BETAR=BETA/DPDEG
              RELAZR=REL_AZ/DPDEG
              DENOM=COS(BETAR)*SNDELO*COS(RELAZR)-SIN(BETAR)*COS(DELOR)
              IF(ABS(DENOM).LE.EPSILN)THEN

!                 DENOMINATOR CAN STILL GO TO ZERO:
                  IF(REL_AZ.LT.0.)THEN
                      SOLAZM=-90.D0
                  ELSE
                      SOLAZM=90.D0
                  ENDIF
              ELSE

!                 GENERAL CASE
                  SOLAZM=DPDEG*ATAN(SNDELO*SIN(RELAZR)/DENOM)

!                 |SOLAZM|<90, BUT SOLAZM & REL_AZ SHOULD HAVE SAME SIGN
                  IF(REL_AZ.GT.0.D0 .AND. SOLAZM.LT.0.D0)THEN
                      SOLAZM=SOLAZM+180
                  ELSEIF(REL_AZ.LT.0.D0 .AND. SOLAZM.GT.0.D0)THEN
                      SOLAZM=SOLAZM-180
                  ENDIF
              ENDIF
          ENDIF

!     SPECIAL CASES WHEN REL_AZ IS ARBITRARY
      ELSEIF(IARBO.EQ.1 .OR. IARBO.EQ.3 .OR. BETA.LE.EPSILN)THEN
          IARB=IARBO
          SOLAZM=REL_AZ
      ELSE

!         SOLAZM=180. (MOVED OUT FROM UNDER THE SUN)
          IARB=0
          SOLAZM=180.D0
      ENDIF
      RETURN
      END
