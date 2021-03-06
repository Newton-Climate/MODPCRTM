      SUBROUTINE SMPREP(REARTH,H1,H2,ANGLE,HRANGE,BETA,ISLCT)

!     THIS SUBROUTINE PREPS OR PREPROCESSES PATH GEOMETRY
!     PARAMETERS TO SEE IS IF THE HRANGE IS SMALL.  ALL SMALL PATH
!     CASES ARE CAST INTO THE EQUIVALENT CASE 2C (H1,H2,HRANGE).
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       REARTH   RADIUS OF THE EARTH [KM].
      DOUBLE PRECISION REARTH,H1,H2,ANGLE,HRANGE,BETA
      INTEGER ISLCT

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD,ATMCON
      REAL AORIG
      COMMON/SMALL5/AORIG

!     DECLARE LOCAL VARIABLES
      DOUBLE PRECISION R1,STOR

!     AORIG IS THE ORIGINAL INPUT ANGLE IN SINGLE PRECISION
      AORIG=SNGL(ANGLE)

!     IF ALREADY IN CASE 2C (H1,H2,HRANGE) RETURN.
      IF(ISLCT.EQ.23)RETURN

!     DETERMINE HRANGE WITHOUT REFRACTION
!     (REFRACTION IS NOT IMPORTANT FOR SHORT PATHS).
      R1=H1+REARTH
      IF(ISLCT.EQ.24)THEN

!         CASE 2D (H1,H2,BETA)
          HRANGE=SQRT((H1-H2)**2                                        &
     &      +4*R1*(H2+REARTH)*SIN(BETA/(2*DPDEG))**2)
      ELSEIF(ISLCT.EQ.22)THEN

!         CASE 2B (H1,ANGLE,HRANGE)
          STOR=(HRANGE-H1)*(HRANGE+H1)+2*R1*(H1+HRANGE*COS(ANGLE/DPDEG))
          H2=STOR/(SQRT(REARTH**2+STOR)+REARTH)
      ENDIF

!     CHECK FOR SMALL HRANGE (JMASS TREATS MESSAGE AS WARNING):
      IF(HRANGE.LE.RSMALL .AND. HRANGE.GT.0.D0)THEN
          ANGLE=0.D0
          BETA=0.D0
          WRITE(IPR,'(/2A,F12.7,A)')' FROM SMPREP:  HRANGE IS',         &
     &      ' BELOW THE SMALL HRANGE CUTOFF (',RSMALL,'KM).'
          WRITE(IPR,'(/14X,A,/14X,4(A,F12.7))')                         &
     &      ' THE INPUT GEOMETRY HAS BEEN CONVERTED TO CASE 2C WITH',   &
     &      ' H1 =',H1,'KM, H2 =',H2,'KM, AND HRANGE =',HRANGE,'KM.'
      ELSEIF(HRANGE.LT.0.)THEN
          WRITE(IPR,'(/A)')' FROM SMPREP:  HRANGE IS LESS THAN 0.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'FROM SMPREP:  HRANGE IS LESS THAN 0.'
      ENDIF
      RETURN
      END
