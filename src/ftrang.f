      SUBROUTINE FTRANG(REARTH,H1,H2,HRANGE,ANGLE,PHI,LENN,HMIN,IERROR)

!     FTRANG CALCULATES THE ZENITH ANGLE AT H1 (ANGLE)
!     AND AT H2 (PHI) GIVEN H1, H2 AND HRANGE.
      IMPLICIT NONE

!     PARAMETERS:
!       ITERMX   MAXIMUM NUMBER OF ITERATIONS.
!       TOLRNC   MAXIMUM ALLOWED RANGE ABSOLUTE ERROR [KM].
      INTEGER ITERMX
      REAL TOLRNC
      PARAMETER(ITERMX=35,TOLRNC=.001)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       REARTH   RADIUS OF THE EARTH [KM].
!       H1       OBSERVER ALTITUDE [KM].
!       H2       FINAL ALTITUDE [KM].
!       HRANGE   SLANT PATH RANGE [KM].

!     OUTPUT ARGUMENTS:
!       HRANGE   CONVERGED SLANT PATH RANGE [KM].
!       ANGLE    ZENITH ANGLE AT H1 TOWARDS H2 [DEG].
!       PHI      ZENITH ANGLE AT H2 TOWARDS H1 [DEG].
!       LENN     PATH LENGTH SWITCH (=0 SHORT PATH,
!                  =1 FOR PATH THROUGH TANGENT POINT).
!       HMIN     PATH MINIMUM ALTITUDE [KM].
!       IERROR   ERROR FLAG (=0 FOR SUCCESS, =1 FOR FAILURE).
      DOUBLE PRECISION REARTH,H1,H2,HRANGE,ANGLE,PHI,HMIN
      INTEGER LENN,IERROR

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       LOWER    FLAG, TRUE WHEN ANGLE INCREMENT NEEDS TO BE DECREASED.
!       DANGLE   ANGLE INCREMENT [DEG].
!       HA       THE LOWER OF H1 AND H2 [KM].
!       HB       THE HIGHER OF H1 AND H2 [KM].
!       HBSAV    SAVED VALUE OF HB [KM].
!       ANGLE1   CURRENT GUESS FOR ANGLE [DEG].
!       ITER     ITERATION COUNTER.
      LOGICAL LOWER
      DOUBLE PRECISION DANGLE,HA,HB,HBSAV,RA,COEF,STORE,ANGLE1,         &
     &  BETA,RANGE1,BEND,TERM1,TERM2,ADDRAN,ADDANG,ARG
      INTEGER ITER

!     INITIALIZE PARAMETERS
      LOWER=.FALSE.
      DANGLE=.1D0
      IF(H1.GT.H2)THEN
          HA=H2
          HB=H1
      ELSE
          HA=H1
          HB=H2
      ENDIF

!     HBSAV IS DEFINED TO PROTECT AGAINST RPATH REDEFINING HB
      HBSAV=HB

!     GUESS AT ANGLE, INTEGRATE OVER PATH TO FIND RANGE, TEST
!     FOR CONVERGENCE, AND ITERATE ANGLE IF NECESSARY.
!     CALCULATE THE NO REFRACTION ZENITH ANGLE AS AN INITIAL GUESS
      RA=REARTH+HA
      COEF=.5D0/RA
      STORE=(HB-HA)*(REARTH+HB+RA)
      ARG=COEF*(STORE/HRANGE-HRANGE)
      IF(ABS(ARG).LT.1.D0)THEN

!         COSINE WELL DEFINED:
          ANGLE=DPDEG*ACOS(ARG)
      ELSEIF(ABS(ARG).LE.1.00001D0)THEN

!         POSSIBLE NUMERICAL INACCURACY: SET TO NADIR OR ZENITH.
          IF(H1.GT.H2)THEN
              HRANGE=H1-H2
              ANGLE=180.D0
              PHI=0.D0
              HMIN=H2
          ELSE
              HRANGE=H2-H1
              ANGLE=0.D0
              PHI=180.D0
              HMIN=H1
          ENDIF
          LENN=0
          RETURN
      ELSE

!         INPUT RANGE TOO SHORT:
          WRITE(IPR,'(/A,F12.5,A,/3(A,F12.5))')'Error in FTRANG: '      &
     &      //' Input slant range (',HRANGE,'km) less than altitude'    &
     &      //' difference',' between sensor (',H1,'km) and final (',   &
     &      H2,'km) altitudes.'
          IERROR=1
          RETURN
      ENDIF
      ANGLE1=ANGLE

!     TABLE HEADER:
      IF(.NOT.LJMASS)WRITE(IPR,'(///A,3F10.5,A,//A,//(A))')             &
     &  ' CASE 2C: GIVEN H1, H2, HRANGE:  (',H1,H2,HRANGE,' )',         &
     &  ' ITERATE AROUND ANGLE UNTIL RANGE CONVERGES',                  &
     &  ' ITER       ANGLE         RANGE        DRANGE     '//          &
     &  '     BETA          HMIN           PHI       BENDING',          &
     &  '            (DEG)          (KM)          (KM)     '//          &
     &  '    (DEG)          (KM)         (DEG)         (DEG)'

!     BEGIN ITERATIVE PROCEDURE
      ITER=0
   10 CONTINUE
      ITER=ITER+1
          IF(ITER.GT.ITERMX)THEN

!             JMASS TREATS FOLLOWING MESSAGE AS WARNING
              WRITE(IPR,'(/A,//10X,3(A,F13.6,4X),A,I5,                  &
     &          2(//10X,A),F16.9,4X,A,F16.9)')'FTRANG,'                 &
     &          //' CASE 2C (H1,H2,HRANGE): SOLUTION DID NOT CONVERGE', &
     &          'H1 =',H1,'H2 =',H2,'BETA  =',BETA, 'ITERATIONS =',ITER,&
     &         'LAST ITERATION','ANGLE =',ANGLE1,'HRANGE =',RANGE1
              IERROR=1
              RETURN
          ENDIF

!         DETERMINE RANGE1, THE RANGE CORRESPONDING TO ANGLE1
          HB=HBSAV
          IF(ANGLE1.LE.90.D0)THEN

!             SHORT UPWARD PATH
              LENN=0
              CALL FINDMN(REARTH,HA,ANGLE1,HB,                          &
     &          LENN,ITER,HMIN,PHI,IERROR,.TRUE.)
              CALL RFPATH(REARTH,HA,HB,ANGLE1,LENN,HMIN,.FALSE.,        &
     &          PHI,BETA,RANGE1,BEND)
          ELSE

!             LONG PATH
              LENN=1
              CALL FINDMN(REARTH,HA,ANGLE1,HB,                          &
     &          LENN,ITER,HMIN,PHI,IERROR,.TRUE.)
              IF(HA.LE.HMIN)THEN

!                 BECAUSE ANGLE1 EXCEEDS 90 DEGREES BY SO
!                 LITTLE, HA DOES NOT EXCEED HMIN.  CORRECT
!                 BY DECREASING ANGLE1 TO 90 DEGREES, SETTING
!                 HMIN TO HA, AND CHANGING LENN FROM 1 TO 0.
                  HMIN=HA
                  ANGLE1=90.D0
                  LENN=0
                  CALL RFPATH(REARTH,HA,HB,ANGLE1,LENN,HMIN,.FALSE.,    &
     &              PHI,BETA,RANGE1,BEND)
              ELSE
                  CALL RFPATH(REARTH,HA,HB,ANGLE1,LENN,HMIN,.FALSE.,    &
     &              PHI,BETA,RANGE1,BEND)
                  IF(LENN.EQ.0)THEN

!                     PATH INTERSECTED THE EARTH.  DECREASE ANGLE BY .1
                      IF(.NOT.LJMASS)WRITE(IPR,'(I3,7F14.7)')           &
     &                   ITER,ANGLE1,RANGE1,HRANGE-RANGE1,              &
     &                   BETA,HMIN,PHI,BEND
                      LOWER=.TRUE.
                      ANGLE1=ANGLE1-DANGLE
                      GOTO 10
                  ENDIF
              ENDIF
          ENDIF

!         IF THE FINAL ALTITUDE HB HAS BEEN LOWERED (BECAUSE
!         HB WAS ABOVE THE TOP OF THE ATMOSPHERE), ADD ON THE
!         REMAINDER OF PATH LENGTH ASSUMING NO REFRACTION.
          IF(HBSAV.GT.HB)THEN
              TERM1=(REARTH+HB)*COS(PHI/DPDEG)
              TERM2=(HBSAV-HB)*(2*REARTH+HBSAV+HB)
              ADDRAN=TERM2/(SQRT(TERM1*TERM1+TERM2)-TERM1)
              RANGE1=RANGE1+ADDRAN
              ADDANG=DPDEG*ASIN(ADDRAN*SIN(PHI/DPDEG)/(REARTH+HBSAV))
              BETA=BETA+ADDANG
              PHI=PHI+ADDANG
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,'(I3,7F14.7)')                       &
     &      ITER,ANGLE1,RANGE1,HRANGE-RANGE1,BETA,HMIN,PHI,BEND

!         CHECK FOR CONVERGENCE
          IF(ABS(HRANGE-RANGE1).LT.TOLRNC .OR.                          &
     &      ABS(1-RANGE1/HRANGE).LT.2.D-6)THEN
              HRANGE=RANGE1
              IF(H1.LE.H2)THEN
                  ANGLE=ANGLE1
              ELSE
                  ANGLE=PHI
                  PHI=ANGLE1
              ENDIF
              IF(HMIN.GE.GNDALT)RETURN

!             JMASS TREATS FOLLOWING MESSAGE AS WARNING
              WRITE(IPR,'(/A,//8X,2A)')                                 &
     &          'FTRANG, CASE 2C (H1,H2,HRANGE):  HRANGE IS TOO LARGE.',&
     &          ' REFRACTED PATH TANGENT HEIGHT IS LESS THAN',          &
     &          ' GROUND ALTITUDE; THE PATH INTERSECTS THE EARTH.'
              IERROR=1
              RETURN
          ENDIF

!         DETERMINE NEW VALUE FOR ANGLE (.6 IS A FUDGE FACTOR INTRODUCED
!         TO AVOID OVER CORRECTING AND SPEED UP CONVERGENCE).
          IF(RANGE1.LE.0.)THEN
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'FTRANG error:  HRange is zero or less.'
          ENDIF
          ARG=COEF*(STORE/RANGE1-RANGE1)
          IF(ARG.GE.1.D0)THEN
              ANGLE1=ANGLE1+.6D0*ANGLE
          ELSEIF(ARG.GT.-1.D0)THEN
              ANGLE1=ANGLE1+.6D0*(ANGLE-DPDEG*ACOS(ARG))
          ELSE
              ANGLE1=ANGLE1+.6D0*(ANGLE-180)
          ENDIF

!         CHECK IF ANGLE INCREMENT MUST BE LOWERED
          IF(LOWER)THEN
              DANGLE=DANGLE/5
              LOWER=.FALSE.
          ENDIF
      GOTO 10
      END
