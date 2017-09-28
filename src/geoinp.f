      SUBROUTINE GEOINP(REARTH,ZMAX,H1ALT,H2ALT,OBSZEN,HRANGE,BETA,     &
     &  ITYPE,LENN,HMIN,BCKZEN,IERROR,ISLCT)

!     GEOINP INTERPRETS ALLOWABLE COMBINATIONS OF INPUT PATH PARAMETERS
!     INTO THE STANDARD SET:  H1ALT, H2ALT, OBSZEN, BCKZEN, HMIN AND
!     LENN.  THE ALLOWABLE COMBINATIONS OF INPUT PARAMETERS ARE

!          FOR ITYPE = 2  (SLANT PATH H1ALT TO H2ALT)
!              A.  H1ALT, H2ALT AND OBSZEN;
!              B.  H1ALT, OBSZEN AND HRANGE;
!              C.  H1ALT, H2ALT AND HRANGE; OR
!              D.  H1ALT, H2ALT, AND BETA.

!          FOR ITYPE = 3  (SLANT PATH H1ALT TO SPACE/GROUND)
!              A. H1ALT AND OBSZEN; OR
!              B. H1ALT AND HMIN (INPUT AS H2ALT).

!     THIS ROUTINE ALSO DETECTS BAD INPUT (IMPOSSIBLE GEOMETRY) AND
!     ITYPE = 2 CASES WHICH INTERSECT THE EARTH, AND RETURNS THESE
!     CASES WITH ERROR FLAGS.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       REARTH   RADIUS OF THE EARTH [KM].
!       ZMAX     MAXIMUM ATMOSPHERIC PROFILE ALTITUDE [KM].
!       H1ALT    ALTITUDE OF OBSERVER [KM].
!       H2ALT    FINAL OR TANGENT ALTITUDE [KM].
!       OBSZEN   PATH ZENITH ANGLE FROM OBSERVER TOWARDS SPH2 [DEG].
!       HRANGE   SLANT PATH RANGE [KM].
!       BETA     SLANT PATH EARTH CENTER ANGLE [DEG].
!       ITYPE    FLAG FOR GEOMETRY TYPE.
!       LENN     CASE 2A LENGTH FLAG USED WHEN H2ALT<H1ALT AND ANGLE>90
!                (=0 FOR SHORT PATH, =1 FOR PATH THROUGH TANGENT POINT).
!       HMIN     SLANT PATH MINIMUM ALTITUDE [KM].
!       BCKZEN   PATH ZENITH ANGLE ALONG BACKWARD PATH FROM FINAL
!                ALTITUDE TOWARDS OBSERVER [DEG].
!       IERROR   ERROR FLAG (=0 FOR ACCEPTABLE GEOMETRY INPUTS).
!       ISLCT    PATH TYPE LABEL FOR OPTICAL PATH WITH ITYPE=2.

!     OUTPUT ARGUMENTS:
!       H1ALT    SET TO MIN(H1ALT,ZMAX) [KM].
!       H2ALT    DEFINED FOR CASE 2B (H1ALT, OBSZEN, RANGE) OR SET TO
!                ZMAX IF ITYPE=3 OR H2ALT>ZMAX [KM].
!       OBSZEN   DEFINED FOR CASES 2C (H1ALT, H2ALT, HRANGE),
!                2D (H1ALT, H2ALT, BETA) AND 3B (H1ALT AND HMIN)
!                OR IF H1ALT HAS TO BE REDUCED TO ZMAX.
!       LENN     CASE 2A LENGTH FLAG USED WHEN H2ALT<H1ALT AND ANGLE>90
!                (=0 FOR SHORT PATH, =1 FOR PATH THROUGH TANGENT POINT).
!       HMIN     SLANT PATH MINIMUM ALTITUDE [KM].
!       BCKZEN   PATH ZENITH ANGLE ALONG BACKWARD PATH FROM FINAL
!                ALTITUDE TOWARDS OBSERVER [DEG].
      DOUBLE PRECISION REARTH,ZMAX,H1ALT,H2ALT,OBSZEN,HRANGE,BETA,HMIN, &
     &  BCKZEN
      INTEGER ITYPE,LENN,IERROR,ISLCT

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      DOUBLE PRECISION H2ST

!     BRANCH BASED ON PATH TYPE.
      IF(ITYPE.EQ.3)THEN

!         SLANT PATH TO SPACE.
          IF(H2ALT.EQ.0.D0)THEN

!             CASE 3A:  H1ALT AND OBSZEN.
              IF(.NOT.LJMASS .AND. NPR.LT.1)WRITE(IPR,'(//A)')          &
     &          ' CASE 3A:  GIVEN H1ALT,H2ALT=SPACE,OBSZEN'
              H2ALT=ZMAX
              IF(OBSZEN.GT.90.)LENN=1
              CALL FINDMN(REARTH,H1ALT,OBSZEN,H2ALT,                    &
     &          LENN,0,HMIN,BCKZEN,IERROR,.TRUE.)
          ELSE

!             CASE 3B:  H1ALT AND HMIN
              IF(.NOT.LJMASS .AND. NPR.LT.1)WRITE(IPR,'(//A)')          &
     &          ' CASE 3B:  GIVEN H1ALT, HMIN, H2ALT=SPACE'
              HMIN=H2ALT
              H2ALT=ZMAX
              IF(H1ALT.LT.HMIN)THEN
                  WRITE(IPR,'(/2A,//10X,2(A,F13.6))')' GEOINP,',        &
     &              ' CASE 3B (H1ALT,HMIN,SPACE):  ERROR IN INPUT DATA',&
     &              'H1ALT =',H1ALT,'    IS LESS THAN HMIN =',HMIN
                  IERROR=1
              ELSE
                  CALL FINDMN(REARTH,HMIN,90.D0,H1ALT,                  &
     &              LENN,0,HMIN,OBSZEN,IERROR,.TRUE.)
                  CALL FINDMN(REARTH,HMIN,90.D0,H2ALT,                  &
     &              LENN,0,HMIN,BCKZEN,IERROR,.TRUE.)
                  LENN=1
                  IF(HMIN.EQ.H1ALT)LENN=0
              ENDIF
          ENDIF
      ELSEIF(ITYPE.EQ.2)THEN

!         SLANT PATH TO BETWEEN ALTITUDES.
          IF(ISLCT.EQ.1 .OR.                                            &
     &      (ISLCT.EQ.0 .AND. HRANGE.LE.0.D0 .AND. BETA.LE.0.D0))THEN

!             CASE 2A:  H1ALT, H2ALT, OBSZEN
              IF(.NOT.LJMASS .AND. NPR.LT.1)WRITE(IPR,'(//A)')          &
     &          ' CASE 2A:  GIVEN H1ALT, H2ALT AND OBSZEN'
              IF(H1ALT.GE.H2ALT .AND. OBSZEN.LE.90.)THEN
                  WRITE(IPR,'(/A,//(10X,2(A,F13.6),A))')' GEOINP, CASE' &
     &              //' 2A (H1ALT,H2ALT,OBSZEN):  ERROR IN INPUT DATA', &
     &              'H1ALT (=',H1ALT,'KM) IS GREATER THAN OR EQUAL TO'  &
     &              //' H2ALT (= ',H2ALT,'KM)','AND OBSZEN (=',OBSZEN,  &
     &              'DEG) DOES NOT EXCEED 90DEG.'
                  IERROR=1
              ELSEIF(H1ALT.LE.GNDALT .AND. OBSZEN.GT.90.D0)THEN
                  WRITE(IPR,'(/2A)')                                    &
     &              ' GEOINP, ITYPE = 2: SLANT PATH INTERSECTS',        &
     &              ' THE EARTH OR GNDALT, AND CANNOT REACH H2ALT.'
                  IERROR=1
              ELSE
                  IF(.NOT.LJMASS .AND. H2ALT.LT.H1ALT                   &
     &              .AND. OBSZEN.GT.90.D0 .AND. NPR.LT.1)               &
     &              WRITE(IPR,'(//3A,I3)')' Either a short path',       &
     &              ' (LENN=0) or a long path through a tangent',       &
     &              ' height (LENN=1) is possible:  LENN = ',LENN
                  H2ST=H2ALT
                  CALL FINDMN(REARTH,H1ALT,OBSZEN,H2ALT,                &
     &              LENN,0,HMIN,BCKZEN,IERROR,.TRUE.)
                  IF(H2ALT.NE.H2ST)THEN
                      WRITE(IPR,'(/2A)')                                &
     &                  ' GEOINP, ITYPE = 2: SLANT PATH INTERSECTS',    &
     &                  ' THE EARTH OR GNDALT, AND CANNOT REACH H2ALT.'
                      IERROR=1
                  ENDIF
              ENDIF
          ELSEIF(ISLCT.EQ.4 .OR. (ISLCT.EQ.0 .AND. BETA.GT.0.D0))THEN

!             CASE 2D:  H1ALT, H2ALT, BETA
              CALL FDBETA(REARTH,H1ALT,H2ALT,BETA,OBSZEN,BCKZEN,        &
     &          LENN,HMIN,IERROR)
          ELSEIF(ISLCT.EQ.2 .OR. (ISLCT.EQ.0 .AND. OBSZEN.GT.0.D0))THEN

!             CASE 2B:  H1ALT, OBSZEN, HRANGE
              CALL NEWH2(REARTH,ZMAX,                                   &
     &          H1ALT,H2ALT,OBSZEN,HRANGE,BETA,LENN,HMIN,BCKZEN)
              IF(OBSZEN.GT.90.D0 .AND. BCKZEN.GT.90.D0)LENN=1
              CALL FINDMN(REARTH,H1ALT,OBSZEN,H2ALT,                    &
     &          LENN,0,HMIN,BCKZEN,IERROR,.TRUE.)
          ELSE

!             CASE 2C:  H1ALT, H2ALT, HRANGE
              CALL FTRANG(REARTH,H1ALT,H2ALT,HRANGE,                    &
     &          OBSZEN,BCKZEN,LENN,HMIN,IERROR)
          ENDIF
      ELSE
          WRITE(IPR,'(/2A,I10)')' GEOINP:  ERROR IN INPUT DATA,',       &
     &      ' ITYPE NOT EQUAL TO 2 OR 3.   ITYPE =',ITYPE
          IERROR=1
      ENDIF

!     TEST IERROR AND RECHECK LENN
      IF(IERROR.EQ.0)THEN
          LENN=0
          IF(ABS(HMIN-MIN(H1ALT,H2ALT)).GT..00005)LENN=1

!         REDUCE PATH END POINTS ABOVE ZMAX TO ZMAX
          IF(HMIN.GE.ZMAX)THEN
              WRITE(IPR,'(/2A,/(19X,A,F11.5,A))')' Error in GEOINP: ',  &
     &          ' The entire path lies above ZMAX, the top of the'//    &
     &          ' atmospheric profile!','ZMAX =',ZMAX,' KM','H1ALT =',  &
     &          H1ALT,' KM','H2ALT =',H2ALT,' KM','HMIN =',HMIN,' KM'
              IERROR=1
          ELSE
              IF(H1ALT.GT.ZMAX .OR. H2ALT.GT.ZMAX)                      &
     &          CALL REDUCE(REARTH,ZMAX,H1ALT,H2ALT,OBSZEN,BCKZEN)
              IF(.NOT.LJMASS .AND. NPR.LT.1)                            &
     &          WRITE(IPR,'(//A,/5(/10X,A,F11.5,A),                     &
     &          /10X,A,I11)')' SLANT PATH PARAMETERS IN STANDARD FORM', &
     &          'H1ALT  =',H1ALT ,' KM' ,'H2ALT  =',H2ALT ,' KM' ,      &
     &          'OBSZEN =',OBSZEN,' DEG','BCKZEN =',BCKZEN,' DEG',      &
     &          'HMIN   =',HMIN  ,' KM' ,'LENN   =',LENN
          ENDIF
      ENDIF
      RETURN
      END
