      LOGICAL FUNCTION CNVPHI(ITYPE,                                    &
     &  REARTH,ZMAX,H1ALT,H2ALT,ANGLE,HRANGE,BETA,LENN,PHI)

!     THIS ROUTINE CONVERTS TARGET-BASED LINE-OF-SIGHT
!     INPUTS INTO MODTRAN STANDARD OBSERVER-BASED INPUTS.
!     IF UNSUCCESSFUL, CNVPHI IS SET TO .FALSE.  THREE
!     TARGET-BASED LINE-OF-SIGHT INPUT OPTIONS ARE AVAILABLE:

!       CASE     INPUTS
!       ----     --------------------
!        2E      H2ALT, PHI, H1ALT and LENN
!        2F      H2ALT, PHI and HRANGE
!        3C      H2ALT, PHI and H1ALT=SPACE

!     BOTH CASE 2E AND CASE 3C ARE CONVERTED TO CASE
!     2A (H1ALT, H2ALT, ANGLE AND LENN), AND CASE 2F IS
!     CONVERTED TO CASE 2C (H1ALT, H2ALT AND HRANGE)
      IMPLICIT NONE

!     ARGUMENTS:
!       ITYPE    PATH TYPE LABEL (1=HORIZONTAL, 2=SLANT
!                BETWEEN ALTITUDES and 3=SLANT TO SPACE).
!       REARTH   RADIUS OF THE EARTH [KM].
!       ZMAX     MAXIMUM ATMOSPHERIC PROFILE ALTITUDE [KM].
!       H1ALT    OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    TARGET (FINAL) ALTITUDE [KM].
!       ANGLE    ZENITH ANGLE AT OBSERVER TOWARDS TARGET [DEG].
!       HRANGE   OBSERVER TO TARGET SLANT HRANGE [KM].
!       BETA     OBSERVER TO TARGET EARTH CENTER ANGLE [DEG].
!       RE       RADIUS OF THE EARTH [KM].
!       LENN     PATH LENGTH SWITCH (0=SHORT PATH and
!                  1=LONG PATH THROUGH TANGENT HEIGHT).
!       PHI      ZENITH ANGLE AT TARGER TOWARDS OBSERVER [DEG].
      INTEGER ITYPE,LENN
      DOUBLE PRECISION REARTH,ZMAX,H1ALT,H2ALT,ANGLE,HRANGE,BETA,PHI

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       IERROR   ERROR FLAG
!       LENNSV   SAVED VALUE OF INPUT LENGTH SWITCH, LENN.
!       H1TEST   COPY OF H1ALT [KM].
!       HMIN     MINIMUM PATH ALTITUDE [KM].
      INTEGER IERROR,LENNSV
      DOUBLE PRECISION H1TEST,HMIN

!     INITIALIZE CNVPHI
      CNVPHI=.FALSE.
      IERROR=0

!     BRANCH BASED ON CASE.
      IF(ITYPE.EQ.3)THEN

!         CASE 3C:  H2ALT, PHI and H1ALT=SPACE.
          H1ALT=ZMAX
          LENN=0
          IF(PHI.GT.90.D0)LENN=1
          CALL FINDMN(REARTH,H2ALT,PHI,H1ALT,LENN,0,                    &
     &      HMIN,ANGLE,IERROR,.TRUE.)
          IF(IERROR.NE.0)RETURN
          ITYPE=2
          HRANGE=0.D0
          IF(.NOT.LJMASS)WRITE(IPR,                                     &
     &      '(2(/A,T50,A),3(/A,F12.5,A,T50,A,F12.5,A),/T50,A,I6)')      &
     &      ' CASE 3C                    CONVERTED TO','CASE 2A',       &
     &      ' -------',                      '-------',                 &
     &      ' H2ALT',H2ALT,' KM',            'H1ALT',H1ALT,' KM',       &
     &      ' PHI  ',PHI,  ' DEG',           'ANGLE',ANGLE,' DEG',      &
     &      ' H1ALT',ZMAX, ' KM',            'H2ALT',H2ALT,' KM',       &
     &                                       'LENN ',LENN
      ELSEIF(HRANGE.GT.0.D0)THEN

!         CASE 2F:  H2ALT, PHI and HRANGE.
          CALL NEWH2(REARTH,ZMAX,H2ALT,H1ALT,PHI,HRANGE,                &
     &      BETA,LENN,HMIN,ANGLE)
          ANGLE=0.D0
          IF(.NOT.LJMASS)WRITE(IPR,                                     &
     &      '(2(/A,T50,A),/(A,F12.5,A,T50,A,F12.5,A))')                 &
     &      ' CASE 2F                    CONVERTED TO','CASE 2C',       &
     &      ' --------',                      '-------',                &
     &      ' H2ALT ',H2ALT,' KM',            'H1ALT ',H1ALT, ' KM',    &
     &      ' PHI   ',PHI,  ' DEG',           'HRANGE',HRANGE,' KM',    &
     &      ' HRANGE',HRANGE,' KM',           'H2ALT ',H2ALT, ' KM'
      ELSE

!         CASE 2E:  H2ALT, PHI, H1ALT and LENN.
          ANGLE=0.D0
          LENNSV=LENN
          IF(H2ALT.GE.H1ALT .AND. PHI.LE.90.D0)THEN
              WRITE(IPR,'(/2A)')' Error in CNVPHI: ',                   &
     &          ' Upward path begins above final altitude.'
              RETURN
          ELSEIF(H2ALT.LE.GNDALT .AND. PHI.GT.90.D0)THEN
              WRITE(IPR,'(/2A)')' Error in CNVPHI: ',                   &
     &          ' Downward path begins at or below the Earth surface.'
              RETURN
          ENDIF
          IF(.NOT.LJMASS .AND. H1ALT.LT.H2ALT .AND. PHI.GT.90.D0)       &
     &      WRITE(IPR,'(/2A,I6)')                                       &
     &      ' Either a short path (LENN=0) or a long path through',     &
     &      ' a tangent height (LENN=1) is possible:  lenn =',LENN
          H1TEST=H1ALT
          CALL FINDMN(REARTH,H2ALT,PHI,H1TEST,                          &
     &      LENN,0,HMIN,ANGLE,IERROR,.TRUE.)
          IF(IERROR.NE.0)RETURN
          IF(H1ALT.NE.H1TEST)THEN
              WRITE(IPR,'(/2A)')' Error in CNVPHI:  Slant path',        &
     &          ' intersects the Earth and cannot reach H1ALT.'
              RETURN
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,                                     &
     &      '(2(/A,T50,A),3(/A,F12.5,A,T50,A,F12.5,A),/A,I6,T50,A,I6)') &
     &      ' CASE 2E                    CONVERTED TO','CASE 2A',       &
     &      ' -------',                      '-------',                 &
     &      ' H2ALT',H2ALT,' KM',            'H1ALT',H1ALT,' KM',       &
     &      ' PHI  ',PHI,  ' DEG',           'ANGLE',ANGLE,' DEG',      &
     &      ' H1ALT',H1ALT,' KM',            'H2ALT',H2ALT,' KM',       &
     &      ' LENN ',LENNSV,                 'LENN ',LENN
      ENDIF
      BETA=0.D0

!     RETURN TO DRIVER WITH SUCCESS.
      CNVPHI=.TRUE.
      RETURN
      END
