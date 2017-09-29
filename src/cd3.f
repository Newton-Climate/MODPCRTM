      BLOCK DATA CD3DAT

!     COMMONS:

!     /JM3A2/
!       IPRMSV   SAVED VALUE OF IPARM INPUT.
!       PARM1    SENSOR LATITUDE IF IPARM=0 OR =1 [DEG N OF EQUATOR].
!                TARGET LATITUDE IF IPARM=10 OR =11 [DEG N OF EQUATOR].
!                SOLAR REL AZIMUTH AT SENSOR IF IPARM=2 [DEG E OF N].
!                SOLAR REL AZIMUTH AT TARGET IF IPARM=12 [DEG E OF N].
!       PARM2    SENSOR LONGITUDE IF IPARM=0 OR =1 [DEG W OF GRNWCH].
!                TARGET LONGITUDE IF IPARM=10 OR =11 [DEG W OF GRNWCH].
!                SOLAR ZENITH AT SENSOR IF IPARM IS 2 [DEG].
!                SOLAR ZENITH AT TARGET IF IPARM IS 12 [DEG].
!       SRCLAT   SOURCE (SUN OR MOON) LATITUDE [DEG NORTH OF EQUATOR].
!       SRCLON   SOURCE (SUN OR MOON) LONGITUDE [DEG WEST OF GREENWICH].
!       GMTIME   GREENWICH MEAN TIME [DECIMAL HOURS].
!       TRUEAZ   PATH TRUE AZIMUTH ANGLE [DEG EAST OF NORTH].
      INTEGER IPRMSV
      DOUBLE PRECISION PARM1,PARM2,SRCLAT,SRCLON,GMTIME,TRUEAZ
      COMMON/JM3A2/PARM1,PARM2,SRCLAT,SRCLON,GMTIME,TRUEAZ,IPRMSV
      DATA IPRMSV/0/
      END

      SUBROUTINE CD3(ILOS,REFPTH)

!     PROCESS CARD3 INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       ILOS     LOOP INDEX FOR OBSERVER ZENITH ANGLES.
!       REFPTH   NAME OF REFRACTIVE PATH INPUT FILE.
      INTEGER ILOS
      CHARACTER REFPTH*(NAMLEN)

!     COMMONS:
      INCLUDE 'IFIL.h'

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

!     /CARD3/
!       H1ALT    OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    FINAL (TANGENT FOR LIMB PATH) ALTITUDE [KM].
!       OBSZEN   OBSERVER ZENITH ANGLE (H1ALT TO H2ALT) [DEG].
!       HRANGE   DISTANCE FROM H1ALT TO H2ALT [KM].
!       BETA     EARTH CENTER ANGLE BETWEEN H1ALT AND H2ALT [DEG].
!       HMIN     PATH MINIMUM ALTITUDE [KM].
!       HMAX     PATH MAXIMUM ALTITUDE [KM].
!       CKRANG   MAXIMUM PATH RANGE FOR K-DISTRIBUTION OUTPUT
!                (=0. FOR TOTAL PATH ONLY; <0. FOR ALL RANGES).
!       BCKZEN   ZENITH ANGLE FOR BACKWARD (H2ALT TO H1ALT) PATH [DEG].
!       REARTH   EARTH RADIUS [KM].
!       ANGLEM   LUNAR PHASE ANGLE [0 TO 180 DEG].
!       LENN     PATH LENGTH SWITCH (0=SHORT, 1=LONG).
!       IDAY     DAY OF YEAR [0-366, DEFAULT (0) IS DAY 91].
!       ISOURC   SOURCE FLAG [0=SUN AND 1=MOON].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
      INTEGER LENN,IDAY,ISOURC,NLOS
      DOUBLE PRECISION H1ALT,H2ALT,OBSZEN,HRANGE,BETA,REARTH,           &
     &  HMIN,HMAX,CKRANG,BCKZEN,ANGLEM
      COMMON/CARD3/H1ALT(MLOS),H2ALT(MLOS),OBSZEN(MLOS),HRANGE(MLOS),   &
     &  BETA(MLOS),HMIN(MLOS),HMAX(MLOS),CKRANG(MLOS),BCKZEN(MLOS),     &
     &  REARTH,ANGLEM,LENN(MLOS),IDAY,ISOURC,NLOS

!     /JM3A2/
!       IPRMSV   SAVED VALUE OF IPARM INPUT.
!       PARM1    SENSOR LATITUDE IF IPARM=0 OR =1 [DEG N OF EQUATOR].
!                TARGET LATITUDE IF IPARM=10 OR =11 [DEG N OF EQUATOR].
!                SOLAR REL AZIMUTH AT SENSOR IF IPARM=2 [DEG E OF N].
!                SOLAR REL AZIMUTH AT TARGET IF IPARM=12 [DEG E OF N].
!       PARM2    SENSOR LONGITUDE IF IPARM=0 OR =1 [DEG W OF GRNWCH].
!                TARGET LONGITUDE IF IPARM=10 OR =11 [DEG W OF GRNWCH].
!                SOLAR ZENITH AT SENSOR IF IPARM IS 2 [DEG].
!                SOLAR ZENITH AT TARGET IF IPARM IS 12 [DEG].
!       SRCLAT   SOURCE (SUN OR MOON) LATITUDE [DEG NORTH OF EQUATOR].
!       SRCLON   SOURCE (SUN OR MOON) LONGITUDE [DEG WEST OF GREENWICH].
!       GMTIME   GREENWICH MEAN TIME [DECIMAL HOURS].
!       TRUEAZ   PATH TRUE AZIMUTH ANGLE [DEG EAST OF NORTH].
      INTEGER IPRMSV
      DOUBLE PRECISION PARM1,PARM2,SRCLAT,SRCLON,GMTIME,TRUEAZ
      COMMON/JM3A2/PARM1,PARM2,SRCLAT,SRCLON,GMTIME,TRUEAZ,IPRMSV

!     LOCAL VARIABLES:
!       IOS      STATUS OF READ STATEMENT.
!       RAD_E    RADIUS OF THE EARTH.
      INTEGER IOS
      DOUBLE PRECISION RAD_E


!     READ CARD3:
      IF(ITYPE.EQ.4)THEN

!         *.pth FILE DEFINES LINE-OF-SIGHT:  ONLY READ RAD_E & CKNAME.
          READ(IRD,'(30X,2(20X,F10.0))')RAD_E,CKRANG(1)
          WRITE(IPR,'(/A,30X,2(20X,F10.5))')                            &
     &      ' CARD 3  *****',RAD_E,CKRANG(1)
          IDAY=-99
          ISOURC=-99
          ANGLEM=-99.D0
      ELSEIF(LJMASS)THEN
          CALL INITCD('CARD3')
          IF(IEMSCT.EQ.3)THEN

!             CARD 3 FOR DIRECTLY TRANSMITTED SOLAR IRRADIANCE:
              ITYPE=3
              HRANGE(1)=0.D0
              BETA(1)=0.D0
              LENN(1)=0
              BCKZEN(1)=0.D0
          ELSE
              IDAY=-99
              ISOURC=-99
              ANGLEM=-99.D0
          ENDIF
      ELSEIF(IEMSCT.EQ.3)THEN

!         CARD 3 FOR DIRECTLY TRANSMITTED SOLAR IRRADIANCE (IEMSCT=3):
          READ(IRD,'(2(3F10.0,I5,5X,F10.0,I5))')                        &
     &      H1ALT(1),H2ALT(1),OBSZEN(1),IDAY,RAD_E,ISOURC,ANGLEM
          WRITE(IPR,'(/A,2(3F10.5,I5,F15.5,I5))')' CARD 3   *****',     &
     &      H1ALT(1),H2ALT(1),OBSZEN(1),IDAY,RAD_E,ISOURC,ANGLEM
          IF(ISOURC.EQ.1)WRITE(IPR,'(/A)')'   LUNAR SOURCE ONLY'
          ITYPE=3
          HRANGE(1)=0.D0
          BETA(1)=0.D0
          LENN(1)=0
          BCKZEN(1)=0.D0
      ELSE
          READ(IRD,'(6F10.0,I5,5X,2F10.0)',IOSTAT=IOS)                  &
     &      H1ALT(ILOS),H2ALT(ILOS),OBSZEN(ILOS),HRANGE(ILOS),          &
     &      BETA(ILOS),RAD_E,LENN(ILOS),BCKZEN(ILOS),CKRANG(ILOS)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/A,//A,5F10.5,F12.5,I3,5X,2F10.5)')           &
     &          ' Error in CD3:  Read of CARD3 was unsuccessful!',      &
     &          ' CARD 3  *****',H1ALT(ILOS),H2ALT(ILOS),               &
     &          OBSZEN(ILOS),HRANGE(ILOS),BETA(ILOS),RAD_E,             &
     &          LENN(ILOS),BCKZEN(ILOS),CKRANG(ILOS)
              STOP 'Read of CARD3 was unsuccessful!'
          ENDIF
          WRITE(IPR,'(/A,5F10.5,F12.5,I3,5X,2F10.5)')' CARD 3  *****',  &
     &      H1ALT(ILOS),H2ALT(ILOS),OBSZEN(ILOS),HRANGE(ILOS),          &
     &      BETA(ILOS),RAD_E,LENN(ILOS),BCKZEN(ILOS),CKRANG(ILOS)
          IF(ILOS.EQ.1)THEN
              IDAY=-99
              ISOURC=-99
              ANGLEM=-99.D0
          ELSEIF(IPRMSV.EQ.2 .OR. IPRMSV.EQ.12)THEN

!             CONVERT RELATIVE AZIMUTH TO TRUE AZIMUTH, HAVING PLACED
!             THE SUN WEST OF H1 (IF IPRMSV=2) OR H2 (IF IPRMSV=12):
              TRUEAZ=270-RAD_E
          ELSE

!             STORE INPUT VALUE OF TRUE AZIMUTH:
              TRUEAZ=RAD_E
          ENDIF
      ENDIF

!     REARTH IS THE RADIUS OF THE EARTH USED BY MODTRAN
      IF(MODEL.NE.8 .AND. ILOS.EQ.1)REARTH=RAD_E
      IF(REARTH.LE.0.D0)THEN
          IF(MODEL.EQ.1)THEN
              REARTH=6378.39D0
          ELSEIF(MODEL.EQ.4)THEN
              REARTH=6356.91D0
          ELSEIF(MODEL.EQ.5)THEN
              REARTH=6356.91D0
          ELSEIF(MODEL.NE.0)THEN
              REARTH=6371.23D0
          ENDIF
      ENDIF

!     READ IN USER DEFINED PATH DATA:
      IF(ITYPE.EQ.4)CALL RDPATH(REFPTH)

!     UNLESS LJMASS, WRITE GEOMETRY INPUT INFORMATION:
      IF(LJMASS)RETURN
      WRITE(IPR,'(/F12.2,A)')REARTH,'  RADIUS OF THE EARTH [KM].'
      IF(ITYPE.EQ.1)THEN
          WRITE(IPR,'(/A,/(10X,A,F11.5,A))')' HORIZONTAL PATH',         &
     &      'ALTITUDE =',H1ALT(1),' KM','HRANGE    =',HRANGE(1),' KM'
      ELSEIF(ITYPE.EQ.2)THEN
          WRITE(IPR,FMT='(/A,I3,A,7(/10X,A,F11.5,1X,A),/10X,A,I7)')     &
     &      ' SLANT PATH No.',ILOS,', H1ALT TO H2ALT',                  &
     &      'H1ALT  =', H1ALT(ILOS), 'KM','H2ALT  =', H2ALT(ILOS), 'KM',&
     &      'OBSZEN =',OBSZEN(ILOS),'DEG','HRANGE =',HRANGE(ILOS), 'KM',&
     &      'BETA   =',  BETA(ILOS),'DEG','BCKZEN =',BCKZEN(ILOS),'DEG',&
     &      'CKRANG =',CKRANG(ILOS), 'KM','LENN   =',  LENN(ILOS)
      ELSEIF(ITYPE.EQ.4)THEN
          WRITE(IPR,'(/A,6(/10X,A,F11.5,1X,A),/10X,A,I7)')              &
     &      ' USER-DEFINED PATH WITH', 'H1ALT  =',H1ALT(1) ,'KM',       &
     &      'H2ALT  =',H2ALT(1) ,'KM' ,'OBSZEN =',OBSZEN(1),'DEG',      &
     &      'HRANGE =',HRANGE(1),'KM' ,'BETA   =',BETA(1)  ,'DEG',      &
     &      'BCKZEN =',BCKZEN(1),'DEG','LENN   =',LENN(1)
      ELSEIF(BCKZEN(1).LE.0.D0)THEN
          WRITE(IPR,'(/A,/(10X,A,F11.5,1X,A))')                         &
     &      ' SLANT PATH TO SPACE (OR GROUND)','H1ALT  =',H1ALT(1),'KM',&
     &      'HMIN   =',H2ALT(1),'KM',    'OBSZEN =',OBSZEN(1),'DEG'
      ELSE
          WRITE(IPR,'(/A,/(10X,A,F11.5,1X,A))')                         &
     &      ' SLANT PATH TO SPACE (OR GROUND)',                         &
     &      'H2ALT  =',H2ALT(1),'KM','BCKZEN =',BCKZEN(1),'DEG'
      ENDIF
!     RETURN TO DRIVER:
      RETURN
      END

      SUBROUTINE RDPATH(REFPTH)

!     RDPATH READS IN REFRACTIVE PATH DATA.
      IMPLICIT NONE

!     PARAMETERS:
!       ALTTOL   ALTITUDE TOLERANCE OF 0.1 MILLIMETERS [KM].
      INCLUDE 'PARAMS.h'
      DOUBLE PRECISION ALTTOL
!old  PARAMETER(ALTTOL=2.D-6)
      PARAMETER(ALTTOL=1.D-7)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /PATH/
!       PTHCOS   COSINE OF PATH ZENITH AT PATH BOUNDARIES.
!       PTHZEN   PATH ZENITH AT PATH BOUNDARIES [DEG].
!       PTHECA   SENSOR TO PATH EARTH CENTER ANGLE [DEG].
!       PTHALT   ALTITUDES AT PATH BOUNDARIES [KM].
!       PTH_MS   ALTITUDES AT PATH BOUNDARIES FOR THE MS PATH.
!       PTHSEG   PATH SEGMENT LENGTH [KM].
!       PTHRNG   SENSOR TO PATH BOUNDARY RANGE [KM].
!       JMAX     NUMBER OF DISTINCT LOS PATH SEGMENT ENDPOINT ALTITUDES.
!       IKHMIN   PATH BOUNDARY INDEX OF PATH MINIMUM ALTITUDE.
!       IKHMAX   PATH BOUNDARY INDEX OF PATH MAXIMUM ALTITUDE.
!       IKOUT    NUMBER OF PATH BOUNDARIES K DATA IS OUTPUT.
!       NTKDIS   RECORD NUMBER FOR K-DISTRIBUTION TRANSMITTANCE FILE.
!       NRKDIS   RECORD NUMBER FOR K-DISTRIBUTION RADIANCE FILE.
!       MAPPTH   MAPPING FROM PATH SEGMENT MIDPOINT TO VERTICAL LAYER.
!       IPTHHT   ALTITUDES (HEIGHTS) AT PATH BOUNDARIES [M].
!       LOWALT   VERTICAL LAYER BOUNDARY INDEX AT OR JUST BELOW PTHALT.
!       FACALT   ALTITUDE INTERPOLATION FRACTION FOR PTHALT
!       PATH_T   TEMPERATURE AT PATH BOUNDARIES [K].
!       PATH_P   PRESSURE AT PATH BOUNDARIES [ATM].
!       PTHRH    RELATIVE HUMIDITY AT PATH BOUNDARIES [K].
!       LSSGEO   LOGICAL FLAG, .TRUE. FOR SOLAR PATHS.
!       LTANMX   LOGICAL FLAG, .TRUE. IF PATH HAS A TANGENT MAXIMUM.
      DOUBLE PRECISION PTHCOS,PTHZEN,PTHECA,PTHALT,PTH_MS,PTHSEG,PTHRNG
      INTEGER JMAX,IKHMIN,IKHMAX,IKOUT,NTKDIS,NRKDIS,MAPPTH,            &
     &  IPTHHT,LOWALT
      REAL FACALT,PATH_T,PATH_P,PTHRH
      LOGICAL LSSGEO,LTANMX
      COMMON/PATH/PTHCOS(0:LAYTWO),PTHZEN(0:LAYTWO),PTHECA(0:LAYTWO),   &
     &  PTHALT(0:LAYTWO,1:MLOS),PTH_MS(0:LAYDIM),PTHSEG(LAYTWO),        &
     &  PTHRNG(0:LAYTWO,1:MLOS),JMAX,IKHMIN(MLOS),IKHMAX(MLOS),         &
     &  IKOUT(MLOS),MAPPTH(LAYTWO,1:MLOS),IPTHHT(0:LAYTWO),NTKDIS,      &
     &  NRKDIS,LOWALT(0:LAYTWO,1:MLOS),FACALT(0:LAYTWO,1:MLOS),         &
     &  PATH_T(0:LAYTWO,1:MLOS),PATH_P(0:LAYTWO,1:MLOS),                &
     &  PTHRH(0:LAYTWO,1:MLOS),LSSGEO,LTANMX

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

!     /CARD3/
!       H1ALT    OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    FINAL (TANGENT FOR LIMB PATH) ALTITUDE [KM].
!       OBSZEN   OBSERVER ZENITH ANGLE (H1ALT TO H2ALT) [DEG].
!       HRANGE   DISTANCE FROM H1ALT TO H2ALT [KM].
!       BETA     EARTH CENTER ANGLE BETWEEN H1ALT AND H2ALT [DEG].
!       REARTH   RADIUS OF THE EARTH [KM].
!       HMIN     PATH MINIMUM ALTITUDE [KM].
!       HMAX     PATH MAXIMUM ALTITUDE [KM].
!       CKRANG   MAXIMUM PATH RANGE FOR K-DISTRIBUTION OUTPUT
!                (=0. FOR TOTAL PATH ONLY; <0. FOR ALL RANGES).
!       BCKZEN   ZENITH ANGLE FOR BACKWARD (H2ALT TO H1ALT) PATH [DEG].
!       REARTH   EARTH RADIUS [KM].
!       ANGLEM   LUNAR PHASE ANGLE [0 TO 180 DEG].
!       LENN     PATH LENGTH SWITCH (0=SHORT, 1=LONG).
!       IDAY     DAY OF YEAR [0-366, DEFAULT (0) IS DAY 91].
!       ISOURC   SOURCE FLAG [0=SUN AND 1=MOON].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
      INTEGER LENN,IDAY,ISOURC,NLOS
      DOUBLE PRECISION H1ALT,H2ALT,OBSZEN,HRANGE,BETA,REARTH,           &
     &  HMIN,HMAX,CKRANG,BCKZEN,ANGLEM
      COMMON/CARD3/H1ALT(MLOS),H2ALT(MLOS),OBSZEN(MLOS),HRANGE(MLOS),   &
     &  BETA(MLOS),HMIN(MLOS),HMAX(MLOS),CKRANG(MLOS),BCKZEN(MLOS),     &
     &  REARTH,ANGLEM,LENN(MLOS),IDAY,ISOURC,NLOS

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

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     INPUT ARGUMENTS:
!       REFPTH   NAME OF REFRACTIVE PATH INPUT FILE.
      CHARACTER REFPTH*(*)

!     FUNCTIONS:
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER LENSTR

!     LOCAL VARIABLES:
!       IOS      STATUS OF READ STATEMENT.
!       IBANG    LOCATION IN CHARACTER STRING OF BANG (EXCLAMATION).
!       JCHAR    LOOP VARIABLE FOR CHARACTER POSITION.
!       ICOMMA   LOCATION IN CHARACTER STRING OF COMMA.
!       IEND     END POSITION OF NON-BLANK STRING.
!       ISTART   START POSITION OF NEXT INPUT VARIABLE.
!       LAY      VERTICAL LAYER INDEX.
!       LAYMIN   INDEX OF LAYER CONTAINING PATH MINIMUM ALTITUDE.
!       LAYMAX   INDEX OF LAYER CONTAINING PATH MAXIMUM ALTITUDE.
!       IK       PATH BOUNDARY INDEX.
!       ESD2CA   EARTH SURFACE DIST TO CENTER ANGLE CONVERSION FACTOR.
!       ESD      PATH BOUNDARY EARTH SURFACE DISTANCE [M].
!       ALT      PATH BOUNDARY ALTITUDE [M].
!       SEG      PATH BOUNDARY SEGMENT LENGTH [M].
!       LOPEN    LOGICAL FLAG, .TRUE. IF FILE IS OPEN.
!       INLINE   INPUT LINE READ IN AS A CHARACTER STRING.
      INTEGER IOS,IBANG,JCHAR,ICOMMA,IEND,ISTART,LAY,LAYMIN,LAYMAX,IK
      DOUBLE PRECISION ESD2CA,ESD,ALT,SEG
      LOGICAL LOPEN
      CHARACTER INLINE*(132)

!     DATA:
!       FRMT     INPUT VARIABLE FORMAT.
      CHARACTER FRMT*(8)
      DATA FRMT/'(F000.0)'/

!     OPEN REFRACTIVE PATH DATA FILE:
      INQUIRE(LPATH,OPENED=LOPEN)
      IF(.NOT.LOPEN)CALL OPNFL(LPATH,0,REFPTH,'OLD','FORMATTED','CD1')

!     INITIALIZE LINE-OF-SIGHT SEGMENT COUNTER, EARTH SURFACE DISTANCE
!     TO CENTER ANGLE CONVERSION FACTOR AND PATH RANGE:
      NSEG(1)=-1
      ESD2CA=DPDEG/(1000*REARTH)
      PTHRNG(0,1)=0.D0

!     READ AND WRITE REFRACTIVE PATH DATA:
      WRITE(IPR,'(/(A))')' REFRACTIVE PATH INPUT DATA:',                &
     &                   ' ---------------------------'
   10 CONTINUE
      READ(LPATH,'(A132)',IOSTAT=IOS)INLINE

!     PROCESS INLINE IF READ WAS SUCCESSFUL:
      IF(IOS.EQ.0)THEN
          WRITE(IPR,'(A)')INLINE

!         REMOVE COMMENTS:
          IBANG=INDEX(INLINE,'!')
          IF(IBANG.GT.0)THEN
              DO JCHAR=IBANG,132
                  INLINE(JCHAR:JCHAR)=' '
              ENDDO
          ENDIF

!         REPLACE ALL COMMAS WITH SPACES:
   20     CONTINUE
          ICOMMA=INDEX(INLINE,',')
          IF(ICOMMA.GT.0)THEN
              INLINE(ICOMMA:ICOMMA)=' '
              GOTO 20
          ENDIF

!         IF INPUT LINE IS BLANK, GO TO NEXT LINE:
          IF(LENSTR(INLINE).LE.0)GOTO 10

!         CHECK THAT STRING BEGINS WITH A DIGIT, SIGN OR DECIMAL; OTHER-
!         WISE, ASSUME IT IS THE BEGINNING OF THE NEXT REFRACTIVE PATH:
          IF(INDEX('+-.0123456789',INLINE(1:1)).GT.0)THEN

!             DETERMINE END POSITION OF FIRST INPUT:
              NSEG(1)=NSEG(1)+1
              IEND=INDEX(INLINE,' ')-1
              IF(IEND.LE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the EARTH SURFACE DISTANCE'     &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             READ EARTH SURFACE DISTANCE:
              WRITE(FRMT(3:5),'(I3.3)')IEND
              READ(INLINE(1:IEND),FRMT,IOSTAT=IOS)ESD
              IF(IOS.NE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the EARTH SURFACE DISTANCE'     &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             DETERMINE START POSITION OF SECOND INPUT:
              ISTART=IEND+1
              IF(LENSTR(INLINE(ISTART:132)).LE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH BOUNDARY ALTITUDE'     &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             DETERMINE END POSITION OF SECOND INPUT:
              IEND=INDEX(INLINE,' ')-1
              IF(IEND.LT.ISTART)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH BOUNDARY ALTITUDE'     &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             READ PATH BOUNDARY ALTITUDE:
              WRITE(FRMT(3:5),'(I3.3)')IEND-ISTART+1
              READ(INLINE(ISTART:IEND),FRMT,IOSTAT=IOS)ALT
              IF(IOS.NE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH BOUNDARY ALTITUDE'     &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             DETERMINE START POSITION OF THIRD INPUT:
              ISTART=IEND+1
              IF(LENSTR(INLINE(ISTART:132)).LE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH BOUNDARY ZENITH ANGLE' &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             DETERMINE END POSITION OF THIRD INPUT:
              IEND=INDEX(INLINE,' ')-1
              IF(IEND.LT.ISTART)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH BOUNDARY ZENITH ANGLE' &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             READ PATH BOUNDARY ZENITH ANGLE:
              WRITE(FRMT(3:5),'(I3.3)')IEND-ISTART+1
              READ(INLINE(ISTART:IEND),FRMT,IOSTAT=IOS)PTHZEN(NSEG(1))
              IF(IOS.NE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH BOUNDARY ZENITH ANGLE' &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             DETERMINE START POSITION OF FOURTH INPUT:
              ISTART=IEND+1
              IF(LENSTR(INLINE(ISTART:132)).LE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH SEGMENT LENGTH ending' &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             DETERMINE END POSITION OF FOURTH INPUT:
              IEND=INDEX(INLINE,' ')-1
              IF(IEND.LT.ISTART)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH SEGMENT LENGTH ending' &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             READ PATH SEGMENT LENGTH:
              WRITE(FRMT(3:5),'(I3.3)')IEND-ISTART+1
              READ(INLINE(ISTART:IEND),FRMT,IOSTAT=IOS)SEG
              IF(IOS.NE.0)THEN

!                 ERROR READING FILE:
                  WRITE(IPR,'(/(2A,I4,A))')' Error in RDPATH:  Unable', &
     &              ' to parse and read the PATH SEGMENT LENGTH ending' &
     &              //' at path boundary',NSEG(1),' in file',REFPTH
                  STOP 'Error reading refractive path input file!'
              ENDIF

!             PROCESS/STORE DATA (CONVERT M TO KM):
              PTHECA(NSEG(1))=ESD2CA*ESD
              PTHALT(NSEG(1),1)=ALT/1000
              PTHCOS(NSEG(1))=COS(PTHZEN(NSEG(1))/DPDEG)
              IF(NSEG(1).EQ.0)THEN

!                 INITIAL PATH ALTITUDE:
                  OBSZEN(1)=PTHZEN(NSEG(1))
                  IKHMIN(1)=0
                  HMIN(1)=PTHALT(0,1)
                  IKHMAX(1)=0
                  HMAX(1)=PTHALT(0,1)
              ELSE
                  PTHSEG(NSEG(1))=SEG/1000
                  PTHRNG(NSEG(1),1)=PTHRNG(NSEG(1)-1,1)+PTHSEG(NSEG(1))
                  IF(PTHALT(NSEG(1),1).LT.HMIN(1))THEN
                      IKHMIN(1)=NSEG(1)
                      HMIN(1)=PTHALT(IKHMIN(1),1)
                  ELSEIF(PTHALT(NSEG(1),1).GT.HMAX(1))THEN
                      IKHMAX(1)=NSEG(1)
                      HMAX(1)=PTHALT(IKHMAX(1),1)
                  ENDIF
              ENDIF

!             GET NEXT DATA LINE:
              GOTO 10
          ENDIF
      ENDIF

!     PATH MUST HAVE AT LEAST ONE SEGMENT:
      IF(NSEG(1).LE.0)THEN

!         PATH DEFINITION ERROR:
          WRITE(IPR,'(/(A))')' Error in RDPATH:  Less than'             &
     &      //' 2 lines of data read from path file',REFPTH
          STOP 'Error processing refractive path input file!'
      ENDIF

!     PATH MUST BE ENCAPSULATED BY DEFINED ATMOSPHERE:
      H1ALT(1)=PTHALT(0,1)
      H2ALT(1)=PTHALT(NSEG(1),1)
      IF(HMAX(1).GT.ZM(ML) .OR. HMIN(1).LT.GNDALT)THEN

!         PATH EXTENDS OUTSIDE OF DEFINED ATMOSPHERE:
          WRITE(IPR,'(/A,/(F20.5,A))')                                  &
     &      ' Error in RDPATH:  Path extends outside of atmosphere',    &
     &      GNDALT,'  Ground altitude [km]',                            &
     &      ZM(ML),'  Top of atmosphere [km]',                          &
     &      H1ALT(1),'  Sensor/Observer altitude [km]',                 &
     &      H2ALT(1),'  Path final altitude [km]',                      &
     &      HMIN(1),'  Path minimum altitude [km]',                     &
     &      HMAX(1),'  Path maximum altitude [km]'
          STOP 'Error in RDPATH:  Path extends outside of atmosphere.'
      ELSEIF((2*IKHMIN(1).EQ.NSEG(1) .OR. 2*IKHMAX(1).EQ.NSEG(1))       &
     &  .AND. ABS(H1ALT(1)-H2ALT(1)).GT.ALTTOL)THEN

!         SINCE THE NUMBER OF PATH SEGMENTS FROM H1ALT TO
!         THE PATH TANGENT EQUALS THE NUMBER FROM THE PATH
!         TANGENT TO H2ALT, H1ALT AND H2ALT MUST BE EQUAL:
          WRITE(IPR,'(/A,F15.8,A,/18X,A,F15.8,A)')                      &
     &      ' Error in RDPATH:  The number of path segments'//          &
     &      ' from the Sensor/Observer (H1ALT =',H1ALT(1),'km) to'//    &
     &      ' the path tangent',' height equals the number from'//      &
     &      ' the tangent height to the final altitude (H2ALT =',       &
     &      H2ALT(1),'km), so H1ALT should have equaled H2ALT!'
          STOP 'Altitude mismatch in RDPATH!'
      ENDIF

!     DEFINE TANGENT MAX LOGICAL FLAG AMD BRANCH ACCORDINGLY:
      LTANMX=IKHMAX(1).NE.0 .AND. IKHMAX(1).NE.NSEG(1)
      IF(LTANMX)THEN

!         PATH HAS A TANGENT MAXIMUM.  DETERMINE THE VERTICAL
!         LAYER CONTAINING THE MAXIMUM PATH ALTITUDE:
          DO LAY=2,ML
              IF(HMAX(1).LE.ZM(LAY)+ALTTOL)GOTO 30
          ENDDO
          WRITE(IPR,'(/(A))')                                           &
     &      ' Error in RDPATH:  Maximum path altitude is above'         &
     &      //' the top-of-atmosphere for the path from file',REFPTH
          STOP 'Maximum path altitude above top-of-atmosphere.'
   30     CONTINUE
          LAYMAX=LAY-1
          LAY=LAYMAX

!         VERIFY THAT ALL PATH BOUNDARIES BETWEEN THE MAXIMUM AND FIRST
!         ALTITUDE (H1ALT) COINCIDE WITH ADJOINING ATMOSPHERIC LEVELS:
          DO IK=IKHMAX(1)-1,1,-1
              IF(2*IKHMAX(1)-IK.EQ.NSEG(1))THEN

!                 VERIFY THAT THE NEXT ALTITIUDE IS H2ALT ( < H1ALT ).
                  IF(ABS(PTHALT(IK,1)-H2ALT(1)).GT.ALTTOL)THEN
                      WRITE(IPR,'(/2A,/(19X,A))')' Error in RDPATH: ',  &
     &                  ' The line-of-sight passes through a maximum'   &
     &                  //' tangent','height with H2ALT exceeding'      &
     &                  //' H1ALT, but the H2ALT altitude is not',      &
     &                  'a path boundary along the upward path in'      &
     &                  //' file',REFPTH
                      STOP 'Altitude mismatch in RDPATH!'
                  ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).LE.ALTTOL)THEN

!                     H2ALT COINCIDES WITH AN ATMOSPHERIC LEVEL:
                      IF(LAY.LE.1)THEN

!                         PATH EXTENDS BELOW BOTTOM OF ATMOSPHERE.
                          WRITE(IPR,'(/A)')' Error in RDPATH:  Bottom'//&
     &                      ' of atmosphere reached before end of path.'
                          STOP                                          &
     &                      'Path below bottom of atmosphere in RDPATH!'
                      ENDIF
                      LAY=LAY-1
                  ENDIF
              ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).GT.ALTTOL)THEN

!                 ALTITUDE MISMATCH:
                  WRITE(IPR,'(/A,I3,A,F15.8,A,/(19X,A,I3,A,F15.8,A))')  &
     &              ' Error in RDPATH:  Path altitude',IK,' (',         &
     &              PTHALT(IK,1),' km) does not equal the next lower',  &
     &              'atmospheric level (Level',LAY,' =',ZM(LAY),        &
     &              ' km) in file',REFPTH
                  STOP 'Altitude mismatch in RDPATH!'
              ELSEIF(LAY.LE.1)THEN

!                 PATH EXTENDS BELOW BOTTOM OF ATMOSPHERE.
                  WRITE(IPR,'(/A)')' Error in RDPATH:  Bottom of'       &
     &              //' atmosphere reached before end of path.'
                  STOP                                                  &
     &              'Path extends below bottom of atmosphere in RDPATH!'
              ELSE
                  LAY=LAY-1
              ENDIF
          ENDDO

!         VERIFY THAT FIRST PATH ALTITUDE IS IN NEXT LAYER:
          IF(PTHALT(0,1).LT.ZM(LAY)-ALTTOL)THEN
              WRITE(IPR,'(/A,F15.8,A,/(18X,A,I3,A,F15.8,A))')           &
     &          ' Error in RDPATH:  Path altitude  0 (',H1ALT(1),       &
     &          ' km) is below',' the next lower atmospheric level'     &
     &          //'(Level ',LAY,' =',ZM(LAY),' km) in file',REFPTH
              STOP 'Altitude mismatch in RDPATH!'
          ENDIF

!         VERIFY THAT ALL PATH BOUNDARIES BETWEEN THE MAXIMUM AND LAST
!         ALTITUDE COINCIDE WITH ADJOINING ATMOSPHERIC LEVELS:
          LAY=LAYMAX
          DO IK=IKHMAX(1)+1,NSEG(1)-1
              IF(2*IKHMAX(1)-IK.EQ.0)THEN

!                 VERIFY THAT THE NEXT ALTITIUDE IS H1ALT ( < H2ALT ).
                  IF(ABS(PTHALT(IK,1)-H1ALT(1)).GT.ALTTOL)THEN
                      WRITE(IPR,'(/2A,/(19X,A))')' Error in RDPATH: ',  &
     &                  ' The line-of-sight passes through a maximum'   &
     &                  //' tangent','height with H1ALT exceeding'      &
     &                  //' H2ALT, but the H1ALT altitude is not',      &
     &                  'a path boundary along the downward path in'    &
     &                  //' file',REFPTH
                      STOP 'Altitude mismatch in RDPATH!'
                  ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).LE.ALTTOL)THEN

!                     H1ALT COINCIDES WITH AN ATMOSPHERIC LEVEL:
                      IF(LAY.LE.1)THEN

!                         PATH EXTENDS ABOVE TOP-OF-ATMOSPHERE.
                          WRITE(IPR,'(/A)')' Error in RDPATH:  Bottom'//&
     &                      ' of Atmosphere reached before end of path.'
                          STOP                                          &
     &                      'Path below bottom of atmosphere in RDPATH!'
                      ENDIF
                      LAY=LAY-1
                  ENDIF
              ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).GT.ALTTOL)THEN

!                 ALTITUDE MISMATCH:
                  WRITE(IPR,'(/A,I3,A,F15.8,A,/(19X,A,I3,A,F15.8,A))')  &
     &              ' Error in RDPATH:  Path altitude',IK,' (',         &
     &              PTHALT(IK,1),' km) does not equal the next lower',  &
     &              'atmospheric level (Level',LAY,' =',ZM(LAY),        &
     &              ' km) in file',REFPTH
                  STOP 'Altitude mismatch in RDPATH!'
              ELSEIF(LAY.LE.1)THEN

!                 PATH EXTENDS ABOVE TOP-OF-ATMOSPHERE.
                  WRITE(IPR,'(/A)')' Error in RDPATH:  Bottom of'       &
     &              //' Atmosphere reached before end of path.'
                  STOP                                                  &
     &              'Path extends below bottom of atmosphere in RDPATH!'
              ELSE
                  LAY=LAY-1
              ENDIF
          ENDDO

!         VERIFY THAT LAST PATH ALTITUDE IS IN NEXT LAYER:
          IF(PTHALT(NSEG(1),1).LT.ZM(LAY)-ALTTOL)THEN
              WRITE(IPR,'(/1X,A,I3,A,F15.8,A,/(18X,A,I3,A,F15.8,A))')   &
     &          'Error in RDPATH:  Path altitude',NSEG(1),' (',H2ALT(1),&
     &          ' km) is below',' the next lower atmospheric level'     &
     &          //' (Level',LAY,' =',ZM(LAY),' km) in file',REFPTH
              STOP 'Altitude mismatch in RDPATH!'
          ENDIF
      ELSE

!         NO TANGENT MAXIMUM.  DETERMINE THE VERTICAL
!         LAYER CONTAINING THE MINIMUM PATH ALTITUDE:
          DO LAY=2,ML
              IF(HMIN(1).LT.ZM(LAY)-ALTTOL)GOTO 40
          ENDDO
          WRITE(IPR,'(/(A))')                                           &
     &      ' Error in RDPATH:  Minimum path altitude is at or above'   &
     &      //' the top-of-atmosphere for the path from file',REFPTH
          STOP 'Minimum path altitude at or above top-of-atmosphere.'
   40     CONTINUE
          LAYMIN=LAY-1

!         VERIFY THAT ALL PATH BOUNDARIES BETWEEN THE MINIMUM AND FIRST
!         ALTITUDE (H1ALT) COINCIDE WITH ADJOINING ATMOSPHERIC LEVELS:
          DO IK=IKHMIN(1)-1,1,-1
              IF(2*IKHMIN(1)-IK.EQ.NSEG(1))THEN

!                 VERIFY THAT THE NEXT ALTITIUDE IS H2ALT ( < H1ALT ).
                  IF(ABS(PTHALT(IK,1)-H2ALT(1)).GT.ALTTOL)THEN
                      WRITE(IPR,'(/2A,/(19X,A))')' Error in RDPATH: ',  &
     &                  ' The line-of-sight passes through a minimum'   &
     &                  //' tangent','height with H2ALT being less'     &
     &                  //' than H1ALT, but the H2ALT altitude is not', &
     &                  'a path boundary along the downward path in'    &
     &                  //' file',REFPTH
                      STOP 'Altitude mismatch in RDPATH!'
                  ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).LE.ALTTOL)THEN

!                     H2ALT COINCIDES WITH AN ATMOSPHERIC LEVEL:
                      IF(LAY.GE.ML)THEN

!                         PATH EXTENDS ABOVE TOP-OF-ATMOSPHERE.
                          WRITE(IPR,'(/A)')' Error in RDPATH:  Top-Of'  &
     &                      //'-Atmosphere reached before end of path.'
                          STOP 'Path above Top-Of-Atmosphere in RDPATH!'
                      ENDIF
                      LAY=LAY+1
                  ENDIF
              ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).GT.ALTTOL)THEN

!                 ALTITUDE MISMATCH:
                  WRITE(IPR,'(/A,I3,A,F15.8,A,/(19X,A,I3,A,F15.8,A))')  &
     &              ' Error in RDPATH:  Path altitude',IK,' (',         &
     &              PTHALT(IK,1),' km) does not equal the next higher', &
     &              'atmospheric level (Level',LAY,' =',ZM(LAY),        &
     &              ' km) in file',REFPTH
                  STOP 'Altitude mismatch in RDPATH!'
              ELSEIF(LAY.GE.ML)THEN

!                 PATH EXTENDS ABOVE TOP-OF-ATMOSPHERE.
                  WRITE(IPR,'(/A)')' Error in RDPATH:  Top-Of'          &
     &              //'-Atmosphere reached before end of path.'
                  STOP 'Path extends above Top-Of-Atmosphere in RDPATH!'
              ELSE
                  LAY=LAY+1
              ENDIF
          ENDDO

!         VERIFY THAT FIRST PATH ALTITUDE IS IN NEXT LAYER:
          IF(PTHALT(0,1).GT.ZM(LAY)+ALTTOL)THEN
              WRITE(IPR,'(/A,F15.8,A,/(18X,A,I3,A,F15.8,A))')           &
     &          ' Error in RDPATH:  Path altitude  0 (',H1ALT(1),       &
     &          ' km) is above',' the next higher atmospheric level'    &
     &          //'(Level',LAY,' =',ZM(LAY),' km) in file',REFPTH
              STOP 'Altitude mismatch in RDPATH!'
          ENDIF

!         VERIFY THAT ALL PATH BOUNDARIES BETWEEN THE MINIMUM AND LAST
!         ALTITUDE COINCIDE WITH ADJOINING ATMOSPHERIC LEVELS:
          LAY=LAYMIN+1
          DO IK=IKHMIN(1)+1,NSEG(1)-1
              IF(2*IKHMIN(1)-IK.EQ.0)THEN

!                 VERIFY THAT THE NEXT ALTITIUDE IS H1ALT ( < H2ALT ).
                  IF(ABS(PTHALT(IK,1)-H1ALT(1)).GT.ALTTOL)THEN
                      WRITE(IPR,'(/2A,/(19X,A))')' Error in RDPATH: ',  &
     &                  ' The line-of-sight passes through a minimum'   &
     &                  //' tangent','height with H1ALT being less'     &
     &                  //' than H2ALT, but the H1ALT altitude is not', &
     &                  'a path boundary along the upward path in'      &
     &                  //' file',REFPTH
                      STOP 'Altitude mismatch in RDPATH!'
                  ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).LE.ALTTOL)THEN

!                     H1ALT COINCIDES WITH AN ATMOSPHERIC LEVEL:
                      IF(LAY.GE.ML)THEN

!                         PATH EXTENDS ABOVE TOP-OF-ATMOSPHERE.
                          WRITE(IPR,'(/A)')' Error in RDPATH:  Top-Of'  &
     &                      //'-Atmosphere reached before end of path.'
                          STOP 'Path above Top-Of-Atmosphere in RDPATH!'
                      ENDIF
                      LAY=LAY+1
                  ENDIF
              ELSEIF(ABS(PTHALT(IK,1)-ZM(LAY)).GT.ALTTOL)THEN

!                 ALTITUDE MISMATCH:
                  WRITE(IPR,'(/A,I3,A,F15.8,A,/(19X,A,I3,A,F15.8,A))')  &
     &              ' Error in RDPATH:  Path altitude',IK,' (',         &
     &              PTHALT(IK,1),' km) does not equal the next higher', &
     &              'atmospheric level (Level',LAY,' =',ZM(LAY),        &
     &              ' km) in file',REFPTH
                  STOP 'Altitude mismatch in RDPATH!'
              ELSEIF(LAY.GE.ML)THEN

!                 PATH EXTENDS ABOVE TOP-OF-ATMOSPHERE.
                  WRITE(IPR,'(/A)')' Error in RDPATH:  Top-Of'          &
     &              //'-Atmosphere reached before end of path.'
                  STOP 'Path extends above Top-Of-Atmosphere in RDPATH'
              ELSE
                  LAY=LAY+1
              ENDIF
          ENDDO

!         VERIFY THAT LAST PATH ALTITUDE IS IN NEXT LAYER:
          IF(PTHALT(NSEG(1),1).GT.ZM(LAY)+ALTTOL)THEN
              WRITE(IPR,'(/1X,A,I3,A,F15.8,A,/(18X,A,I3,A,F15.8,A))')   &
     &          'Error in RDPATH:  Path altitude',NSEG(1),' (',H2ALT(1),&
     &          ' km) is above',' the next higher atmospheric level'    &
     &          //'(Level',LAY,' =',ZM(LAY),' km) in file',REFPTH
              STOP 'Altitude mismatch in RDPATH!'
          ENDIF
      ENDIF

!     DEFINE PATH GEOMETRY PARAMETERS:
      BCKZEN(1)=180-PTHZEN(NSEG(1))
      HRANGE(1)=PTHRNG(NSEG(1),1)
      BETA(1)=PTHECA(NSEG(1))
      IF(PTHALT(NSEG(1),1).LE.PTHALT(0,1) .AND.                         &
     &  PTHALT(NSEG(1)-1,1).LT.PTHALT(NSEG(1),1))THEN
          LENN(1)=1
      ELSE
          LENN(1)=0
      ENDIF

!     RETURN TO CD3:
      RETURN
      END
