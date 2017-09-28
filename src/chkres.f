      SUBROUTINE CHKRES(RESCHR,FWSCAN,DEGALL)

!     CHKRES DETERMINES THE COMPUTATIONAL SPECTRAL RANGE
!     REQUIRED TO OUTPUT THE CONVOLVED SPECTRAL DATA.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       RESCHR   BAND MODEL RESOLUTION CHARACTER STRING.
!       FWSCAN   NUMBER OF FULL-WIDTH AT HALF-MAXIMUM (FWHM)
!                USED TO INTEGRATE SCANNING FUNCTION
!       DEGALL   LOGICAL FLAG, TRUE IF ALL SPECTRAL OUTPUT IS CONVOLVED.
      CHARACTER RESCHR*2
      REAL FWSCAN
      LOGICAL DEGALL

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

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

!     /CARD4/
!       IV1      LOWEST SPECTRAL FREQUENCY OUTPUT [CM-1].
!       IV2      HIGHEST SPECTRAL FREQUENCY OUTPUT [CM-1].
!       IDV      PRINTOUT SPECTRAL FREQUENCY STEP SIZE [CM-1].
!       IFWHM    TRIANGULAR SLIT FULL-WIDTH-HALF-MAXIMUM [CM-1].
!       VBAND    CURRENT COMPUTATION BAND FREQUENCY [CM-1].
!                (EQUALS BAND CENTER FOR 1, 5 & 15 CM-1 BAND MODELS;
!                EQUALS THE MINIMUM BAND VALUE FOR 0.1 CM-1 BAND MODEL)
!       IBINPT   BIN NUMBER OF CURRENT SPECTRAL POINT.
!                (CENTER FREQUENCY = IBINPT * BNDWID + OSHIFT).
!       IBINLO   BIN NUMBER OF (PADDED) SPECTRAL RANGE LOWER BOUND.
!       IBINHI   BIN NUMBER OF (PADDED) SPECTRAL RANGE UPPER BOUND.
!       IBINMN   BIN NUMBER OF MINIMUM COMPUTATION SPECTRAL POINT.
!       IBINMX   BIN NUMBER OF MAXIMUM COMPUTATION SPECTRAL POINT.
!       IBINDL   BIN NUMBER INCREMENT FOR SPECTRAL PRINTOUT.
!       IBINRS   BIN NUMBER INCREMENT EQUAL TO SPECTRAL RESOLUTION.
!       IBINOS   BIN NUMBER OFFSET BETWEEN CURRENT & OUTPUT SPC POINTS.
!       IBINWR   BIN NUMBER OF NEXT SPECTRAL DATA WRITE.
!       MBINPT   BIN NUMBER MAXIMUM FOR CURRENT BAND MODEL RESOLUTION.
!       IDBIN5   SPECTRAL BIN NUMBER STEP SIZE FOR 5 CM-1 GRID.
!       ISTEP5   INCREMENT FOR RETRIEVING 5 CM-1 RESOLUTION DATA [CM-1].
!       NSPCDT   NUMBER OF OUTPUT SPECTRAL DATA POINTS.
      DOUBLE PRECISION IDV
      REAL IV1,IV2,IFWHM,VBAND
      INTEGER IBINPT,IBINLO,IBINHI,IBINMN,IBINMX,IBINDL,                &
     &  IBINRS,IBINOS,IBINWR,MBINPT,IDBIN5,ISTEP5,NSPCDT
      COMMON/CARD4/IDV,IV1,IV2,IFWHM,VBAND,IBINPT,IBINLO,IBINHI,IBINMN, &
     &  IBINMX,IBINDL,IBINRS,IBINOS,IBINWR,MBINPT,IDBIN5,ISTEP5,NSPCDT

!     /SCAN/
!       V1       LOWER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       V2       UPPER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       DV       SPECTRAL STEP SIZE FOR OUTPUT [CHUNIT DEFINES UNIT].
!       FWHM     FULL-WIDTH-AT-HALF-MAXIMUM [CHUNIT DEFINES UNIT].
!       FWHMSQ   TRIANGULAR SLIT NORMALIZATION FACTOR
!                (EQUALS FWHM SQUARED) [CHUNIT DEFINES UNIT].
!       VOUT     CURRENT SPECTRAL OUTPUT [CHUNIT DEFINES UNIT].
!       M2_FAC   FACTOR IN 2ND MOMENT DIFFERENCE INTEGRAL [CHUNIT UNIT].
      DOUBLE PRECISION V1,V2,DV,VOUT,FWHMSQ
      REAL FWHM,M2_FAC
      COMMON/SCANFN/V1,V2,DV,VOUT,FWHMSQ,FWHM,M2_FAC

!     /CFLAGS/
!       YFLAG    Y COORDINATE FLAG FOR PLOT.DAT FILE
!                  = "T" FOR TRANSMITTANCE
!                  = "R" FOR RADIANCE (IRRADIANCE FOR IEMSCT=3)
!                  = "N" FOR NO PLOT.DAT OUTPUT
!       XFLAG    X COORDINATE FLAG FOR PLOT.DAT FILE
!                  = "W" FOR FREQUENCY IN WAVENUMBERS (CM-1) AND
!                        RADIANCE IN W SR-1 CM-2 / CM-1
!                  = "M" FOR WAVELENGTH IN MICRONS AND
!                        RADIANCE IN W SR-1 CM-2 / MICRON
!                  = "N" FOR WAVELENGTH IN NANOMETERS AND
!                        RADIANCE IN MICRO-WATTS SR-1 CM-2 / NANOMETER
!       DLIMIT   DELIMITER CHARACTER STRING BETWEEN MODTRAN RUNS
!       FLAGS    SCANNING FUNCTION FLAGS.
      CHARACTER YFLAG*1,XFLAG*1,DLIMIT*8,FLAGS*7
      COMMON/CFLAGS/YFLAG,XFLAG,DLIMIT,FLAGS

!     /CSCAN/
!       CHUNIT   UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       RELABS   SPECTRAL RESOLUTION FLAG('A'=ABSOLUTE;'R'=RELATIVE[%]).
!       LNFEED   LINE FEED FLAG FOR .FLX FILE ('T' FOR 80 CHARACTER
!                  LINES, 'F' FOR LONG LINES, ' ' FOR NO .FLX FILE).
      CHARACTER CHUNIT*1,RELABS*1,LNFEED*1
      COMMON/CSCAN/CHUNIT,RELABS,LNFEED

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       SPCPAD   FILTER SPECTRAL PADDING [CHUNIT DEFINES UNIT].
!       WNMIN    MINIMUM SPECTRAL FREQUENCY [CM-1].
!       VSTORE   TEMPORARY STORAGE OF SPECTRAL INPUT.
!       VMAX     MAXIMUM MODTRAN SPECTRAL FREQUENCY.
!       WPERCM   WAVELENGTHS PER CENTIMETER.
      DOUBLE PRECISION SPCPAD,WNMIN,VSTORE,VMAX,WPERCM

!     USER INPUT INSTRUMENT BANDPASS FREQUENCY (INPUTS):
!     V1 IS NUMERICALLY MINIMUM IN CM-1, MICRONS OR NM
!     V2 IS NUMERICALLY MAXIMUM IN CM-1, MICRONS OR NM
!     V1 AND V2 MAY BE REDEFINED IF REQUIRED

!     UNITS OF SCANNING FUNCTION, V1 AND V2 (INPUTS):
!     CHUNIT='W' (CM-1)
!     CHUNIT='M' (MICRON)
!     CHUNIT='N' (NM)

!     RELATIVE OR ABSOLUTE FWHM (INPUTS):
!     RELABS='A' OR BLANK; FWHM IS IN CM-1, MICRON OR NM
!     RELABS='R'; FWHM IS IN PERCENT

!     MAXIMUM FREQUENCY AND BIN:
      IF(RESCHR.EQ.'15')THEN
          VMAX=49995.D0
      ELSE
          VMAX=50000.D0
      ENDIF
      MBINPT=NINT(VMAX/DBLE(BNDWID))

!     IF FLAGS(4:4)='A', DEGRADE EVERYTHING IN TAPE7:
      DEGALL=FLAGS(4:4).EQ.'A'

!     ANALYZE VALIDITY AND MEANING OF NUMBERS IN CARD4.
!     CONVERT TO EQUIVALENT INTEGERS IN CM-1 PRIOR TO THE MODS.
!     THE PURPOSE OF THESE MODS ARE ENUMERATED BELOW:

!     1.  RELABS = 'A' OR " ", FULL-WIDTH-AT-HALF-MAX IS ABSOLUTE
!         RELABS='R', FWHM IS IN PERCENT

!     2.  CHUNIT DENOTES THE UNITS OF FREQUENCY AND RADIANCE
!         CHUNIT= 'W', V1, V2, ETC IN CM-1, RADIANCE IN W/CM^2/SR/CM-1
!         CHUNIT= 'M', V1, ETC IN MICRON, RADIANCE IN W/CM^2/SR/MICRON
!         CHUNIT= 'N', V1, ETC IN NM, RADIANCE IN MICROWATT/CM^2/SR/NM

!     CHECK RELATIVE/ABSOLUTE FLAG:
      RELABS=FLAGS(3:3)
      IF(RELABS.EQ.' ')THEN

!         SET RELABS TO DEFAULT VALUE:
          RELABS='A'
      ELSEIF(RELABS.NE.'A' .AND. RELABS.NE.'R')THEN

!         RELABS MUST BE 'R' OR 'A'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'BAD INPUT FOR RELABS IN CARD 4'
      ENDIF

!     CHECK UNIT FLAG:
      CHUNIT=FLAGS(1:1)
      IF(CHUNIT.EQ.' ')CHUNIT='W'
      IF(CHUNIT.NE.'W' .AND. CHUNIT.NE.'M' .AND. CHUNIT.NE.'N')THEN
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'BAD INPUT FOR CHUNIT IN CARD 4'
      ENDIF

!     CARRIAGE RETURN FLAG (FOR .FLX FILE).
      IF(FLAGS(7:7).EQ.'T' .OR. FLAGS(7:7).EQ.'F')THEN
          LNFEED=FLAGS(7:7)
      ELSE
          LNFEED=' '
      ENDIF

!     SWAP V2 AND V1 IF OUT OF ORDER:
      IF(V2.LT.V1)THEN
          VSTORE=V1
          V1=V2
          V2=VSTORE
      ENDIF

!     BRANCH BASED ON SPECTRAL FREQUENCY(CM-1) OR WAVELENGTH(uM or nm):
      IF(CHUNIT.EQ.'W')THEN

!         SPECTRAL RANGE IN CM-1 UNITS:
          IF(RELABS.EQ.'A')THEN

!             IF DEFAULT SPECTRAL SCENARIO, PADDING IS NOT NECESSARY:
              IF(.NOT.DODGRD)THEN
                  IBINLO=NINT(V1/DBLE(BNDWID))
                  IBINHI=NINT(V2/DBLE(BNDWID))
                  IBINDL=NINT(DV/DBLE(BNDWID))
                  IBINRS=NINT(FWHM/BNDWID)

!                 CHECK SPECTRAL INPUTS AND RETURN:
                  IF(IBINDL.LT.1)THEN
                      WRITE(IPR,'(3(A,F12.6))')                         &
     &                  ' Warning from CHKRES:  DV is being'//          &
     &                  ' increased from',DV,' to',BNDWID,' CM-1.'
                      DV=DBLE(BNDWID)
                      IBINDL=1
                  ELSEIF(IBINDL.GT.50)THEN
                      WRITE(IPR,'(3(A,F12.6))')                         &
     &                  ' Warning from CHKRES:  DV is being'//          &
     &                  ' decreased from',DV,' to',50*BNDWID,' CM-1.'
                      DV=50*DBLE(BNDWID)
                      IBINDL=50
                  ENDIF
                  IF(IBINRS.LT.1)THEN
                      WRITE(IPR,'(3(A,F12.6))')                         &
     &                  ' Warning from CHKRES:  FWHM is being'//        &
     &                  ' increased from',FWHM,' to',BNDWID,' CM-1.'
                      FWHM=BNDWID
                      IBINRS=1
                  ELSEIF(IBINRS.GT.50)THEN
                      WRITE(IPR,'(3(A,F12.6))')                         &
     &                  ' Warning from CHKRES:  FWHM is being'//        &
     &                  ' decreased from',FWHM,' to',50*BNDWID,' CM-1.'
                      FWHM=50*BNDWID
                      IBINRS=50
                  ENDIF
                  IF(IBINDL.GT.IBINRS)WRITE(IPR,'(3(A,F12.6))')         &
     &              ' Warning from CHKRES:  Output spectral'            &
     &              //' step size (',DV,' CM-1) exceeds the'            &
     &              //' spectral resolution (',FWHM,'CM-1).'
                  IF(IBINLO.LT.0)THEN
                      WRITE(IPR,'(3(A,F12.6))')                         &
     &                  ' Warning from CHKRES:  V1 is being'//          &
     &                  ' increased from',V1,' to 0. CM-1.'
                      V1=0.D0
                      IBINLO=0
                  ENDIF
                  IF(IBINHI.GT.MBINPT)THEN
                      WRITE(IPR,'(3(A,F12.6))')                         &
     &                  ' Warning from CHKRES:  V2 is being'//          &
     &                  ' decreased from',V2,' to',VMAX,' CM-1.'
                      V2=VMAX
                      IBINHI=MBINPT
                  ENDIF
                  IF(IBINHI.LT.IBINLO+IBINDL)THEN
                      IF(IBINLO+IBINDL.LE.MBINPT)THEN
                          IBINHI=IBINLO+IBINDL
                          WRITE(IPR,'(3(A,F12.6))')' Warning from'//    &
     &                      ' CHKRES:  V2 is being increased from',     &
     &                      V2,' to',IBINHI*BNDWID,' CM-1.'
                          V2=DBLE(IBINHI*BNDWID)
                      ELSE
                          IBINLO=IBINHI-IBINDL
                          WRITE(IPR,'(3(A,F12.6))')' Warning from'//    &
     &                      ' CHKRES:  V1 is being decreased from',     &
     &                      V1,' to',IBINLO*BNDWID,' CM-1.'
                          V1=DBLE(IBINLO*BNDWID)
                      ENDIF
                  ENDIF
                  RETURN
              ENDIF

!             THE FWHM IS AN ABSOLUTE (NOT RELATIVE) RESOLUTION:
              IF(FWHM.LT.BNDWID)THEN
                  WRITE(IPR,'(/A,/18X,A,F12.5,A)')                      &
     &              ' Error in CHKRES:  Specified spectral'//           &
     &              ' resolution (FWHM on CARD4) is too fine.',         &
     &              ' It must be greater or equal to',BNDWID,' CM-1.'
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP 'Error in CHKRES:  Spectral resolution too fine.'
              ENDIF

!             CHECK LOWER SPECTRAL FREQUENCY LIMIT:
              IF(LNFEED.NE.' ' .AND. FWSCAN.LE.1.)THEN

!                 WITH A RECTANGULAR SCANNING FUNCTION, THE TRIANGULAR
!                 SLIT OF THE *.FLX FILE DICTATES THE REQUIRED PADDING:
                  SPCPAD=DBLE(FWHM)
              ELSE
                  SPCPAD=DBLE(FWSCAN*FWHM)
              ENDIF
              IF(V1.LT.SPCPAD)THEN
                  WRITE(IPR,'(3(A,F12.6))')                             &
     &              ' Warning from CHKRES:  V1 is being'//              &
     &              ' increased from',V1,' to',SPCPAD,' CM-1.'
                  V1=SPCPAD
              ENDIF
              IBINLO=NINT((V1-SPCPAD)/DBLE(BNDWID))-1

!             CHECK UPPER SPECTRAL FREQUENCY LIMIT:
              IF(V2.GT.VMAX-SPCPAD)THEN
                  WRITE(IPR,'(3(A,F12.6))')                             &
     &              ' Warning from CHKRES:  V2 is being'//              &
     &              ' reduced from',V2,' to',VMAX-SPCPAD,' CM-1.'
                  V2=VMAX-SPCPAD
                  IF(V2.LE.V1)THEN
                      WRITE(IPR,'(/A)')' Error in CHKRES:  New V2 < V1.'
                      STOP ' Error in CHKRES:  New V2 < V1.'
                  ENDIF
              ENDIF
              IBINHI=NINT(REAL(V2+SPCPAD)/BNDWID)+1
          ELSE

!             THE FWHM(%) IS A RELATIVE (NOT ABSOLUTE) RESOLUTION:
              IF(REAL(V1)*FWHM/100.LT.BNDWID .OR. FWSCAN*FWHM.GE.100.)  &
     &          THEN
                  WRITE(IPR,'(/A,/17X,A,2(F12.5,A))')'Error in CHKRES: '&
     &              //' Specified spectral resolution (FWHM on CARD4)', &
     &              ' is out of range.  It must exceed',                &
     &              100*BNDWID/REAL(V1),'%, but be less than',          &
     &              100/FWSCAN,'%.'
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP 'Error in CHKRES:  Spectral resolution too fine.'
              ENDIF
              IF(LNFEED.NE.' ' .AND. FWSCAN.LE.1.)THEN

!                 WITH A RECTANGULAR SCANNING FUNCTION, THE TRIANGULAR
!                 SLIT OF THE *.FLX FILE DICTATES THE REQUIRED PADDING:
                  SPCPAD=V1*DBLE(FWHM)/100
                  IBINLO=NINT((V1-SPCPAD)/DBLE(BNDWID))-1
                  SPCPAD=V2*DBLE(FWHM)/100
              ELSE
                  SPCPAD=V1*DBLE(FWSCAN*FWHM)/100
                  IBINLO=NINT((V1-SPCPAD)/DBLE(BNDWID))-1
                  SPCPAD=V2*DBLE(FWSCAN*FWHM)/100
              ENDIF

!             CHECK UPPER SPECTRAL FREQUENCY LIMIT:
              IF(V2.GT.VMAX-SPCPAD)THEN
                  WRITE(IPR,'(3(A,F12.6))')                             &
     &              ' Warning from CHKRES:  V2 is being reduced from',  &
     &              V2,' to',VMAX/DBLE(1-FWSCAN*FWHM/100),' CM-1.'
                  V2=VMAX/DBLE(1-FWSCAN*FWHM/100)
                  IF(V2.LE.V1)THEN
                      WRITE(IPR,'(/A)')' Error in CHKRES:  New V2 < V1.'
                      STOP ' Error in CHKRES:  New V2 < V1.'
                  ENDIF
              ENDIF
              IBINHI=NINT(REAL(V2+SPCPAD)/BNDWID)+1
          ENDIF
          WNMIN=V1

      ELSE

!         WAVELENGTH (MICRONS OR NM) SPECTRAL RANGE
          IF(CHUNIT.EQ.'M')THEN

!             MICRON INPUTS:
              WPERCM=1.D4
          ELSE

!             NANOMETER INPUTS:
              WPERCM=1.D7
          ENDIF
          WNMIN=WPERCM/V2

!         RELATIVE OR ABSOLUTE RESOLUTION INPUT?
          IF(RELABS.EQ.'A')THEN

!             FWHM IS AN ABSOLUTE RESOLUTION.  CHECK IF SUFFICIENT:
              IF(DBLE(FWHM)*WNMIN.LT.V2*DBLE(BNDWID))THEN
                  IF(CHUNIT.EQ.'M')THEN
                      WRITE(IPR,'(/2A,/18X,A,F12.5,A)')                 &
     &                  ' Error in CHKRES:  Specified spectral'//       &
     &                  ' resolution (FWHM on CARD4) is too fine for',  &
     &                  ' the chosen band',' model; it must exceed',    &
     &                  V2*DBLE(BNDWID)/WNMIN,' Microns.'
                  ELSE
                      WRITE(IPR,'(/2A,/18X,A,F12.2,A)')                 &
     &                  ' Error in CHKRES:  Specified spectral'//       &
     &                  ' resolution (FWHM on CARD4) is too fine for',  &
     &                  ' the chosen band',' model; it must exceed',    &
     &                  V2*DBLE(BNDWID)/WNMIN,' Nanometers.'
                  ENDIF
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP 'Error in CHKRES:  Spectral resolution too fine.'
              ENDIF
              IF(LNFEED.NE.' ' .AND. FWSCAN.LE.1.)THEN

!                 WITH A RECTANGULAR SCANNING FUNCTION, THE TRIANGULAR
!                 SLIT OF THE *.FLX FILE DICTATES THE REQUIRED PADDING:
                  SPCPAD=DBLE(FWHM)
              ELSE
                  SPCPAD=DBLE(FWSCAN*FWHM)
              ENDIF
              IBINLO=NINT(WPERCM/((V2+SPCPAD)*DBLE(BNDWID)))-1

!             CHECK UPPER SPECTRAL FREQUENCY (LOWER WAVELENGTH) LIMIT:
              VSTORE=WPERCM/VMAX+SPCPAD
              IF(VSTORE.GT.V1)THEN
                  IF(CHUNIT.EQ.'M')THEN
                      WRITE(IPR,'(3(A,F12.6))')                         &
     &                  ' Warning from CHKRES:  V1 is being increased'  &
     &                  //' from',V1,' to',VSTORE,' Microns.'
                  ELSE
                      WRITE(IPR,'(3(A,F12.3))')                         &
     &                  ' Warning from CHKRES:  V1 is being increased'  &
     &                  //' from',V1,' to',VSTORE,' Nanometers.'
                  ENDIF
                  V1=VSTORE
                  IF(V1.GE.V2)THEN
                      WRITE(IPR,'(/A)')' Error in CHKRES:  New V1 > V2.'
                      STOP ' Error in CHKRES:  New V1 > V2.'
                  ENDIF
              ENDIF
              IBINHI=NINT(WPERCM/((V1-SPCPAD)*DBLE(BNDWID)))+1
          ELSE

!             FWHM IS IN PERCENT MICRON OR NM.  CHECK RESOLUTION LIMIT:
              IF(REAL(WNMIN)*FWHM/100.LT.BNDWID                         &
     &          .OR. FWSCAN*FWHM.GE.100.)THEN
                  WRITE(IPR,'(/A,/17X,A,2(F12.5,A))')                   &
     &              ' Error in CHKRES:  Specified spectral'//           &
     &              ' resolution (FWHM on CARD4) is out of range.',     &
     &              ' It must exceed',100*DBLE(BNDWID)/WNMIN,           &
     &              '%, but be less than',100/FWSCAN,'%.'
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP                                                  &
     &              'Error in CHKRES: Spectral resolution out of range.'
              ENDIF
              IF(LNFEED.NE.' ' .AND. FWSCAN.LE.1.)THEN

!                 WITH A RECTANGULAR SCANNING FUNCTION, THE TRIANGULAR
!                 SLIT OF THE *.FLX FILE DICTATES THE REQUIRED PADDING:
                  SPCPAD=V2*DBLE(FWHM)/100
                  IBINLO=NINT(WPERCM/((V2+SPCPAD)*DBLE(BNDWID)))-1
                  SPCPAD=V1*DBLE(FWHM)/100
              ELSE
                  SPCPAD=V2*DBLE(FWSCAN*FWHM)/100
                  IBINLO=NINT(WPERCM/((V2+SPCPAD)*DBLE(BNDWID)))-1
                  SPCPAD=V1*DBLE(FWSCAN*FWHM)/100
              ENDIF

!             CHECK UPPER SPECTRAL FREQUENCY (LOWER WAVELENGTH) LIMIT:
              VSTORE=WPERCM/VMAX+SPCPAD
              IF(VSTORE.GT.V1)THEN
                  WRITE(IPR,'(3(A,F12.6))')                             &
     &              ' Warning from CHKRES:  V1 is being increased from',&
     &              V1,' to',VSTORE,' MICRONS.'
                  V1=VSTORE
                  IF(V1.GE.V2)THEN
                      WRITE(IPR,'(/A)')' Error in CHKRES:  New V1 > V2.'
                      STOP ' Error in CHKRES:  New V1 > V2.'
                  ENDIF
              ENDIF

!             UPPER LIMIT SHOULD BE A MULTIPLE OF 5 CM-1.
              IBINHI=NINT(WPERCM/((V1-SPCPAD)*DBLE(BNDWID)))+1

          ENDIF
      ENDIF

!     CHECK FOR LOWTRAN BAND MODEL:
      IF(WNMIN.GT.DBLE(MXFREQ) .OR. .NOT.MODTRN)WRITE(IPR,'(/A)')       &
     &  'Warning from CHKRES:  LOWTRAN (not MODTRAN) band model.'

!     MAKE MINIMUM & MAXIMUM SPECTRAL FREQUENCIES MULTIPLES OF 5 CM-1.
      IF(RESCHR.EQ.'01')THEN
          IBINLO=5*(IBINLO/5)
          IBINHI=5*((IBINHI+4)/5)
      ELSEIF(RESCHR.EQ.'p1')THEN
          IBINLO=50*(IBINLO/50)
          IBINHI=50*((IBINHI+49)/50)
      ENDIF
      IF(IBINLO.LE.0)IBINLO=0
      IF(IBINHI.GT.MBINPT)IBINHI=MBINPT
      IF(IBINLO.GE.IBINHI)THEN
          WRITE(IPR,'(/A)')'Error in CHKRES: Need more spectral points.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error in CHKRES: Need more spectral points.'
      ENDIF
      IBINDL=1
      IBINRS=1

      RETURN
      END

      SUBROUTINE REBIN(LSKIP,IEMSCT,IPUSCN)

!     TO REUSE TAPE7.SCR, READ OLD HEADER FROM TAPE7.SCN (OR *.7SC).
!     IV1 AND IV2 INPUTS TO THIS ROUTINE ARE FOR CURRENT CARD 4.
!     BUT THESE MAY NOT MATCH THE PARAMETERS OF THE EARLIER RUN.
!     SO RESET IV1 AND IV2 TO PARAMETERS OF EARLIER RUN.
!     TAPE7B.SCR (<ROOTNAME>B.7SR) CONTAINS VALUES FROM IV1 TO IV2 CM-1 IN
!     INCREMENTS OF 1
!     UNLESS IV1 AND IV2 ARE KNOWN THIS FILE CANNOT BE READ PROPERLY.

!     IDV INPUT MUST MATCH WHAT IS IN TAPE7B.SCR.
!     IF NOT, RUN MODTRAN AS THE OLD TAPE7B.SCR IS NO GOOD.

!     All format/text messages in the file are preceeded by two records,
!     each with -9999.  All text messages are 80 byte long.
      IMPLICIT NONE

!     ARGUMENTS:
      INTEGER IEMSCT,IPUSCN
      LOGICAL LSKIP

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'BMHEAD.h'

!     /CARD4/
!       IV1      LOWEST SPECTRAL FREQUENCY OUTPUT [CM-1].
!       IV2      HIGHEST SPECTRAL FREQUENCY OUTPUT [CM-1].
!       IDV      PRINTOUT SPECTRAL FREQUENCY STEP SIZE [CM-1].
!       IFWHM    TRIANGULAR SLIT FULL-WIDTH-HALF-MAXIMUM [CM-1].
!       VBAND    CURRENT COMPUTATION BAND FREQUENCY [CM-1].
!                (EQUALS BAND CENTER FOR 1, 5 & 15 CM-1 BAND MODELS;
!                EQUALS THE MINIMUM BAND VALUE FOR 0.1 CM-1 BAND MODEL)
!       IBINPT   BIN NUMBER OF CURRENT SPECTRAL POINT.
!                (CENTER FREQUENCY = IBINPT * BNDWID + OSHIFT).
!       IBINLO   BIN NUMBER OF (PADDED) SPECTRAL RANGE LOWER BOUND.
!       IBINHI   BIN NUMBER OF (PADDED) SPECTRAL RANGE UPPER BOUND.
!       IBINMN   BIN NUMBER OF MINIMUM COMPUTATION SPECTRAL POINT.
!       IBINMX   BIN NUMBER OF MAXIMUM COMPUTATION SPECTRAL POINT.
!       IBINDL   BIN NUMBER INCREMENT FOR SPECTRAL PRINTOUT.
!       IBINRS   BIN NUMBER INCREMENT EQUAL TO SPECTRAL RESOLUTION.
!       IBINOS   BIN NUMBER OFFSET BETWEEN CURRENT & OUTPUT SPC POINTS.
!       IBINWR   BIN NUMBER OF NEXT SPECTRAL DATA WRITE.
!       MBINPT   BIN NUMBER MAXIMUM FOR CURRENT BAND MODEL RESOLUTION.
!       IDBIN5   SPECTRAL BIN NUMBER STEP SIZE FOR 5 CM-1 GRID.
!       ISTEP5   INCREMENT FOR RETRIEVING 5 CM-1 RESOLUTION DATA [CM-1].
!       NSPCDT   NUMBER OF OUTPUT SPECTRAL DATA POINTS.
      DOUBLE PRECISION IDV
      REAL IV1,IV2,IFWHM,VBAND
      INTEGER IBINPT,IBINLO,IBINHI,IBINMN,IBINMX,IBINDL,                &
     &  IBINRS,IBINOS,IBINWR,MBINPT,IDBIN5,ISTEP5,NSPCDT
      COMMON/CARD4/IDV,IV1,IV2,IFWHM,VBAND,IBINPT,IBINLO,IBINHI,IBINMN, &
     &  IBINMX,IBINDL,IBINRS,IBINOS,IBINWR,MBINPT,IDBIN5,ISTEP5,NSPCDT

!     LOCAL VARIABLES:
      CHARACTER CTEMP*80
      INTEGER K,IEM1,INMX,INMY,IIMULT,IDIS,ICOUNT
      REAL FREQ1,FREQ2

!     FIND WHERE TAPE7.SCN OR <ROOTNAME>.7SC HEADER ENDS (UNIT=IPUSCN):
   10 CONTINUE
      READ(IPUSCN,END=40)CTEMP
      IF (CTEMP(1:7).NE.TXTSPR)GOTO 10
      READ(IPUSCN,END=40)CTEMP
      IF (CTEMP(1:7).NE.TXTSPR)GOTO 10
      READ(IPUSCN,ERR=40,END=40)IEM1,INMX,INMY,IIMULT,IDIS
      DO K=1,INMX
          READ(IPUSCN)
      ENDDO
      DO K=1,INMY
          READ(IPUSCN)
      ENDDO

      IF(IEMSCT.EQ.0)READ(IPUSCN)

      READ(IPUSCN)FREQ1
      ICOUNT=1
   20 CONTINUE
      READ(IPUSCN,ERR=30,END=30)FREQ2
      ICOUNT=ICOUNT+1
      GOTO 20
   30 CONTINUE
      REWIND(IPUSCN)
      IF(ABS(IDV-DBLE(FREQ2-FREQ1)/ICOUNT).LT..033D0)THEN
          LSKIP=.TRUE.
          IBINLO=NINT(FREQ1/BNDWID-.1)
          IBINHI=NINT(FREQ2/BNDWID-.1)
          IV1=IBINLO*BNDWID
          IV2=IBINHI*BNDWID
          RETURN
      ENDIF

!     BYPASS WHEN "FREQ" CANNOT BE FOUND:
   40 CONTINUE
      LSKIP=.FALSE.
      RETURN
      END
