      SUBROUTINE FLXTBL

!     FLXTBL PERFORMS INITIALIZATIONS AND
!     WRITES HEADER FOR THE FLUX TABLE.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     /WTFLX/
!       UPDIFF  BOUNDARY UPWARD DIFFUSE SPECTRAL FLUX [W CM-2 / CM-1].
!       DNDIFF  BOUNDARY DOWNWARD DIFFUSE SPECTRAL FLUX [W CM-2 / CM-1].
!       DNDRCT  BOUNDARY DIRECT SOLAR SPECTRAL FLUX [W CM-2 / CM-1].
!       SMUPDF  LAYER BOUNDARY UPWARD DIFFUSE IN-BAND FLUX [W CM-2].
!       SMDNDF  LAYER BOUNDARY DOWNWARD DIFFUSE IN-BAND FLUX [W CM-2].
!       SMDNDR  LAYER BOUNDARY DIRECT SOLAR IN-BAND FLUX [W CM-2].
!       NFLUX   SPECTRAL BIN COUNTER FOR FLUX TABLE.
!       NTERMS  NUMBER OF TERMS IN FLUX SPECTRAL SUM.
      DOUBLE PRECISION UPDIFF,DNDIFF,DNDRCT,SMUPDF,SMDNDF,SMDNDR
      INTEGER NFLUX,NTERMS
      COMMON/WTFLX/UPDIFF(1:LAYDIM,-1:MWGT),DNDIFF(1:LAYDIM,-1:MWGT),   &
     &  DNDRCT(1:LAYDIM,-1:MWGT),SMUPDF(LAYDIM),SMDNDF(LAYDIM),         &
     &  SMDNDR(LAYDIM),NFLUX,NTERMS

!     /WTFLXC/
!       FRMT    FORMAT USED IN FLUX TABLE.
      CHARACTER FRMT*50
      COMMON/WTFLXC/FRMT

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

!     /TRIANG/
!       NWGT     NUMBER OF SPECTRAL BINS CONTRIBUTING TO SLIT FUNCTION.
!       IBIN0    SPECTRAL BIN OFFSET (1 FOR .2 BAND MODEL, 0 OTHERWISE).
!       WNORM    TRIANGULAR SLIT FUNCTION NORMALIZATION FACTOR.
!       WGT      NORMALIZED WEIGHTS USED TO DEFINE THE SLIT FUNCTION.
      INTEGER NWGT,IBIN0
      REAL WNORM,WGT
      COMMON/TRIANG/NWGT,IBIN0,WNORM,WGT(0:MWGT)
      SAVE /TRIANG/

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

!     /CSCAN/
!       CHUNIT   UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       RELABS   SPECTRAL RESOLUTION FLAG('A'=ABSOLUTE;'R'=RELATIVE[%]).
!       LNFEED   LINE FEED FLAG FOR .FLX FILE ('T' FOR 80 CHARACTER
!                  LINES, 'F' FOR LONG LINES, ' ' FOR NO .FLX FILE).
      CHARACTER CHUNIT*1,RELABS*1,LNFEED*1
      COMMON/CSCAN/CHUNIT,RELABS,LNFEED

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

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       IK       LAYER INDEX.
!       IBIN     SPECTRAL INDEX FOR SLIT FUNCTION.
!       UNIT     UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       WRTFLG   WRITE OUTPUT FLAG.
      INTEGER IK,IBIN
      CHARACTER UNIT*1
      LOGICAL WRTFLG

!     BRANCH BASED ON SLIT FUNCTION TYPE:
      NTERMS=0
      WRTFLG=.NOT.LJMASS .AND. NOPRNT.LE.1
      IF(.NOT.DODGRD)THEN

!         INITIALIZE SPECTRAL BIN COUNTER:
          NFLUX=0

!         INITIALIZE DOWNWARD DIRECT SOLAR FLUX ARRAY:
          DO IK=1,ML
              DO IBIN=0,NWGT
                  DNDRCT(IK,IBIN)=0.D0
              ENDDO
              SMUPDF(IK)=0.D0
              SMDNDF(IK)=0.D0
              SMDNDR(IK)=0.D0
          ENDDO

!         RETURN IF FLUX FILE IS NOT TO BE WRITTEN:
          IF(LNFEED.EQ.' ')RETURN
          IF(WRTFLG .AND. .NOT.BINOUT)THEN

!             WRITE TABLE HEADER:
              IF(XFLAG.EQ.'M')THEN
                  WRITE(IFLUX,'(///(A))')                               &
     &              ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / MICRON)',  &
     &              ' ----------------------------------------------'
                  UNIT='M'
              ELSEIF(XFLAG.EQ.'N')THEN
                  WRITE(IFLUX,'(///(A))')                               &
     &              ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / NM)',      &
     &              ' ------------------------------------------'
                  UNIT='N'
              ELSE
                  WRITE(IFLUX,'(///(A))')                               &
     &              ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / CM-1)',    &
     &              ' --------------------------------------------'
                  UNIT='W'
              ENDIF
              WRITE(IFLUX,'(/2(I3,A))')                                 &
     &          MLFLX+1,' LEVELS (',3*MLFLX+4,' TABLE COLUMNS)'
              IF(IFWHM.LE.1.5*BNDWID)THEN
                  WRITE(IFLUX,'(/A,F10.2,A,/)')' RECTANGULAR SLIT'      &
     &              //' FULL-WIDTH-AT-HALF-MAXIMUM:',IFWHM,' CM-1.'
              ELSE
                  WRITE(IFLUX,'(/A,F10.2,A,/)')' TRIANGULAR  SLIT'      &
     &              //' FULL-WIDTH-AT-HALF-MAXIMUM:',IFWHM,' CM-1.'
              ENDIF
          ENDIF

!     NON-DEFAULT SPECTRAL SCENARIO.
      ELSE

!         CURRENT FLUX DATA IS STORED IN FIRST ELEMENT:
          NFLUX=-1

!         INITIALIZE ENTIRE FLUX ARRAYS:
          DO IK=1,ML
              DO IBIN=0,MWGT
                  UPDIFF(IK,IBIN)=0.D0
                  DNDIFF(IK,IBIN)=0.D0
                  DNDRCT(IK,IBIN)=0.D0
              ENDDO
              SMUPDF(IK)=0.D0
              SMDNDF(IK)=0.D0
              SMDNDR(IK)=0.D0
          ENDDO

!         RETURN IF FLUX FILE IS NOT TO BE WRITTEN:
          IF(LNFEED.EQ.' ')RETURN

          IF(WRTFLG)THEN

!             WRITE TABLE HEADER:
              IF(CHUNIT.EQ.'M')THEN
                  IF(.NOT.BINOUT)WRITE(IFLUX,'(///(A))')                &
     &              ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / MICRON)',  &
     &              ' ----------------------------------------------'
                  UNIT='M'
                  VOUT=V2
              ELSEIF(CHUNIT.EQ.'N')THEN
                  IF(.NOT.BINOUT)WRITE(IFLUX,'(///(A))')                &
     &              ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / NM)',      &
     &              ' ------------------------------------------'
                  UNIT='N'
                  VOUT=V2
              ELSE
                  IF(.NOT.BINOUT)WRITE(IFLUX,'(///(A))')                &
     &              ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / CM-1)',    &
     &              ' --------------------------------------------'
                  UNIT='W'
                  VOUT=V1
              ENDIF
              IF(.NOT.BINOUT)WRITE(IFLUX,'(/2(I3,A))')                  &
     &          MLFLX+1,' LEVELS (',3*MLFLX+4,' TABLE COLUMNS)'
              IF(RELABS.EQ.'R')THEN
                  IF(.NOT.BINOUT)WRITE(IFLUX,'(/,A,F10.5,A,/)')         &
     &              ' TRIANGULAR  SLIT FULL-WIDTH-AT-HALF-MAXIMUM:',    &
     &              FWHM,' PERCENT.'
              ELSE
                  FWHMSQ=DBLE(FWHM)**2
                  IF(CHUNIT.EQ.'M')THEN
                      IF(.NOT.BINOUT)WRITE(IFLUX,'(/A,F10.6,A,/)')      &
     &                  ' TRIANGULAR  SLIT FULL-WIDTH-AT-HALF-MAXIMUM:',&
     &                  FWHM,' MICRONS.'
                  ELSEIF(CHUNIT.EQ.'N')THEN
                      IF(.NOT.BINOUT)WRITE(IFLUX,'(/A,F10.3,A,/)')      &
     &                  ' TRIANGULAR  SLIT FULL-WIDTH-AT-HALF-MAXIMUM:',&
     &                  FWHM,' NANOMETERS.'
                  ELSE
                      IF(.NOT.BINOUT)WRITE(IFLUX,'(/A,F10.2,A,/)')      &
     &                  ' TRIANGULAR  SLIT FULL-WIDTH-AT-HALF-MAXIMUM:',&
     &                  FWHM,' CM-1.'
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
      IF(.NOT.WRTFLG)RETURN
      IF(LNFEED.EQ.'T')THEN
          IF(.NOT.BINOUT)THEN
              WRITE(IFLUX,'(/A,F17.5,A,F33.5,A:,/(F28.5,A,F33.5,A))')   &
     &          ' ALTITUDES:',(ZM(IK),' KM',IK=1,MLFLX),ZM(ML),' KM'
              WRITE(IFLUX,'(A8,2A36)')'        ',                       &
     &          ('  ----------------------------------',IK=1,2)
              IF(UNIT.EQ.'M')THEN
                  WRITE(IFLUX,'(A8,2A36)')'  WAVLEN',                   &
     &              ('      UPWARD    DOWNWARD      DIRECT',IK=1,2)
                  WRITE(IFLUX,'(A8,2A36)')'(MICRON)',                   &
     &              ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,2)
              ELSEIF(UNIT.EQ.'N')THEN
                  WRITE(IFLUX,'(A8,2A36)')'  WAVLEN',                   &
     &              ('      UPWARD    DOWNWARD      DIRECT',IK=1,2)
                  WRITE(IFLUX,'(A8,2A36)')'   (NM) ',                   &
     &              ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,2)
              ELSE
                  WRITE(IFLUX,'(A8,2A36)')'   FREQ ',                   &
     &              ('      UPWARD    DOWNWARD      DIRECT',IK=1,2)
                  WRITE(IFLUX,'(A8,2A36)')'  (CM-1)',                   &
     &              ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,2)
              ENDIF
              WRITE(IFLUX,'(A8,2A36)')' -------',                       &
     &          (' ----------- ----------- -----------',IK=1,2)
          ENDIF
          FRMT='(F8.2,1P,  2(3(1X,E11.5)):,/(8X,  2(3(1X,E11.5))))'
      ELSE
          IF(.NOT.BINOUT)THEN
              WRITE(IFLUX,'(/A,F17.5,A,124(F33.5,A):,/(125(F28.5,A:,5X))&
     &          )')' ALTITUDES:',(ZM(IK),' KM',IK=1,MLFLX),ZM(ML),' KM'
              WRITE(IFLUX,'((8X,125A36))')                              &
     &          ('  ----------------------------------',IK=0,MLFLX)
              IF(UNIT.EQ.'M')THEN
                  WRITE(IFLUX,'(A8,125A:,/(8X,125A))')'  WAVLEN',       &
     &              ('      UPWARD    DOWNWARD      DIRECT',IK=0,MLFLX)
                  WRITE(IFLUX,'(A8,125A:,/(8X,125A))')'(MICRON)',       &
     &              ('     DIFFUSE     DIFFUSE       SOLAR',IK=0,MLFLX)
              ELSEIF(UNIT.EQ.'N')THEN
                  WRITE(IFLUX,'(A8,125A:,/(8X,125A))')'  WAVLEN',       &
     &              ('      UPWARD    DOWNWARD      DIRECT',IK=0,MLFLX)
                  WRITE(IFLUX,'(A8,125A:,/(8X,125A))')'   (NM) ',       &
     &              ('     DIFFUSE     DIFFUSE       SOLAR',IK=0,MLFLX)
              ELSE
                  WRITE(IFLUX,'(A8,125A:,/(8X,125A))')'   FREQ ',       &
     &              ('      UPWARD    DOWNWARD      DIRECT',IK=0,MLFLX)
                  WRITE(IFLUX,'(A8,125A:,/(8X,125A))')'  (CM-1)',       &
     &              ('     DIFFUSE     DIFFUSE       SOLAR',IK=0,MLFLX)
              ENDIF
              WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')' -------',       &
     &          ('  ----------  ----------  ----------',IK=0,MLFLX)
          ENDIF
          FRMT='(F8.2,1P,125(3(1X,E11.5)):,/(8X,125(3(1X,E11.5))))'
      ENDIF

!     WRITE FLUX TABLE DATA FORMAT:
      IF(UNIT.EQ.'M')THEN
          IF(IV1.GT.100.)THEN
              FRMT(5:5)='5'
          ELSEIF(IV1.GT.10.)THEN
              FRMT(5:5)='4'
          ELSE
              FRMT(5:5)='3'
          ENDIF
      ELSEIF(UNIT.EQ.'N')THEN
          IF(IV1.GT.100.)THEN
              FRMT(5:5)='2'
          ELSEIF(IV1.GT.10.)THEN
              FRMT(5:5)='1'
          ELSE
              FRMT(5:5)='0'
          ENDIF
      ENDIF
      IF(BINOUT)CALL BNFLX3(MLFLX,ML,NFLUX,FWHM,IV1,ZM,                 &
     &  IFWHM,BNDWID,XFLAG,CHUNIT,RELABS,LNFEED,FRMT)

!     RETURN TO TRANS.
      RETURN
      END
