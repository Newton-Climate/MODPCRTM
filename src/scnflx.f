      SUBROUTINE SCNFLX(VCEN)

!     SCNFLX WRITES OUT SPECTRAL FLUX, DEGRADED WITH
!     A SCANNING FUNCTION, TO UNIT IFLUX.

!     INPUT ARGUMENT:
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
      REAL VCEN

!     PARAMETERS:
      DOUBLE PRECISION SMALL
      PARAMETER(SMALL=1.D-7)
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'BMHEAD.h'

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

!     FUNCTIONS:
!       WTRGT    RIGHT SIDE WAVELENGTH TRIANGULAR SLIT INTEGRATION.
!       WTLFT    LEFT SIDE WAVELENGTH TRIANGULAR SLIT INTEGRATION.
      DOUBLE PRECISION WTRGT,WTLFT

!     LOCAL VARIABLES:
!       JFLUX    SPECTRAL INDEX FOR OUTPUT FLUX VALUES.
!       JOFF     SPECTRAL INDEX OFFSET (EQUALS NUMBER OF SLIT
!                FUNCTIONS TERMINATING IN CURRENT SPECTRAL BIN).
!       LAST     LOGICAL FLAG, TRUE IF CURRENT SLIT FUNCTION
!                TERMINATES IN CURRENT SPECTRAL BIN.
!       CONVRT   SPECTRAL CONVERSION FACTOR [WAVELENGTHS / CM].
!       VOUTCN   SLIT FUNCTION CENTRAL WAVELENGTH [MICRONS OR NM].
!       VOUTMN   SLIT FUNCTION MINIMUM WAVELENGTH [MICRONS OR NM].
!       VOUTMX   SLIT FUNCTION MAXIMUM WAVELENGTH [MICRONS OR NM].
!       WNCEN    SLIT FUNCTION CENTRAL FREQUENCY [CM-1].
!       WNMIN    SLIT FUNCTION MINIMUM FREQUENCY [CM-1].
!       WNMAX    SLIT FUNCTION MAXIMUM FREQUENCY [CM-1].
!       WNFWHM   FULL-WIDTH-AT-HALF-MAXIMUM [CM-1].
!       BINMIN   SPECTRAL BIN MINIMUM FREQUENCY [CM-1].
!       BINMAX   SPECTRAL BIN MAXIMUM FREQUENCY [CM-1].
!       CONFLX   FLUX SPECTRAL CONVERSION FACTOR [CM-1 / SPECTRAL UNIT].
!       DIFF     FREQUENCY DIFFERENCE OF INTEGRAL LIMITS [CM-1].
!       COEF     COEFFICIENT USED TO EVALUATE INTEGRAL [WAVELENGTH CM].
!       AVG      AVERAGE OF INTEGRAL LIMITS [CM-1].
!       WT       WEIGHT OF SPECTRAL BIN IN EVALUATING SLIT FUNCTION.
      INTEGER JFLUX,JOFF
      LOGICAL LAST
      DOUBLE PRECISION VOUTCN,WNCEN,CONVRT,BINMIN,BINMAX,WNMIN,WNMAX,   &
     &  DIFF,VOUTMN,VOUTMX,WNFWHM,CONFLX,COEF,AVG,WT

!     THERE ARE 10 POSSIBLE OVERLAP SCENARIOS BETWEEN THE TRIANGULAR
!     SLIT AND THE SPECTRAL BIN (THE ABSCISSA IS WAVELENGTH):

!                                  / \
!                                 / | \
!                                /  |  \
!                               /   |   \
!                              /    |    \
!                             /     |     \
!                            /      |      \
!                           /       |       \
!                          /        |        \
!                         /         |         \

!     CASE A:                                      |__|

!     CASE B:                            |____________|

!     CASE C:               |_________________________|

!     CASE D:       |_________________________________|

!     CASE E:                            |__|

!     CASE F:               |_______________|

!     CASE C:       |_______________________|

!     CASE E:               |__|

!     CASE I:       |__________|

!     CASE J:       |__|

!     SKIP WRITING FLUX FILE IF LNFEED IS BLANK:
      IF(LNFEED.EQ.' ')RETURN

!     BRANCH BASED ON UNITS:
      IF(CHUNIT.EQ.'W')THEN

!         INITIALIZE WAVENUMBER OUTPUT:
          JFLUX=0
          JOFF=0
          LAST=.TRUE.
          WNCEN=VOUT
          BINMIN=DBLE(VCEN-HBNDWD)
          BINMAX=DBLE(VCEN+HBNDWD)
          CONFLX=1.D0
   10     CONTINUE
          IF(WNCEN.LE.V2)THEN
              IF(RELABS.EQ.'R')THEN

!                 FWHM IS A RELATIVE [%] RESOLUTION.
                  WNFWHM=DBLE(FWHM)*WNCEN/100
                  FWHMSQ=WNFWHM**2
              ELSE

!                 FWHM IS AN ABSOLUTE RESOLUTION.
                  WNFWHM=DBLE(FWHM)
              ENDIF
              WNMIN=WNCEN-WNFWHM
              WNMAX=WNCEN+WNFWHM

!             BRANCH BASED ON SPECTRAL OVERLAP:
              IF(BINMAX.LE.WNMIN)THEN

!                 CASE A:  BIN DOES NOT CONTRIBUTE TO SLIT.
                  RETURN
              ELSEIF(BINMIN.LE.WNMIN)THEN

!                 INTEGRATION BEGINS AT WNMIN:
                  IF(BINMAX.LE.WNCEN)THEN

!                     CASE B:  WNMIN TO BINMAX INTEGRAL.
                      DIFF=BINMAX-WNMIN
                      AVG=DIFF/2
                      WT=DIFF*AVG/FWHMSQ
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSEIF(BINMAX.LT.WNMAX)THEN

!                     CASE C:  WNCEN TO BINMAX INTEGRAL + 1/2.
                      DIFF=BINMAX-WNCEN
                      AVG=WNFWHM-DIFF/2
                      WT=.5D0+DIFF*AVG/FWHMSQ
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSE

!                     CASE D:  WNMIN TO WNMAX INTEGRAL.
                      WT=1.D0
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                      JOFF=JFLUX
                      VOUT=VOUT+DV
                  ENDIF
              ELSEIF(BINMIN.LT.WNCEN)THEN

!                 INTEGRATION BEGINS AT BINMIN:
                  IF(BINMAX.LE.WNCEN)THEN

!                     CASE E:  BINMIN TO BINMAX INTEGRAL - RIGHT SIDE.
                      DIFF=DBLE(BNDWID)
                      AVG=DBLE(VCEN)-WNMIN
                      WT=DIFF*AVG/FWHMSQ
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSEIF(BINMAX.LT.WNMAX)THEN

!                     CASE F:  BINMIN TO BINMAX INTEGRAL - BOTH SIDES.
                      DIFF=WNCEN-BINMIN
                      AVG=WNFWHM-DIFF/2
                      WT=DIFF*AVG/FWHMSQ
                      DIFF=BINMAX-WNCEN
                      AVG=WNFWHM-DIFF/2
                      WT=WT+DIFF*AVG/FWHMSQ
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSE

!                     CASE G:  BINMIN TO WNCEN INTEGRAL + 1/2.
                      DIFF=WNCEN-BINMIN
                      AVG=WNFWHM-DIFF/2
                      WT=.5D0+DIFF*AVG/FWHMSQ
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                      JOFF=JFLUX
                      VOUT=VOUT+DV
                  ENDIF
              ELSEIF(BINMIN.LT.WNMAX)THEN

!                 INTEGRATION BEGINS AT BINMIN:
                  IF(BINMAX.LT.WNMAX)THEN

!                     CASE H:  BINMIN TO BINMAX INTEGRAL - LEFT SIDE.
                      DIFF=DBLE(BNDWID)
                      AVG=WNMAX-DBLE(VCEN)
                      WT=DIFF*AVG/FWHMSQ
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                   ELSE

!                     CASE I:  BINMIN TO WNMAX INTEGRAL
                      DIFF=WNMAX-BINMIN
                      AVG=DIFF/2
                      WT=DIFF*AVG/FWHMSQ
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                      JOFF=JFLUX
                      VOUT=VOUT+DV
                   ENDIF
              ELSE

!                 CASE J:  BIN DOES NOT CONTRIBUTE TO SLIT.
                  RETURN
              ENDIF
              WNCEN=WNCEN+DV
              GOTO 10
          ENDIF
      ELSE

!         INITIALIZE MICRON OR NANOMETER OUTPUT:
          JFLUX=0
          JOFF=0
          LAST=.TRUE.
          VOUTCN=VOUT
          BINMIN=DBLE(VCEN-HBNDWD)
          BINMAX=DBLE(VCEN+HBNDWD)
   20     CONTINUE
          IF(VOUTCN.GE.V1-SMALL*(VOUTCN+V1))THEN
              IF(CHUNIT.EQ.'M')THEN
                  CONVRT=1.D4
              ELSE
                  CONVRT=1.D7
              ENDIF
              IF(RELABS.EQ.'R')THEN

!                 FWHM IS A RELATIVE [%] RESOLUTION.
                  FWHMSQ=(DBLE(FWHM)*VOUTCN/100)**2
                  VOUTMN=(1-DBLE(FWHM)/100)*VOUTCN
                  VOUTMX=(1+DBLE(FWHM)/100)*VOUTCN
              ELSE

!                 FWHM IS A ABSOLUTE.
                  VOUTMN=VOUTCN-DBLE(FWHM)
                  VOUTMX=VOUTCN+DBLE(FWHM)
              ENDIF
              WNCEN=CONVRT/VOUTCN
              WNMIN=CONVRT/VOUTMX
              WNMAX=CONVRT/VOUTMN
              CONFLX=WNCEN/VOUTCN

!             BRANCH BASED ON SPECTRAL OVERLAP:
              IF(BINMAX.LE.WNMIN)THEN

!                 CASE A:  BIN DOES NOT CONTRIBUTE TO SLIT.
                  RETURN
              ELSEIF(BINMIN.LE.WNMIN)THEN

!                 INTEGRATION BEGINS AT WNMIN:
                  IF(BINMAX.LE.WNCEN)THEN

!                     CASE B:  WNMIN TO BINMAX INTEGRAL.
                      DIFF=BINMAX-WNMIN
                      COEF=VOUTMX/BINMAX
                      AVG=(WNMIN+BINMAX)/2
                      WT=WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSEIF(BINMAX.LT.WNMAX)THEN

!                     CASE C:  WNCEN TO BINMAX INTEGRAL + 1/2.
                      DIFF=BINMAX-WNCEN
                      COEF=VOUTCN/BINMAX
                      AVG=(WNCEN+BINMAX)/2
                      WT=.5D0+WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSE

!                     CASE D:  WNMIN TO WNMAX INTEGRAL.
                      WT=1.D0
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                      JOFF=JFLUX
                      VOUT=VOUT-DV
                  ENDIF
              ELSEIF(BINMIN.LT.WNCEN)THEN

!                 INTEGRATION BEGINS AT BINMIN:
                  IF(BINMAX.LE.WNCEN)THEN

!                     CASE E:  BINMIN TO BINMAX INTEGRAL - RIGHT SIDE.
                      DIFF=DBLE(BNDWID)
                      COEF=CONVRT/(BINMIN*BINMAX)
                      AVG=DBLE(VCEN)
                      WT=WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSEIF(BINMAX.LT.WNMAX)THEN

!                     CASE F:  BINMIN TO BINMAX INTEGRAL - BOTH SIDES.
                      DIFF=WNCEN-BINMIN
                      COEF=VOUTCN/BINMIN
                      AVG=(BINMIN+WNCEN)/2
                      WT=WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                      DIFF=BINMAX-WNCEN
                      COEF=VOUTCN/BINMAX
                      AVG=(WNCEN+BINMAX)/2
                      WT=WT+WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSE

!                     CASE G:  BINMIN TO WNCEN INTEGRAL + 1/2.
                      DIFF=WNCEN-BINMIN
                      COEF=VOUTCN/BINMIN
                      AVG=(BINMIN+WNCEN)/2
                      WT=.5D0+WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                      JOFF=JFLUX
                      VOUT=VOUT-DV
                  ENDIF
              ELSEIF(BINMIN.LT.WNMAX)THEN

!                 INTEGRATION BEGINS AT BINMIN:
                  IF(BINMAX.LT.WNMAX)THEN

!                     CASE H:  BINMIN TO BINMAX INTEGRAL - LEFT SIDE.
                      DIFF=DBLE(BNDWID)
                      COEF=CONVRT/(BINMIN*BINMAX)
                      AVG=DBLE(VCEN)
                      WT=WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                      LAST=.FALSE.
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                  ELSE

!                     CASE I:  BINMIN TO WNMAX INTEGRAL
                      DIFF=WNMAX-BINMIN
                      COEF=VOUTMN/BINMIN
                      AVG=(BINMIN+WNMAX)/2
                      WT=WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                      CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)
                      JOFF=JFLUX
                      VOUT=VOUT-DV
                  ENDIF
              ELSE

!                 CASE J:  BIN DOES NOT CONTRIBUTE TO SLIT.
                  RETURN
              ENDIF
              VOUTCN=VOUTCN-DV
              GOTO 20
          ENDIF
      ENDIF
      RETURN
      END
