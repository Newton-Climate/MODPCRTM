      LOGICAL FUNCTION WTCHAN(LSTBIN,SPC_LO,WGTLO,SPC_UP,WGTUP)

!     WTCHAN DETERMINES THE RESPONSE FUNCTION CONTRIBUTIONS (WTCHN)
!     TO BAND MODEL SPECTRAL BINS.  A VALUE OF FALSE IS RETURNED
!     IF THE CALCULATION IS NOT SUCCESSFUL.  THE CONTRIBUTIONS ARE
!     CALCULATED AS FREQUENCY INTEGRALS FROM V0 TO V1 OVER A
!     TABULATED CHANNEL RESPONSE FUNCTION, W(V):

!                   / V1
!         WTCHN  =  |      W(V)  DV
!                   / V0

!     IF SPC_LO AND SPC_UP ARE IN WAVENUMBERS [CM-1],

!                          V - SPC_LO
!        W(V) = WGTLO + --------------- (WGTUP - WGTLO)   .
!                       SPC_UP - SPC_LO

!     IF SPC_LO AND SPC_UP ARE IN WAVELENGTH UNITS,

!                       CONVRT/V - SPC_UP
!        W(V) = WGTUP + ----------------- (WGTUP - WGTLO)   .
!                        SPC_LO - SPC_UP

!     WHERE CONVRT IS THE CONVERSION FACTOR BETWEEN THE WAVELENGTH
!     UNIT AND WAVENUMBERS.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'CHANLS.h'

!     ARGUMENTS:
!       LSTBIN   LAST WAVENUMBER SPECTRAL BIN FOR WHICH A CHANNEL
!                CONTRIBUTION WAS CALCULATED [=-1 FOR NEW CHANNELS].
!       SPC_LO   SPECTRAL GRID LOWER VALUE [UNITS FROM UNTFLG].
!       WGTLO    FILTER RESPONSE WEIGHT AT SPC_LO.
!       SPC_UP   SPECTRAL GRID UPPER VALUE [UNITS FROM UNTFLG].
!       WGTUP    FILTER RESPONSE WEIGHT AT SPC_UP.
        INTEGER LSTBIN
        REAL SPC_LO,WGTLO,SPC_UP,WGTUP

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     /RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /JM1A1/
!       RESCHR   BAND MODEL RESOLUTION CHARACTER STRING.
!       DISSTR   CHARACTER STRING USED TO READ IN DISORT LOGICALS.
!       H2OSTR   VERTICAL WATER COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE WATER
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE WATER COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS A WATER COLUMN SCALING FACTOR).
!       O3STR    VERTICAL OZONE COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE OZONE
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE OZONE COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS AN OZONE COLUMN SCALING FACTOR).
!       USRSUN   USER-SPECIFIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       BMROOT   PREFIX OF MOLECULAR BAND MODEL PARAMETERS FILE.
!       FILTNM   NAME OF FILTER RESPONSE FUNCTION FILE.
!       H2OAER   FLAG, TRUE IF DEFAULT AEROSOL PROPERTIES ARE REVISED
!                BASED ON WATER COLUMN SCALING.
!       DATDIR   NAME OF THE MODTRAN DATA DIRECTORY.
      CHARACTER RESCHR*2,DISSTR*3,H2OSTR*10,O3STR*10,USRSUN*(NAMLEN),   &
     &  FILTNM*(NAMLEN),BMROOT*(NAMLEN),H2OAER*1,DATDIR*(NAMLEN-LENSUN)
      COMMON/JM1A1/RESCHR,DISSTR,H2OSTR,O3STR,USRSUN,                   &
     &  FILTNM,BMROOT,H2OAER,DATDIR

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       NBINLO   LOWER MOST BIN NUMBER OVERLAPPING SPECTRAL CHANNEL.
!       NBINUP   UPPER MOST BIN NUMBER OVERLAPPING SPECTRAL CHANNEL.
!       NBIN     CURRENT BAND MODEL BIN NUMBER.
!       SPEC_L   SPECTRAL GRID LOWER VALUE [UNITS FROM UNTFLG].
!       WGT_LO   FILTER RESPONSE WEIGHT AT SPEC_L.
!       SPEC_U   SPECTRAL GRID UPPER VALUE [UNITS FROM UNTFLG].
!       WGT_UP   FILTER RESPONSE WEIGHT AT SPEC_U.
!       BANDWD   SPECTRAL WIDTH OF BAND MODEL [CM-1].
!       WGT      BAND MODEL BIN CONTRIBUTION TO WTCHN ARRAY.
!       DFREQ    SPECTRAL FREQUENCY RANGE OF CURRENT "WGT" CONTRIBUTION.
!       WN_RGT   RIGHT SIDE EDGE OF CURRENT SPECTRAL BIN [CM-1].
!       DWGT     WGT_UP MINUS WGT_LO.
!       DSPEC    INTEGRAL RANGE (SPEC_L - SPEC_U) [UNITS FROM UNTFLG].
!       SLOPE    DWGT OVER DSPEC [UNITS FROM UNTFLG].
!       BN_CEN   CENTER OF BAND MODEL BIN IS (NBIN+BN_CEN)*BANDWD.
!       BN_RGT   RIGHT EDGE OF BAND MODEL BIN IS (NBIN+BN_RGT)*BANDWD.
!       COEFLN   COEFFICIENT OF NATURAL LOGARITHM.
!       SLOPWD   PRODUCT OF SLOPE AND BAND WIDTH (BANDWD).
!       CONVRT   CONVERSION FACTOR BETWEEN WAVELENGTH AND FREQUENCY.
!       DRATIO   SPECTRAL DIFFERENCE OVER SPECTRAL VALUE.
      INTEGER NBINLO,NBINUP,NBIN
      DOUBLE PRECISION SPEC_L,WGT_LO,SPEC_U,WGT_UP,BANDWD,WGT,WN_LO,    &
     &  WN_UP,DFREQ,WN_RGT,DWGT,DSPEC,SLOPE,COEFLN,WGT0,WGT0WD,SLOPWD,  &
     &  CONVRT,DRATIO,BN_CEN,BN_RGT

!     INTERNAL FUNCTIONS:
!       EXPAND   SMALL X EXPANSION OF   - [ 1 + ln( 1 - X ) / X ] / X.
!       X        ARGUMENT OF FUNCTION EXPAND.
      DOUBLE PRECISION EXPAND,X
      EXPAND(X)=(3+X*(2+X*(1.5D0+X*(1.2D0+X))))/6

!     SWITCH TO DOUBLE PRECISION VARIABLES:
      SPEC_L=DBLE(SPC_LO)
      WGT_LO=DBLE(WGTLO)
      SPEC_U=DBLE(SPC_UP)
      WGT_UP=DBLE(WGTUP)
      BANDWD=DBLE(BNDWID)

!     INITIALIZE WTCHAN TO .FALSE. AND CHECK INPUTS:
      WTCHAN=.FALSE.
      IF(SPEC_L.GT.SPEC_U)THEN
          WRITE(IPR,'(/1X,A)')                                          &
     &      'WARNING:  SPC_UP is less than SPC_LO in function WTCHAN'
          RETURN
      ELSEIF(WGT_LO.LT.0.D0 .OR. WGT_UP.LT.0.D0)THEN
          WRITE(IPR,'(/1X,A)')                                          &
     &      'WARNING:  Negative weight input to function WTCHAN'
          RETURN
      ENDIF

!     BRANCH BASED ON UNTFLG:
      IF(UNTFLG.EQ.'W')THEN

!         FREQUENCY GRID [CM-1] - DETERMINE MINIMUM
!         AND MAXIMUM SPECTRAL WAVENUMBER BINS:
          IF(RESCHR.EQ.'p1')THEN
              NBINLO=INT(SPEC_L/BANDWD)
              NBINUP=INT(SPEC_U/BANDWD)
              BN_CEN=.5D0
              BN_RGT=1.D0
          ELSE
              NBINLO=NINT(SPEC_L/BANDWD)
              NBINUP=NINT(SPEC_U/BANDWD)
              BN_CEN=.0D0
              BN_RGT=.5D0
          ENDIF
          IF(NBINLO.LT.MNBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MNBIN (=',MNBIN,') must be decreased.',     &
     &          ' Current NBINLO in function WTCHAN equals',NBINLO
              RETURN
          ELSEIF(NBINUP.GT.MXBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MXBIN (=',MXBIN,') must be increased.',     &
     &          ' Current NBINUP in function WTCHAN equals',NBINUP
              RETURN
          ENDIF

!         INCREMENT INTEGRATED WIDTH OF CHANNEL:
          WDFREQ(NCHAN)                                                 &
     &      =WDFREQ(NCHAN)+SNGL((WGT_UP+WGT_LO)*(SPEC_U-SPEC_L)/2)
          SLOPE=(WGT_UP-WGT_LO)/(SPEC_U-SPEC_L)
          IF(SPEC_L.GT.0.D0)THEN
              WDWAVE(NCHAN)=WDWAVE(NCHAN)+SNGL(10000*                   &
     &          (WGT_LO/SPEC_L-WGT_UP/SPEC_U+SLOPE*LOG(SPEC_U/SPEC_L)))
          ELSE
              WDWAVE(NCHAN)=BIGNUM
          ENDIF

!         DETERMINE CONTRIBUTIONS TO EACH SPECTRAL BIN.
          IF(NBINLO.EQ.NBINUP)THEN

!             ONE BIN ONLY; CALCULATE WEIGHT.
              WGT=(WGT_LO+WGT_UP)*(SPEC_U-SPEC_L)/2
              IF(NBINLO.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINLO),NBINLO)=                         &
     &              WTCHN(NUMCHN(NBINLO),NBINLO)+SNGL(WGT)
              ELSEIF(NUMCHN(NBINLO).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINLO CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINLO)=NUMCHN(NBINLO)+1
                  LSTCHN(NUMCHN(NBINLO),NBINLO)=NCHAN
                  WTCHN(NUMCHN(NBINLO),NBINLO)=SNGL(WGT)
              ENDIF
          ELSE

!             MULTIPLE BINS; DEFINE FIRST WEIGHT.
              DFREQ=BANDWD*(NBINLO+BN_RGT)-SPEC_L
              WGT=DFREQ*(WGT_LO+DFREQ*SLOPE/2)
              IF(NBINLO.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINLO),NBINLO)=                         &
     &              WTCHN(NUMCHN(NBINLO),NBINLO)+SNGL(WGT)
              ELSEIF(NUMCHN(NBINLO).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINLO CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINLO)=NUMCHN(NBINLO)+1
                  LSTCHN(NUMCHN(NBINLO),NBINLO)=NCHAN
                  WTCHN(NUMCHN(NBINLO),NBINLO)=SNGL(WGT)
              ENDIF

!             LOOP OVER INTERMEDIATE BINS (IF ANY):
              SLOPWD=BANDWD*SLOPE
              WGT0=BANDWD*WGT_LO
              DO NBIN=NBINLO+1,NBINUP-1
                  IF(NUMCHN(NBIN).GE.MXNCHN)THEN

!                     SPECTRAL BIN NBIN CONTRIBUTES TO TOO MANY CHANNELS
                      WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',      &
     &                  ' Parameter MXNCHN must be increased. ',        &
     &                  ' Single spectral',' bins contribute',          &
     &                  ' to more than',MXNCHN,' filter channels.'
                      RETURN
                  ENDIF

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBIN)=NUMCHN(NBIN)+1
                  LSTCHN(NUMCHN(NBIN),NBIN)=NCHAN
                  WTCHN(NUMCHN(NBIN),NBIN)                              &
     &              =SNGL(WGT0+(BANDWD*(NBIN+BN_CEN)-SPEC_L)*SLOPWD)
              ENDDO

!             LAST BIN:  INCREMENT COUNTER, LIST CHANNEL
!             AND DEFINE WEIGHT AND VARIABLE LSTBIN.
              IF(NUMCHN(NBINUP).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINUP CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ENDIF
              DFREQ=SPEC_U-(NBINUP+BN_RGT-1)*BANDWD
              NUMCHN(NBINUP)=NUMCHN(NBINUP)+1
              LSTCHN(NUMCHN(NBINUP),NBINUP)=NCHAN
              WTCHN(NUMCHN(NBINUP),NBINUP)                              &
     &          =SNGL(DFREQ*(WGT_UP-DFREQ*SLOPE/2))
          ENDIF

!         SET LSTBIN TO THE MAXIMUM BIN FOR FREQUENCY FILTERS.
          LSTBIN=NBINUP
      ELSE

!         WAVELENGTH GRID: NANOMETERS OR MICRONS?
          IF(UNTFLG.EQ.'N')THEN

!             NANOMETERS:
              CONVRT=1.D7
          ELSEIF(UNTFLG.EQ.'M')THEN

!             MICRONS:
              CONVRT=1.D4
          ELSE
              WRITE(IPR,'(/2A,2(/11X,A),2A)')' WARNING:  ',             &
     &          'Input UNTFLG to routine WTCHAN must equal "W" for',    &
     &          'WAVENUMBERS, "M" for Microns, or "N" for Nanometers.', &
     &          'Currently, UNTFLG = "',UNTFLG,'"'
              RETURN
          ENDIF

!         DETERMINE BAND MODEL BIN RANGE:
          WN_LO=CONVRT/SPEC_U
          WN_UP=CONVRT/SPEC_L
          IF(RESCHR.EQ.'p1')THEN
              NBINLO=INT(WN_LO/BANDWD)
              NBINUP=INT(WN_UP/BANDWD)
              BN_RGT=1.D0
          ELSE
              NBINLO=NINT(WN_LO/BANDWD)
              NBINUP=NINT(WN_UP/BANDWD)
              BN_RGT=.5D0
          ENDIF
          IF(NBINLO.LT.MNBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MNBIN (=',MNBIN,') must be decreased.',     &
     &          ' Current NBINLO in function WTCHAN equals',NBINLO
              RETURN
          ELSEIF(NBINUP.GT.MXBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MXBIN (=',MXBIN,') must be increased.',     &
     &          ' Current NBINUP in function WTCHAN equals',NBINUP
              RETURN
          ENDIF

!         INCREMENT INTEGRATED WIDTH OF CHANNEL:
          WDWAVE(NCHAN)                                                 &
     &      =WDWAVE(NCHAN)+SNGL((WGT_UP+WGT_LO)*(SPEC_U-SPEC_L)/2)

!         THE WEIGHTING FUNCTION SLOPE:
          DSPEC=SPEC_U-SPEC_L
          DWGT=WGT_UP-WGT_LO
          SLOPE=DWGT/DSPEC
          COEFLN=CONVRT*SLOPE
          IF(NBINLO.EQ.NBINUP)THEN

!             ONE BIN ONLY; CALCULATE WEIGHT.
              DRATIO=DSPEC/SPEC_U
              IF(DRATIO.GT..01D0)THEN
                  WGT=(WN_UP*WGT_LO-WN_LO*WGT_UP)-COEFLN*LOG(1-DRATIO)
              ELSE
                  WGT=DRATIO*(WGT_LO*WN_UP+DWGT*WN_LO*EXPAND(DRATIO))
              ENDIF
              IF(NBINLO.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINLO),NBINLO)=                         &
     &              WTCHN(NUMCHN(NBINLO),NBINLO)+SNGL(WGT)
              ELSEIF(NUMCHN(NBINLO).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINLO CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINLO)=NUMCHN(NBINLO)+1
                  LSTCHN(NUMCHN(NBINLO),NBINLO)=NCHAN
                  WTCHN(NUMCHN(NBINLO),NBINLO)=SNGL(WGT)
              ENDIF
              WDFREQ(NCHAN)=WDFREQ(NCHAN)+SNGL(WGT)
          ELSE

!             MULTIPLE BINS; DEFINE FIRST WEIGHT.
              WN_RGT=(NBINLO+BN_RGT)*BANDWD
              DFREQ=WN_RGT-WN_LO
              DRATIO=DFREQ/WN_LO
              WGT0=(WGT_LO*SPEC_U-WGT_UP*SPEC_L)/DSPEC
              IF(DRATIO.GT..01D0)THEN
                  WGT=DFREQ*WGT0+COEFLN*LOG(1+DRATIO)
              ELSE
                  WGT=DFREQ*(WGT_UP-SPEC_U*SLOPE*DRATIO*EXPAND(-DRATIO))
              ENDIF
              IF(NUMCHN(NBINLO).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINLO CONTRIBUTES TO TOO MANY CHANNELS:
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ENDIF

!             INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
              NUMCHN(NBINLO)=NUMCHN(NBINLO)+1
              LSTCHN(NUMCHN(NBINLO),NBINLO)=NCHAN
              WTCHN(NUMCHN(NBINLO),NBINLO)=SNGL(WGT)
              WDFREQ(NCHAN)=WDFREQ(NCHAN)+SNGL(WGT)

!             LOOP OVER INTERMEDIATE BINS:
              WGT0WD=BANDWD*WGT0
              DO NBIN=NBINLO+1,NBINUP-1
                  IF(NUMCHN(NBIN).GE.MXNCHN)THEN

!                     SPECTRAL BIN NBIN CONTRIBUTES TO TOO MANY CHANNELS
                      WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',      &
     &                  ' Parameter MXNCHN must be increased. ',        &
     &                  ' Single spectral',' bins contribute',          &
     &                  ' to more than',MXNCHN,' filter channels.'
                      RETURN
                  ENDIF

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBIN)=NUMCHN(NBIN)+1
                  LSTCHN(NUMCHN(NBIN),NBIN)=NCHAN
                  WN_RGT=(NBIN+BN_RGT)*BANDWD
                  DRATIO=BANDWD/WN_RGT
                  IF(DRATIO.GT..01D0)THEN
                      WGT=WGT0WD-COEFLN*LOG(1-DRATIO)
                  ELSE
                      WGT=WGT0WD+COEFLN*DRATIO*(1+DRATIO*EXPAND(DRATIO))
                  ENDIF
                  WTCHN(NUMCHN(NBIN),NBIN)=SNGL(WGT)
                  WDFREQ(NCHAN)=WDFREQ(NCHAN)+SNGL(WGT)
              ENDDO

!             LAST BIN:  INCREMENT COUNTER, LIST CHANNEL
!             AND DEFINE WEIGHT AND VARIABLE LSTBIN.
              DFREQ=WN_UP-WN_RGT
              DRATIO=DFREQ/WN_UP
              IF(DRATIO.GT..01)THEN
                  WGT=DFREQ*WGT0-COEFLN*LOG(1-DRATIO)
              ELSE
                  WGT=DFREQ*(WGT_LO+SPEC_L*SLOPE*DRATIO*EXPAND(DRATIO))
              ENDIF
              IF(NBINUP.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINUP),NBINUP)=                         &
     &              WTCHN(NUMCHN(NBINUP),NBINUP)+SNGL(WGT)
              ELSEIF(NUMCHN(NBINUP).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINUP CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINUP)=NUMCHN(NBINUP)+1
                  LSTCHN(NUMCHN(NBINUP),NBINUP)=NCHAN
                  WTCHN(NUMCHN(NBINUP),NBINUP)=SNGL(WGT)
              ENDIF
              WDFREQ(NCHAN)=WDFREQ(NCHAN)+SNGL(WGT)
          ENDIF

!         SET LSTBIN TO THE MINIMUM BIN FOR WAVELENGTH FILTERS.
          LSTBIN=NBINLO
      ENDIF

!     CHANGE WTCHAN TO .TRUE. AND RETURN:
      WTCHAN=.TRUE.
      RETURN
      END
