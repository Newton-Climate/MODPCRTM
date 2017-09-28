      SUBROUTINE WTSUM(JFLUX,JOFF,CONFLX,WT,LAST)

!     WTSUM WRITES OUT SPECTRAL FLUX, DEGRADED WITH
!     A TRIANGULAR SLIT FUNCTION, TO UNIT IFLUX.

!     ARGUMENTS:
!       JFLUX    SPECTRAL INDEX FOR OUTPUT FLUX VALUES.
!       JOFF     SPECTRAL INDEX OFFSET (EQUALS NUMBER OF SLIT
!                FUNCTIONS TERMINATING IN CURRENT SPECTRAL BIN).
!       CONFLX   FLUX SPECTRAL CONVERSION FACTOR [CM-1 / SPECTRAL UNIT].
!       WT       WEIGHT OF SPECTRAL BIN IN EVALUATING SLIT FUNCTION.
!       LAST     LOGICAL FLAG, TRUE IF CURRENT SLIT FUNCTION
!                TERMINATES IN CURRENT SPECTRAL BIN.
      INTEGER JFLUX,JOFF
      LOGICAL LAST
      DOUBLE PRECISION CONFLX,WT

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

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
!        FRMT    FORMAT USED IN FLUX TABLE.
      CHARACTER FRMT*50
      COMMON/WTFLXC/FRMT

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

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       IK       LAYER INDEX.
!       KFLUX    SPECTRAL SLIT INDEX.
!       CONWT    PRODUCT OF FLUX SPECTRAL CONVERSION AND SLIT
!                FUNCTION WEIGHTING FACTORS [CM-1 / SPECTRAL UNIT].
!       UDIFF    BOUNDARY UPWARD DIFFUSE SPECTRAL FLUX [W CM-2 / CM-1].
!       DDIFF    BOUNDARY DOWNWARD DIFFUSE SPECTRAL FLUX [W CM-2 / CM-1]
!       DDRCT    BOUNDARY DIRECT SOLAR SPECTRAL FLUX [W CM-2 / CM-1].
      INTEGER IK,KFLUX
      DOUBLE PRECISION CONWT,UDIFF(LAYDIM),DDIFF(LAYDIM),DDRCT(LAYDIM)

!     COMBINE FLUX SPECTRAL CONVERSION AND WEIGHTING FACTOR.
      CONWT=CONFLX*WT
      IF(LAST)THEN

!         WRITE OUT FLUXES:
          IF(JFLUX.GE.MWGT)THEN
              WRITE(IPR,'(/A,F8.1,A,/17X,A,I3,A,/17X,2A)')              &
     &          ' Error in WTSUM:  The spectral bin at',VBAND,          &
     &          ' CM-1 contributes to',' more than MWGT (=',MWGT,       &
     &          ') spectral slit functions.',' Increase parameter',     &
     &          ' MWGT or modify CARD4 spectral inputs.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'Error in WTSUM:  See tape6 (or .tp6) for details.'
          ENDIF
          JFLUX=JFLUX+1
          DO IK=1,ML
              UDIFF(IK)=UPDIFF(IK,JFLUX)+CONWT*UPDIFF(IK,-1)
              DDIFF(IK)=DNDIFF(IK,JFLUX)+CONWT*DNDIFF(IK,-1)
              DDRCT(IK)=DNDRCT(IK,JFLUX)+CONWT*DNDRCT(IK,-1)
          ENDDO
          IF(.NOT.LJMASS)THEN
              IF(BINOUT)THEN
                  CALL BNFLX4(1.D0,VOUT,                                &
     &              UDIFF,DDIFF,DDRCT,LAYDIM,ML,MLFLX)
              ELSE
                  WRITE(IFLUX,FMT=FRMT)VOUT,(UDIFF(IK),DDIFF(IK),       &
     &              DDRCT(IK),IK=1,MLFLX),UDIFF(ML),DDIFF(ML),DDRCT(ML)
              ENDIF
          ENDIF
          NTERMS=NTERMS+1
          DO IK=1,ML
              SMUPDF(IK)=SMUPDF(IK)+UDIFF(IK)
              SMDNDF(IK)=SMDNDF(IK)+DDIFF(IK)
              SMDNDR(IK)=SMDNDR(IK)+DDRCT(IK)
              UPDIFF(IK,JFLUX)=0.D0
              DNDIFF(IK,JFLUX)=0.D0
              DNDRCT(IK,JFLUX)=0.D0
          ENDDO


      ELSEIF(JOFF.EQ.0)THEN

!         INCREMENT WEIGHTED SUMS:
          IF(JFLUX.GE.MWGT)THEN
              WRITE(IPR,'(/A,F8.1,A,/17X,A,I3,A,/17X,2A)')              &
     &          ' Error in WTSUM:  The spectral bin at',VBAND,          &
     &          ' CM-1 contributes to',' more than MWGT (=',MWGT,       &
     &          ') spectral slit functions.',' Increase parameter',     &
     &          ' MWGT or modify CARD4 spectral inputs.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'Error in WTSUM:  See tape6 (or .tp6) for details.'
          ENDIF
          JFLUX=JFLUX+1
          DO IK=1,ML
              UPDIFF(IK,JFLUX)=UPDIFF(IK,JFLUX)+CONWT*UPDIFF(IK,-1)
              DNDIFF(IK,JFLUX)=DNDIFF(IK,JFLUX)+CONWT*DNDIFF(IK,-1)
              DNDRCT(IK,JFLUX)=DNDRCT(IK,JFLUX)+CONWT*DNDRCT(IK,-1)
          ENDDO
      ELSEIF(JFLUX.GE.MWGT)THEN

!         TRANSLATE AND INITIALIZE WEIGHTED SUMS:
          JFLUX=JFLUX+1-JOFF
          JOFF=0
          DO IK=1,ML
              UPDIFF(IK,JFLUX)=CONWT*UPDIFF(IK,-1)
              DNDIFF(IK,JFLUX)=CONWT*DNDIFF(IK,-1)
              DNDRCT(IK,JFLUX)=CONWT*DNDRCT(IK,-1)
              DO KFLUX=JFLUX+1,MWGT
                  UPDIFF(IK,KFLUX)=0.D0
                  DNDIFF(IK,KFLUX)=0.D0
                  DNDRCT(IK,KFLUX)=0.D0
              ENDDO
          ENDDO
      ELSE

!         TRANSLATE, INCREMENT, AND REINITIALIZE "JFLUX" WEIGHTED SUMS:
          JFLUX=JFLUX+1
          DO IK=1,ML
              UPDIFF(IK,JFLUX-JOFF)=UPDIFF(IK,JFLUX)+CONWT*UPDIFF(IK,-1)
              DNDIFF(IK,JFLUX-JOFF)=DNDIFF(IK,JFLUX)+CONWT*DNDIFF(IK,-1)
              DNDRCT(IK,JFLUX-JOFF)=DNDRCT(IK,JFLUX)+CONWT*DNDRCT(IK,-1)
              UPDIFF(IK,JFLUX)=0.D0
              DNDIFF(IK,JFLUX)=0.D0
              DNDRCT(IK,JFLUX)=0.D0
          ENDDO
      ENDIF
      RETURN
      END
