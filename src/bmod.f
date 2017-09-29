      SUBROUTINE BMOD(ILOS,ILAY,NSEG,IPATH,IBINPT,LBMCEN,LBMFLG)

!     BMOD RETURNS BAND MODEL MOLECULAR TRANSMITTANCES:
!       IMOL = 1      H2O LINE (CENTERS, TAILS)
!       IMOL = 2      CO2 LINE (CENTERS, TAILS)
!       IMOL = 3       O3 LINE (CENTERS, TAILS)
!       IMOL = 4      N2O LINE (CENTERS, TAILS)
!       IMOL = 5       CO LINE (CENTERS, TAILS)
!       IMOL = 6      CH4 LINE (CENTERS, TAILS)
!       IMOL = 7       O2 LINE (CENTERS, TAILS)
!       IMOL = 8       NO LINE (CENTERS, TAILS)
!       IMOL = 9      SO2 LINE (CENTERS, TAILS)
!       IMOL =10      NO2 LINE (CENTERS, TAILS)
!       IMOL =11      NH3 LINE (CENTERS, TAILS)
!       IMOL =12     HNO3 LINE (CENTERS, TAILS)
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ONLY ARGUMENTS:
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
!       ILAY     LAYER INDEX.
!       NSEG     NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       IPATH    PATH TYPE SWITCH.
!       IBINPT   BIN NUMBER OF CURRENT SPECTRAL POINT.
!                (CENTER FREQUENCY = IBINPT * BNDWID + OSHIFT).
      INTEGER ILOS,ILAY,NSEG,IPATH,IBINPT

!     INPUT/OUTPUT ARGUMENTS:
!       LBMCEN   FLAG, TRUE IF S/d & 1/d DATA FOR CURRENT FREQUENCY BIN.
!       LBMFLG   FLAG, TRUE IF CENTER OR TAIL DATA FOR CURRENT FREQ BIN.
      LOGICAL LBMCEN,LBMFLG

!     COMMONS:
      INCLUDE 'YPROP.h'
      INCLUDE 'BMHEAD.h'
      INCLUDE 'BMPATH.h'
      INCLUDE 'BMPTHS.h'
      INCLUDE 'BASE.h'
      INCLUDE 'SEGDAT.h'
      INCLUDE 'WSOL.h'

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

!     /JM5/
!       IRPT     REPEAT INPUT FLAG (0=NONE, 1=ALL, 3=GEOM, 4=SPEC).
!       IFAC     CURRENT COLUMN SCALING FACTOR INDEX.
!       NFACMN   NUMBER OF COLUMN SCALING FACTOR LESS THAN 1.
!       NFACMX   NUMBER OF COLUMN SCALING FACTOR GREATER THAN 1.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       SCALMN   MINIMUM COLUMN SCALING FACTOR.
!       SCALMX   MAXIMUM COLUMN SCALING FACTOR.
!       LBMRES   LOGICAL FLAG, .TRUE. FOR BAND MODEL AND .FALSE.
!                FOR 1 NM SPECTRAL RESOLUTION OUTPUT.
      INTEGER IRPT,IFAC,NFACMN,NFACMX
      REAL FACMC
      DOUBLE PRECISION SCALMN,SCALMX
      LOGICAL LBMRES
      COMMON/JM5/SCALMN,SCALMX,IRPT,IFAC,NFACMN,NFACMX,FACMC,LBMRES

!     /ACTIVE/
!       NACT     NUMBER OF ACTIVE MOLECULES FOR CURRENT FREQUENCY BIN.
!       NACTBM   NUMBER OF ACTIVE BAND MODEL MOLECULES FOR CURRENT FREQ.
!       NACTX    NUMBER OF ACTIVE X & Y CROSS-SECTION MOLECULES.
!       MACTBM   LIST OF ACTIVE BAND MODEL MOLECULES FOR CURRENT FREQ.
!       MACTX    LIST OF ACTIVE X & Y CROSS-SECTION MOLECULES.
      INTEGER NACT,NACTBM,NACTX,MACTBM,MACTX
      COMMON/ACTIVE/NACT,NACTBM,NACTX,MACTBM(MMOLYT),MACTX(NMOLX)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       IMOL     MOLECULAR INDEX.
!       LMAX     LAYER LOOP MAXIMUM INDEX.
!       LSOLAR   SOLAR PATH SEGMENT INDEX.
      INTEGER IMOL,LMAX,LSOLAR

!     SAVED DATA:
      INTEGER IPATH0
      SAVE IPATH0
      DATA IPATH0/0/

!     INITIALIZE TRANSMITTANCES AND ABSORPTION COEFFICIENTS SKIPPING
!     H2-H2 & H2-HE DIMERS:
      DO IMOL=1,NMOLXT-2
          TX(MPOINT(IMOL))=1.
      ENDDO
      DO IMOL=NMOLXT+1,NMOLT
          TX(MPOINT(IMOL))=1.
      ENDDO

!     DOES IBINPT EXCEED BAND MODEL DATA RANGE?
      IF(IBINPT.GT.MBINBM)THEN

!         SET NUMBER OF ACTIVE MOLECULES TO ZERO AND RETURN:
          NACTBM=0
          NACTX=0
          LBMCEN=.FALSE.
          LBMFLG=.FALSE.
          RETURN
      ENDIF

!     RETURN TO LOOP IF NO ACTIVE MOLECULES:
      IF(NACT.EQ.0)RETURN

!     DETERMINE PATH TYPE:
      IF(ILOS.EQ.MLOSP1)THEN

!         MULTIPLE SCATTERING VERTICAL PATH:
          IF(IEMSCT.EQ.1)THEN

!             LOOP OVER SINGLE LAYER FOR THERMAL RADIANCE CALCULATION:
              LMAX=ILAY
          ELSEIF(IPATH.LT.3)THEN

!             IEMSCT=2 & (IPATH=1 OR IPATH=2): SKIP LAYER LOOP.
              LMAX=0
              IPATH0=IPATH
              LSOLAR=ILAY-1+IPATH
          ELSE

!             IEMSCT=2 & IPATH=3: LOOP OVER SINGLE LAYER.
              LMAX=ILAY
              IPATH0=IPATH
          ENDIF

!         START BAND MODEL MOLECULE LOOP:
          CALL BMLOOP(IEMSCT,IPATH,ILAY,LMAX,LSOLAR,.TRUE.,             &
     &      FACMC,JTMS(1),KTMS(1),FTMS(1),GTMS(1),WPTHMS(0,1),          &
     &      T5MS(1),PTMS(1),PATMMS(1),P2MS(1),JTMSS(1,1),               &
     &      KTMSS(1,1),FTMSS(1,1),GTMSS(1,1),WSPTHM(0,1),               &
     &      T5MSS(1,1),PTMSS(1,1),PSATMM(1,1),P2MSS(1,1))
      ELSE

!         LINE-OF-SIGHT PATH:
          IF(IEMSCT.EQ.0 .OR. IEMSCT.EQ.3)THEN

!             LOOP OVER ALL LAYERS FOR TRANSMITTANCE CALCULATIONS:
              LMAX=NSEG
          ELSEIF(IEMSCT.EQ.1)THEN

!             LOOP OVER SINGLE LAYER FOR THERMAL RADIANCE CALCULATION:
              LMAX=ILAY
          ELSEIF(IPATH.EQ.1)THEN

!             IEMSCT=2 & IPATH=1 (FIRST SOLAR PATH): SKIP LAYER LOOP.
              LMAX=0
              LSOLAR=1
              IPATH0=IPATH
          ELSEIF(IPATH.EQ.2)THEN

!             IEMSCT=2 & IPATH=2: SINGLE LAYER ONLY FOR "L-PATH".
              LMAX=ILAY
              LSOLAR=ILAY+1
              IPATH0=IPATH
          ELSEIF(IPATH0.EQ.2)THEN

!             IEMSCT=2 & IPATH=3: SKIP SINGLE LAYER (DONE WHEN IPATH=2).
              LMAX=0
              IPATH0=IPATH
          ELSE

!             IEMSCT=2 & IPATH=3: SINGLE LAYER (SHADE, IPATH=2 SKIPPED).
              LMAX=ILAY
              IPATH0=IPATH
          ENDIF

!         START BAND MODEL MOLECULE LOOP:
          CALL BMLOOP(IEMSCT,IPATH,ILAY,LMAX,LSOLAR,.FALSE.,FACMC,      &
     &      JTLOS(1,ILOS),KTLOS(1,ILOS),FTLOS(1,ILOS),GTLOS(1,ILOS),    &
     &      WPTH(0,1,ILOS),T5LOS(1,ILOS),PTLOS(1,ILOS),PATM(1,ILOS),    &
     &      P2LOS(1,ILOS),JTLOSS(1,1,ILOS),KTLOSS(1,1,ILOS),            &
     &      FTLOSS(1,1,ILOS),GTLOSS(1,1,ILOS),WSPTH(0,1,ILOS),          &
     &      T5LOSS(1,1,ILOS),PTLOSS(1,1,ILOS),                          &
     &      PSATM(1,1,ILOS),P2LOSS(1,1,ILOS))
      ENDIF

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE BMLOOP(IEMSCT,IPATH,ILAY,LMAX,LSOLAR,                  &
     &  MSPATH,FACMC,JT,KT,FT,GT,WPTH,T5,PTM75,PATM,P2,                 &
     &  JTS,KTS,FTS,GTS,WSPTH,T5S,PTM75S,PSATM,P2S)

!     PERFORMS BAND MODEL SPECIES LOOP OVER ACTIVE MOLECULES:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       IEMSCT   RADIATIVE TRANSFER MODE (0=TRANSMITTANCE, 1=THERMAL,
!                  2=THERMAL+SOLAR, 3=SOLAR_IRRADIANCE, 4=SOLAR)
!       IPATH    PATH TYPE SWITCH.
!       ILAY     LAYER INDEX.
!       LMAX     LAYER LOOP MAXIMUM INDEX.
!       LSOLAR   SOLAR PATH SEGMENT INDEX.
!       MSPATH   LOGICAL, .TRUE. FOR MULTIPLE SCATTERING VERTICAL PATH.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       JT       BAND MODEL DATA TEMPERATURE INTERPOLATION UPPER INDEX.
!       FT       TEMPERATURE INTERPOLATION FRACTION USED WITH JT.
!       WPTH     INCREMENTAL MOLECULAR COLUMN AMOUNTS [ATM CM].
!       T5       SQUARE-ROOT OF SEGMENT TEMPERATURE TO 273.15K RATIO.
!       PTM75    PRESSURE OVER (T / 273.15K) RAISED TO THE 3/4 POWER.
!       PATM     SEGMENT DENSITY WEIGHTED PRESSURE [ATM].
!       P2       PRESSURE SQUARED INTERPOLATION FACTOR FOR SEGMENT.
!       JTS      SOLAR PATH TEMPERATURE INTERPOLATION UPPER INDEX.
!       FTS      TEMPERATURE INTERPOLATION FRACTION USED WITH JTS.
!       WSPTH    SCATTERING POINT TO SUN MOLECULAR COLUMN [ATM CM].
!       T5S      SQUARE-ROOT OF SOLAR PATH TEMPERATURE TO 273.15K RATIO.
!       PTM75S   PRESSURE OVER (T/273.15K) RAISED TO 3/4 FOR SOLAR PATH.
!       PSATM    SOLAR PATH DENSITY WEIGHTED PRESSURE [ATM].
!       P2S      PRESSURE SQUARED INTERPOLATION FACTOR FOR SOLAR PATH.
      LOGICAL MSPATH
      INTEGER IEMSCT,IPATH,ILAY,LMAX,LSOLAR,JT(*),KT(*),                &
     &  JTS(MMOLT,*),KTS(MMOLT,*)
      REAL FACMC,FT(*),GT(*),WPTH(0:MEXTXY,1:*),T5(*),PTM75(*),PATM(*), &
     &  P2(*),FTS(MMOLT,*),GTS(MMOLT,*),WSPTH(0:MEXTXY,1:*),            &
     &  T5S(MMOLT,*),PTM75S(MMOLT,*),PSATM(MMOLT,*),P2S(MMOLT,*)

!     COMMONS:
      INCLUDE 'BASE.h'
      INCLUDE 'BMDAT.h'

!     /ACTIVE/
!       NACT     NUMBER OF ACTIVE MOLECULES FOR CURRENT FREQUENCY BIN.
!       NACTBM   NUMBER OF ACTIVE BAND MODEL MOLECULES FOR CURRENT FREQ.
!       NACTX    NUMBER OF ACTIVE X & Y CROSS-SECTION MOLECULES.
!       MACTBM   LIST OF ACTIVE BAND MODEL MOLECULES FOR CURRENT FREQ.
!       MACTX    LIST OF ACTIVE X & Y CROSS-SECTION MOLECULES.
      INTEGER NACT,NACTBM,NACTX,MACTBM,MACTX
      COMMON/ACTIVE/NACT,NACTBM,NACTX,MACTBM(MMOLYT),MACTX(NMOLX)

!     /BMOD_0/
!       COL_O3   2ND ORDER PATH LORENTZ WIDTH SUM FOR OZONE [CM2 ATM2].
!       SDSUM    WEAK-LINE OPTICAL DEPTH SUM.
!       ODSUM    WEAK-LINE OPTICAL DEPTH WEIGHTED LINE-SPACING SUM [CM].
!       DOPSUM   PATH DOPPLER WIDTH SUM [CM].
!       COLSUM   PATH COLLISION (LORENTZ) WIDTH SUM [CM ATM].
      REAL COL_O3,SDSUM,ODSUM,DOPSUM,COLSUM,TAILSM
      COMMON/BMOD_0/COL_O3,SDSUM(MMOLT),ODSUM(MMOLT),DOPSUM(MMOLT),     &
     &  COLSUM(MMOLT),TAILSM(MMOLT,0:NTLSUB)
      SAVE /BMOD_0/

!     LOCAL FUNCTIONS:
      REAL BMTRAN,BMTRN

!     LOCAL VARIABLES:
!       IACT     LOOP INDEX FOR ACTIVE MOLECULES.
!       IMOL     MOLECULAR INDEX.
!       MOLPNT   COLUMN AMOUNT MOLECULAR INDEX.
!       LLAY     LAYER INDEX.
!       JUPTMP   TEMPERATURE INTERPOLATION UPPER INDEX.
!       JLOTMP   TEMPERATURE INTERPOLATION LOWER INDEX.
!       ITLSUB   LOOP INDEX FOR LINE TAIL SUB-INTERVALS.
!       COLDEN   MOLECULAR COLUMN DENSITY [ATM CM].
!       F        TEMPERATURE INTERPOLATION FRACTION.
!       ABSM     MOLECULAR ABSORPTION COEFFICIENT [CM-1/ATM].
!       DINV     LINE SPACING BAND MODEL PARAMETER [CM].
!       STORE    WORKING VARIABLE STORAGE.
!       TAIL     UPPER TEMPERATURE LINE TAIL ABS COEF [CM-1/ATM^2]
!       TAILM1   LOWER TEMPERATURE LINE TAIL ABS COEF [CM-1/ATM^2]
!       DEPTH    WEAK-LINE OPTICAL DEPTH CURTIS-GODSON SUM.
!       ODBAR    LINE SPACING CURTIS-GODSON SUM [CM].
!       ADBAR    DOPPLER WIDTH CURTIS-GODSON SUM [CM].
!       ACBAR    COLLISION (LORENTZ) WIDTH CURTIS-GODSON SUM [CM ATM].
!       ACBAR2   2ND ORDER COLLISION WIDTH CURTIS-GODSON SUM [CM2 ATM2].
!       TRANSM   MOLECULAR TRANSMITTANCE.
!       PRIOR    PRIOR FREQUENCY LINE-TAIL TRANSMITTANCE.
!       CURREN   CURRENT FREQUENCY LINE-TAIL TRANSMITTANCE.
!       DELTA    LINE-TAIL OPTICAL DEPTH DIFFERENCE.
!       TDEPTH   LINE-TAIL OPTICAL DEPTH.
!       LBMWID   FLAG, FALSE TO NUMERICALLY INTEGRATE EQUIVALENT WIDTH.
      INTEGER IACT,IMOL,MOLPNT,LLAY,JUPTMP,JLOTMP,ITLSUB
      REAL COLDEN,F,ABSM,DINV,STORE,TAIL,TAILM1,DEPTH,ODBAR,ADBAR,      &
     &  ACBAR,ACBAR2,TRANSM,PRIOR,CURREN,DELTA,TDEPTH(0:NTLSUB)
      LOGICAL LBMWID
      DO 10 IACT=1,NACTBM
          IMOL=ABS(MACTBM(IACT))
          MOLPNT=MPOINT(IMOL)

!         CHECK IF LINE CENTER CONTRIBUTES TO ABSORPTION
          IF(MACTBM(IACT).LT.0)THEN

!             LOOP OVER LAYERS
              DO LLAY=ILAY,LMAX

!                 PATH AMOUNT
                  COLDEN=FACMC*WPTH(MOLPNT,LLAY)
                  IF(COLDEN.GT.0.)THEN

!                     INTERPOLATE BAND MODEL PARAMETERS OVER TEMPERATURE
                      JUPTMP=JT(LLAY)
                      JLOTMP=JUPTMP-1
                      F=FT(LLAY)
                      ABSM=SDCN(JUPTMP,IMOL)
                      ABSM=ABSM+F*(SDCN(JLOTMP,IMOL)-ABSM)
                      DINV=ODCN(JUPTMP,IMOL)
                      DINV=DINV+F*(ODCN(JLOTMP,IMOL)-DINV)

!                     PERFORM CURTIS-GODSON SUMS
                      STORE=ABSM*COLDEN
                      SDSUM(IMOL)=SDSUM(IMOL)+STORE
                      STORE=DINV*STORE
                      ODSUM(IMOL)=ODSUM(IMOL)+STORE
                      DOPSUM(IMOL)=DOPSUM(IMOL)+STORE*T5(LLAY)
                      STORE=STORE*PTM75(LLAY)
                      COLSUM(IMOL)=COLSUM(IMOL)+STORE
                      IF(IMOL.EQ.3)COL_O3=COL_O3+STORE*DINV*PTM75(LLAY)
                      COLDEN=COLDEN*PATM(LLAY)

                      DO ITLSUB=0,NTLSUB
                          TAIL=SDTL(ITLSUB,JUPTMP,1,IMOL)
                          TAIL=TAIL+P2(LLAY)                            &
     &                      *(SDTL(ITLSUB,JUPTMP,2,IMOL)-TAIL)
                          TAILM1=SDTL(ITLSUB,JLOTMP,1,IMOL)
                          TAILM1=TAILM1+P2(LLAY)                        &
     &                      *(SDTL(ITLSUB,JLOTMP,2,IMOL)-TAILM1)
                          TAILSM(IMOL,ITLSUB)=TAILSM(IMOL,ITLSUB)       &
     &                      +COLDEN*(TAIL+F*(TAILM1-TAIL))
                      ENDDO
                  ENDIF
              ENDDO
              DEPTH=SDSUM(IMOL)
              ODBAR=ODSUM(IMOL)
              ADBAR=DOPSUM(IMOL)
              ACBAR=COLSUM(IMOL)
              ACBAR2=0.
              IF(IMOL.EQ.3)ACBAR2=COL_O3
              DO ITLSUB=0,NTLSUB
                  TDEPTH(ITLSUB)=TAILSM(IMOL,ITLSUB)
              ENDDO

!             IF SOLAR PATH, CALCULATE ADDITIONAL LAYER:
              IF(IEMSCT.EQ.2 .AND. IPATH.NE.3)THEN
                  COLDEN=FACMC*WSPTH(MOLPNT,LSOLAR)
                  IF(COLDEN.GT.0.)THEN
                      JUPTMP=JTS(IMOL,LSOLAR)
                      JLOTMP=JUPTMP-1
                      F=FTS(IMOL,LSOLAR)
                      ABSM=SDCN(JUPTMP,IMOL)
                      ABSM=ABSM+F*(SDCN(JLOTMP,IMOL)-ABSM)
                      DINV=ODCN(JUPTMP,IMOL)
                      DINV=DINV+F*(ODCN(JLOTMP,IMOL)-DINV)
                      IF(MSPATH)THEN

!                         SOLAR PATH ONLY
                          DEPTH=ABSM*COLDEN
                          ODBAR=DINV*DEPTH
                          ADBAR=ODBAR*T5S(IMOL,LSOLAR)
                          ACBAR=ODBAR*PTM75S(IMOL,LSOLAR)
                          IF(IMOL.EQ.3)                                 &
     &                      ACBAR2=ACBAR*DINV*PTM75S(IMOL,LSOLAR)
                          COLDEN=COLDEN*PSATM(IMOL,LSOLAR)
                          DO ITLSUB=0,NTLSUB
                              TAIL=SDTL(ITLSUB,JUPTMP,1,IMOL)
                              TAIL=TAIL+P2S(IMOL,LSOLAR)                &
     &                          *(SDTL(ITLSUB,JUPTMP,2,IMOL)-TAIL)
                              TAILM1=SDTL(ITLSUB,JLOTMP,1,IMOL)
                              TAILM1=TAILM1+P2S(IMOL,LSOLAR)            &
     &                          *(SDTL(ITLSUB,JLOTMP,2,IMOL)-TAILM1)
                              TDEPTH(ITLSUB)                            &
     &                          =COLDEN*(TAIL+F*(TAILM1-TAIL))
                          ENDDO
                      ELSE

!                         L-SHAPED PATH:
                          STORE=ABSM*COLDEN
                          DEPTH=DEPTH+STORE
                          STORE=DINV*STORE
                          ODBAR=ODBAR+STORE
                          ADBAR=ADBAR+STORE*T5S(IMOL,LSOLAR)
                          STORE=STORE*PTM75S(IMOL,LSOLAR)
                          ACBAR=ACBAR+STORE
                          IF(IMOL.EQ.3)                                 &
     &                      ACBAR2=ACBAR2+STORE*DINV*PTM75S(IMOL,LSOLAR)
                          COLDEN=COLDEN*PSATM(IMOL,LSOLAR)
                          DO ITLSUB=0,NTLSUB
                              TAIL=SDTL(ITLSUB,JUPTMP,1,IMOL)
                              TAIL=TAIL+P2S(IMOL,LSOLAR)                &
     &                          *(SDTL(ITLSUB,JUPTMP,2,IMOL)-TAIL)
                              TAILM1=SDTL(ITLSUB,JLOTMP,1,IMOL)
                              TAILM1=TAILM1+P2S(IMOL,LSOLAR)            &
     &                          *(SDTL(ITLSUB,JLOTMP,2,IMOL)-TAILM1)
                              TDEPTH(ITLSUB)=TDEPTH(ITLSUB)             &
     &                          +COLDEN*(TAIL+F*(TAILM1-TAIL))
                          ENDDO
                      ENDIF
                  ELSEIF(MSPATH)THEN
                      TX(MOLPNT)=1.
                      GOTO 10
                  ENDIF
              ENDIF

!             CHECK FOR WEAK LINE, NARROW HALF-WIDTH LIMIT:
              ADBAR=DOP0(IMOL)*ADBAR
              ACBAR=ALF0(IMOL)*ACBAR
              IF(DEPTH.LT..001 .AND. 10*(ADBAR+ACBAR).LE.ODBAR)THEN
!2BIN             IF(EDGENR.GT.0.)THEN
                      TRANSM=1-DEPTH
!2BIN             ELSE
!2BIN                 TRANSM=1-DEPTH/2
!2BIN             ENDIF
              ELSEIF(DEPTH.GT.0.)THEN

!                 CALCULATE EQUIVALENT WIDTH TRANSMITTANCE
                  ODBAR=ODBAR/DEPTH
                  ADBAR=ADBAR/DEPTH
                  ACBAR=ACBAR/DEPTH
                  LBMWID=.TRUE.
                  IF(ACBAR2.GT.0.)THEN
                      ACBAR2=ACBAR2/DEPTH*ALF0(IMOL)**2
                      TRANSM                                            &
     &                  =BMTRAN(DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2,LBMWID)
                  ELSE
                      TRANSM=BMTRN(DEPTH,ODBAR,ADBAR,ACBAR,LBMWID)
                  ENDIF
              ELSE
                  TRANSM=1.
              ENDIF

!             ADD LINE TAIL CONTRIBUTIONS.
              IF(TRANSM.GT.0.)THEN
                  TAIL=0.
                  PRIOR=EXP(-TDEPTH(0))
                  DO ITLSUB=1,NTLSUB
                      CURREN=EXP(-TDEPTH(ITLSUB))
                      DELTA=TDEPTH(ITLSUB)-TDEPTH(ITLSUB-1)
                      IF(ABS(DELTA).GE..01)THEN
                          TAIL=TAIL+(PRIOR-CURREN)/DELTA
                      ELSE
                          TAIL=TAIL+PRIOR*(1-DELTA*(1-DELTA/3)/2)
                      ENDIF
                      PRIOR=CURREN
                  ENDDO
                  TRANSM=TRANSM*TAIL/NTLSUB
              ENDIF
              TX(MOLPNT)=TRANSM

!         ONLY LINE TAILS CONTRIBUTE TO ABSORPTION
          ELSE

!             LOOP OVER LAYERS
              DO LLAY=ILAY,LMAX

!                 PERFORM CURTIS-GODSON SUMS
                  COLDEN=FACMC*WPTH(MOLPNT,LLAY)
                  IF(COLDEN.GT.0.)THEN

!                     LINE TAILS ARE SCALED BY PRESSURE;
!                     EXTINCTION CROSS-SECTIONS ARE NOT.
                      COLDEN=COLDEN*PATM(LLAY)

!                     INTERPOLATE BAND MODEL PARAMETERS OVER TEMPERATURE
                      JUPTMP=JT(LLAY)
                      JLOTMP=JUPTMP-1
                      F=FT(LLAY)
                      DO ITLSUB=0,NTLSUB
                          TAIL=SDTL(ITLSUB,JUPTMP,1,IMOL)
                          TAIL=TAIL+P2(LLAY)                            &
     &                      *(SDTL(ITLSUB,JUPTMP,2,IMOL)-TAIL)
                          TAILM1=SDTL(ITLSUB,JLOTMP,1,IMOL)
                          TAILM1=TAILM1+P2(LLAY)                        &
     &                      *(SDTL(ITLSUB,JLOTMP,2,IMOL)-TAILM1)
                          TAILSM(IMOL,ITLSUB)=TAILSM(IMOL,ITLSUB)       &
     &                      +COLDEN*(TAIL+F*(TAILM1-TAIL))
                      ENDDO
                  ENDIF
              ENDDO
              DO ITLSUB=0,NTLSUB
                  TDEPTH(ITLSUB)=TAILSM(IMOL,ITLSUB)
              ENDDO

!             IF SOLAR PATH, CALCULATE ADDITIONAL LAYER:
              IF(IEMSCT.EQ.2 .AND. IPATH.NE.3)THEN
                  COLDEN=FACMC*WSPTH(MOLPNT,LSOLAR)
                  IF(COLDEN.GT.0.)THEN

!                     LINE TAILS ARE SCALED BY PRESSURE;
!                     EXTINCTION CROSS-SECTIONS ARE NOT.
                      COLDEN=COLDEN*PSATM(IMOL,LSOLAR)
                      JUPTMP=JTS(IMOL,LSOLAR)
                      JLOTMP=JUPTMP-1
                      F=FTS(IMOL,LSOLAR)
                      IF(MSPATH)THEN

!                         SOLAR PATH ONLY:
                          DO ITLSUB=0,NTLSUB
                              TAIL=SDTL(ITLSUB,JUPTMP,1,IMOL)
                              TAIL=TAIL+P2S(IMOL,LSOLAR)*               &
     &                          (SDTL(ITLSUB,JUPTMP,2,IMOL)-TAIL)
                              TAILM1=SDTL(ITLSUB,JLOTMP,1,IMOL)
                              TAILM1=TAILM1+P2S(IMOL,LSOLAR)*           &
     &                          (SDTL(ITLSUB,JLOTMP,2,IMOL)-TAILM1)
                              TDEPTH(ITLSUB)                            &
     &                          =COLDEN*(TAIL+F*(TAILM1-TAIL))
                          ENDDO
                      ELSE

!                         L-SHAPED PATH:
                          DO ITLSUB=0,NTLSUB
                              TAIL=SDTL(ITLSUB,JUPTMP,1,IMOL)
                              TAIL=TAIL+P2S(IMOL,LSOLAR)*               &
     &                          (SDTL(ITLSUB,JUPTMP,2,IMOL)-TAIL)
                              TAILM1=SDTL(ITLSUB,JLOTMP,1,IMOL)
                              TAILM1=TAILM1+P2S(IMOL,LSOLAR)*           &
     &                          (SDTL(ITLSUB,JLOTMP,2,IMOL)-TAILM1)
                              TDEPTH(ITLSUB)=TDEPTH(ITLSUB)             &
     &                          +COLDEN*(TAIL+F*(TAILM1-TAIL))
                          ENDDO
                      ENDIF
                  ELSEIF(MSPATH)THEN
                      TX(MOLPNT)=1.
                      GOTO 10
                  ENDIF
              ENDIF

!             CALCULATE LINE TAIL TRANSMITTANCE
              PRIOR=EXP(-TDEPTH(0))
              TAIL=0.
              DO ITLSUB=1,NTLSUB
                  CURREN=EXP(-TDEPTH(ITLSUB))
                  DELTA=TDEPTH(ITLSUB)-TDEPTH(ITLSUB-1)
                  IF(ABS(DELTA).GE..01)THEN
                      TAIL=TAIL+(PRIOR-CURREN)/DELTA
                  ELSE
                      TAIL=TAIL+PRIOR*(1-DELTA*(1-DELTA/3)/2)
                  ENDIF
                  PRIOR=CURREN
              ENDDO
              TX(MOLPNT)=TAIL/NTLSUB
          ENDIF
   10 CONTINUE

!     START CROSS-SECTION MOLECULE LOOP:
      DO 20 IACT=1,NACTX
          IMOL=MACTX(IACT)
          MOLPNT=MPOINT(IMOL)

!         LOOP OVER LAYERS
          DO LLAY=ILAY,LMAX

!             PERFORM CURTIS-GODSON SUMS
              COLDEN=FACMC*WPTH(MOLPNT,LLAY)
              IF(COLDEN.GT.0.)THEN

!                 INTERPOLATE BAND MODEL PARAMETERS OVER TEMPERATURE
                  JUPTMP=KT(LLAY)
                  JLOTMP=JUPTMP-1
                  F=GT(LLAY)
                  TAIL=SDTL(0,JUPTMP,1,IMOL)
                  TAIL=TAIL+F*(SDTL(0,JLOTMP,1,IMOL)-TAIL)
                  TAILSM(IMOL,0)=TAILSM(IMOL,0)+COLDEN*TAIL
              ENDIF
          ENDDO
          TDEPTH(0)=TAILSM(IMOL,0)

!         IF SOLAR PATH, CALCULATE ADDITIONAL LAYER:
          IF(IEMSCT.EQ.2 .AND. IPATH.NE.3)THEN
              COLDEN=FACMC*WSPTH(MOLPNT,LSOLAR)
              IF(COLDEN.GT.0.)THEN

!                 LINE TAILS ARE SCALED BY PRESSURE;
!                 EXTINCTION CROSS-SECTIONS ARE NOT.
                  JUPTMP=KTS(IMOL,LSOLAR)
                  JLOTMP=JUPTMP-1
                  F=GTS(IMOL,LSOLAR)
                  TAIL=SDTL(0,JUPTMP,1,IMOL)
                  TAIL=TAIL+F*(SDTL(0,JLOTMP,1,IMOL)-TAIL)
                  IF(MSPATH)THEN

!                     SOLAR PATH ONLY:
                      TDEPTH(0)=COLDEN*TAIL
                  ELSE

!                     L-SHAPED PATH:
                      TDEPTH(0)=TDEPTH(0)+COLDEN*TAIL
                  ENDIF
              ELSEIF(MSPATH)THEN
                  TX(MOLPNT)=1.
                  GOTO 20
              ENDIF
          ENDIF

!         CALCULATE LINE TAIL TRANSMITTANCE
          TX(MOLPNT)=EXP(-TDEPTH(0))
   20 CONTINUE

!     RETURN TO BMOD:
      RETURN
      END
