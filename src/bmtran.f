      REAL FUNCTION BMTRAN(DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2,LBMWID)

!     THIS FUNCTION RETURNS IN-BAND PATH TRANSMITTANCE.
      IMPLICIT NONE

!     OUTPUTS:
!       BMTRAN  =   PATH TRANSMITTANCE [DIMENSIONLESS]

!     INPUTS:
!                    _
!       DEPTH   =   >_  SU/D,     THE ABSORPTION COEFFICIENT (S/D)
!                                 AND COLUMN DENSITY (U) PRODUCT,
!                                 SUMMED OVER THE PATH [DIMENSIONLESS]

!       ODBAR   =   <1/D>,        THE RECIPROCAL OF THE LINE
!                                 SPACING, (D), PATH AVERAGED [CM]

!       ADBAR   =   <ALFDOP/d>,   THE DOPPLER HALF-WIDTH (ALFDOP)
!                                 OVER THE AVERAGE LINE SPACING (D),
!                                 PATH AVERAGED [DIMENSIONLESS]

!       ACBAR   =   <ALFCOL/d>,   THE COLLISION (LORENTZ) HALF-WIDTH
!                                 (ALFCOL) OVER THE AVERAGE LINE SPACING
!                                 (D), PATH AVERAGED [DIMENSIONLESS]

!                          2  2
!       ACBAR2  =   <ALFCOL /d >, THE COLLISION (LORENTZ) HALF-WIDTH
!                                 (ALFCOL) SQUARED OVER THE AVERAGE
!                                 LINE SPACING (d) SQUARED, PATH
!                                 AVERAGED [DIMENSIONLESS] (OZONE ONLY)

!       LBMWID   FLAG, FALSE TO NUMERICALLY INTEGRATE EQUIVALENT WIDTH.

!     ARGUMENTS:
      REAL DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2
      LOGICAL LBMWID

!     PARAMETERS:
!       ONE6TH   ONE-SIXTH.
!       PT5LN2   LN(2)/2
!       P5OLN2   .5 / LN(2)
!       REC_PI   RECIPROCAL OF PI.
!       PI       3.14159...               2
!       PADED    PI (95 PI - 208) / (32 PI  - 124 PI + 96)
!       PADEB    4 D - 2 PI
!       PADEC    8 PI D - 3 B - 32 PI
!       PADEA    2 PI (C - B)
      REAL ONE6TH,PT5LN2,PI,TWOPI,PADEA,PADEB,PADEC,PADED
      DOUBLE PRECISION P5OLN2,REC_PI
      PARAMETER(PI=3.141592654,TWOPI=2*PI,REC_PI=.318309886183791D0,    &
     &  ONE6TH=1./6,PT5LN2=.34657359,P5OLN2=.721347520444482D0,         &
     &  PADED=PI*(95*PI-208)/((32*PI-124)*PI+96),PADEB=4*PADED-2*PI,    &
     &  PADEC=8*PI*PADED-3*PADEB-32*PI,PADEA=TWOPI*(PADEC-PADEB))
      INCLUDE 'PARAMS.h'

!     COMMONS:
!     /BMHEAD/
!       BNDWID   SPECTRAL BIN WIDTH [CM-1].
!       EDGENR   LINE CENTER TO SPECTRAL BIN NEARER EDGE
!                DISTANCE OVER SPECTRAL BIN WIDTH.
!       EDGEFR   LINE CENTER TO SPECTRAL BIN FURTHER EDGE
!                DISTANCE OVER SPECTRAL BIN WIDTH.
!       EDGEDF   EDGEFR MINUS EDGENR.
!                |                            |
!                |<-------------------------->|                 d
!                |<------>|   BNDWID          |      EDGENR = ------
!                |    d   |                   |               BNDWID
!                |      _/ \_                 |
!                |_____/     \________________|    EDGEFR = 1 - EDGENR
      INCLUDE 'BMHEAD.h'

!     LOCAL VARIABLES:
!       BLINES   EFFECTIVE NUMBER OF LINES IN THE SPECTRAL BIN.
!       HLFWDW   HALF THE WEAK-LINE EQUIVALENT WIDTH OVER THE BIN WIDTH.
!       SUBYGC   LINE STRENGTH S TIMES COLUMN DENSITY U OVER
!                LORENTZ (COLLISION) HALF WIDTH.
!       WDDENL   LORENTZ EQUIVALENT WIDTH PADE APPROXIMATE DENOMINATOR.
!       DOPTRM   TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH.
!       WDRATD   DOPPLER TO WEAK-LINE EQUIVALENT WIDTH RATIO SQUARED.
!       WDRATL   LORENTZ TO WEAK-LINE EQUIVALENT WIDTH RATIO SQUARED.
!       ONEML2   ONE MINUS THE LORENTZ TO WEAK-LINE EQUIVALENT
!                WIDTH RATIO SQUARED (1 - WDRATL).
!       VGTTRM   TERM IN EXPRESSION FOR VOIGT EQUIVALENT WIDTH.
!       WDRATD   DOPPLER TO WEAK-LINE EQUIVALENT WIDTH RATIO SQUARED.
!       HLFWDV   HALF TOTAL VOIGT EQUIVALENT WIDTH OVER BIN WIDTH.
!       GAMC2    SQUARE OF THE RATIO OF THE LORENTZ (COLLISION)
!                HALF-WIDTH TO LINE SPACING PARAMETER.
!       GAMD2    SQUARE OF THE RATIO OF THE DOPPLER HALF-WIDTH TO
!                LINE SPACING PARAMETER, ALL DIVIDED BY (2 LN2).
!       SUGCPI   PRODUCT OF [1] ABSORPTION COEFFICIENT (S/d, LINE
!                STRENGTH OVER LINE SPACING), [2] COLUMN AMOUNT U
!                AND [3] GAMC, I.E., THE RATIO OF LORENTZ HALF-WIDTH
!                TO LINE SPACING, ALL DIVIDED BY PI.
!       VGTRAT   VOIGT EQUIVALENT WIDTH IN A FINITE SPECTRAL INTERVAL
!                OVER THE INTERVAL WIDTH COMPUTED ANALYTICALLY.
      REAL SUBYGC,WDDENL,DOPTRM,WDRATD,WDRATL,ONEML2,VGTTRM,            &
     &  SRAT,RHO,RHOP1,RHO2,FAC,FRHO,STR2PI,SRAT2
      DOUBLE PRECISION BLINES,HLFWDW,HLFWDV,GAMC2,GAMD2,SUGCPI,VGTRAT

!     FUNCTIONS:
!       BMWID    FINITE BIN VOIGT EQUIVALENT WIDTH OVER FULL BIN WIDTH.
!       BMWID2   SUB-BIN VOIGT EQUIVALENT WIDTH OVER FULL BIN WIDTH.
!2BIN REAL BMWID,BMWID2
      DOUBLE PRECISION BMWID2

!     THE RODGERS AND WILLIAMS APPROXIMATION IS USED TO
!     CALCULATE THE TOTAL VOIGT EQUIVALENT WIDTH (WDV):

!                2            2            2            2          2
!       (WDV/WDW)  = (WDL/WDW)  + (WDD/WDW)  - (WDL/WDW)  (WDD/WDW)

!     WHERE WDW IS THE WEAK-LINE EQUIVALENT WIDTH,
!           WDL IS THE LORENTZ EQUIVALENT WIDTH, AND
!           WDD IS THE DOPPLER EQUIVALENT WIDTH.

!     THE APPROXIMATIONS FOR THE EQUIVALENT WIDTHS FOLLOW:
!                                       _
!          WDW   =   <SU>         =  [ >_ SU/D ] / <1/D>

!                             2        A + X (B + 4 X)
!       WDRATL   =   (WDL/WDW)   =  ---------------------
!                                   A + X [C + X (D + X)]

!                             2
!       WDRATD   =   (WDD/WDW)   =  LN (1 + DOPTRM) / DOPTRM

!     WHERE

!                    X  =  WDW/<ALFCOL>
!                                                   2
!                    D  =  PI (95 PI - 208) / (32 PI  - 124 PI + 96)

!                    B  =  4 D - 2 PI

!                    C  =  8 PI D - 3 B - 32 PI

!                    A  =  2 PI (C - B)

!                                                   2
!               DOPTRM  =  [LN(2)/2]  (WDW/<ALFDOP>)

!     THE HALF-WIDTHS ARE CALCULATED FROM THE FOLLOWING EXPRESSIONS:
!                                                                 _
!       <ALFCOL>/WDW = [<ALFCOL/D> / <1/D>] / WDW = <ALFCOL/D> / >_ SU/D
!                                                                 _
!       <ALFDOP>/WDW = [<ALFDOP/D> / <1/D>] / WDW = <ALFDOP/D> / >_ SU/D

!     FOR OZONE AN IMPROVED APPROXIMATION FOR LORENTZ EQUIVALENT WIDTH
!     IS REQUIRED BASED ON THE WORK OF R. M. GOODY, WITH REFINEMENTS
!     BY L. S. BERNSTEIN TO ELIMINATE DECREASES IN THE CURVE-OF-GROWTH

!     EFFECTIVE NUMBER OF LINES:
      BLINES=DBLE(BNDWID)*DBLE(ODBAR)

!     HALF THE WEAK-LINE EQUIVALENT WIDTH OVER THE BIN WIDTH:
      HLFWDW=DBLE(DEPTH)/(2*BLINES)

!     DENOMINATOR IN THE LORENTZ EQUIVALENT WIDTH EXPRESSION:
      SUBYGC=DEPTH/ACBAR
      SRAT=1.
      IF(ACBAR2.NE.0.)THEN

!         REFINEMENT FOR OZONE ONLY
          RHO=ACBAR2/ACBAR**2
          RHOP1=RHO+1
          RHO2=RHO**2
          FAC=RHO*(RHO2+1)/(RHOP1*(RHO2+RHOP1))

!         THE FRHO FIT HAS BEEN IMPROVED.  GOODY USED THE FORM:
!         FRHO=1.0088+(RHO-3)*(.4524+(RHO-3.25)*(.0028-(RHO-3.5)*.0064))
!         FRHO=FRHO/2
!         THE CURVE-OF-GROWTH DECREASES FOR SOME DOWN LOOKING
!         PATHS ON THE WING OF THE 9.6 MICRON O3 BAND.  WITH THE
!         NEW PARAMETERIZATION, THE PROBLEM IS LESS SEVERE:
          FRHO=.36*RHO-.2657
!         FRHO=RHO/2-.25
          STR2PI=SUBYGC/TWOPI
          SRAT=1-STR2PI*(RHO-1)*FAC                                     &
     &      /(2+4*STR2PI*(FRHO+STR2PI*RHO*FAC/(RHOP1**2)))
          SUBYGC=SRAT*SUBYGC
      ENDIF
      WDDENL=PADEA+SUBYGC*(PADEC+SUBYGC*(PADED+SUBYGC))
      SRAT2=SRAT**2

!     CALCULATE TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH.
      DOPTRM=PT5LN2*(DEPTH/ADBAR)**2

!     VOIGT TOTAL EQUIVALENT WIDTH (RODGERS-WILLIAMS APPROXIMATION):
      IF(DOPTRM.GT..01)THEN

!         STANDARD EXPRESSIONS
          WDRATD=LOG(1+DOPTRM)/DOPTRM
          WDRATL=SRAT2*(PADEA+SUBYGC*(PADEB+4*SUBYGC))/WDDENL
          HLFWDV=HLFWDW*SQRT(DBLE(WDRATD+WDRATL-WDRATD*WDRATL))
      ELSE

!         WDRATD IS NEAR ONE.  CALCULATE ONE MINUS WDRATL.
          ONEML2=(PADEA*(1-SRAT)*(1+SRAT)+SUBYGC*                       &
     &      (PADEC-SRAT2*PADEB+SUBYGC*(PADED-4*SRAT2+SUBYGC)))/WDDENL
          IF(DOPTRM.GT..0001)THEN

!             IF DOPTRM IS SMALL, WDRATD IS NEAR ONE.  REPLACE
!             THE LOG AND THE SQRT WITH POWER SERIES EXPANSIONS.
              VGTTRM=DOPTRM*(.25-DOPTRM*(ONE6TH-DOPTRM/8))*ONEML2
              HLFWDV=HLFWDW*DBLE(1-VGTTRM*(1+VGTTRM*(1+VGTTRM)/2))
          ELSE

!             IF DOPTRM IS VERY SMALL, TRUNCATED EXPANSION SUFFICES.
              HLFWDV=HLFWDW*DBLE(1-DOPTRM*ONEML2/4)
          ENDIF
      ENDIF

!     SUBTRACT LINE TAIL CONTRIBUTIONS ASSUMING HALF-WIDTHS
!     SMALL COMPARED TO LINE CENTER TO LINE TAIL DISTANCE:
      GAMC2=DBLE(ACBAR)**2
      GAMD2=P5OLN2*DBLE(ADBAR)**2
      SUGCPI=DBLE(DEPTH)*DBLE(ACBAR)*REC_PI
!2BIN IF(EDGENR.EQ.0.)THEN
!2BIN     VGTRAT=BMWID(HLFWDV,SUGCPI,ACBAR,GAMC2,GAMD2,BLINES)
!2BIN ELSE
      VGTRAT=BMWID2(DBLE(EDGEFR),HLFWDV,SUGCPI,                         &
     &              DBLE(ACBAR),GAMC2,GAMD2,BLINES,LBMWID)              &
     &      +BMWID2(DBLE(EDGENR),HLFWDV,SUGCPI,                         &
     &              DBLE(ACBAR),GAMC2,GAMD2,BLINES,LBMWID)
!2BIN ENDIF

!     CALCULATE TRANSMITTANCE USING POWER LAW EXPRESSION.
      IF(VGTRAT.GE.1.D0)THEN
          BMTRAN=0.
      ELSE
          BMTRAN=SNGL((1-VGTRAT)**BLINES)
      ENDIF

!     RETURN TO BMOD:
      RETURN
      END

      REAL FUNCTION BMTRN(DEPTH,ODBAR,ADBAR,ACBAR,LBMWID)

!     THIS FUNCTION RETURNS IN-BAND PATH TRANSMITTANCE.

!     OUTPUTS:
!       BMTRN   =   PATH TRANSMITTANCE [DIMENSIONLESS]

!     INPUTS:
!                    _
!       DEPTH   =   >_  SU/d,     THE ABSORPTION COEFFICIENT (S/d)
!                                 AND COLUMN DENSITY (U) PRODUCT,
!                                 SUMMED OVER THE PATH [DIMENSIONLESS]

!       ODBAR   =   <1/d>,        THE RECIPROCAL OF THE LINE
!                                 SPACING, (d), PATH AVERAGED [CM]

!       ADBAR   =   <ALFDOP/d>,   THE DOPPLER HALF-WIDTH (ALFDOP)
!                                 OVER THE AVERAGE LINE SPACING (D),
!                                 PATH AVERAGED [DIMENSIONLESS]

!       ACBAR   =   <ALFCOL/d>,   THE COLLISION (LORENTZ) HALF-WIDTH
!                                 (ALFCOL) OVER THE AVERAGE LINE SPACING
!                                 (d), PATH AVERAGED [DIMENSIONLESS]

!       LBMWID   FLAG, FALSE TO NUMERICALLY INTEGRATE EQUIVALENT WIDTH.

!     ARGUMENTS:
      REAL DEPTH,ODBAR,ADBAR,ACBAR
      LOGICAL LBMWID

!     PARAMETERS:
!       ONE6TH   ONE-SIXTH.
!       PT5LN2   LN(2)/2
!       P5OLN2   .5 / LN(2)
!       REC_PI   RECIPROCAL OF PI.
!       PI       3.14159...               2
!       PADED    PI (95 PI - 208) / (32 PI  - 124 PI + 96)
!       PADEB    4 D - 2 PI
!       PADEC    8 PI D - 3 B - 32 PI
!       PADEA    2 PI (C - B)
!       PADECB   PADEC - PADEB
!       PADED4   PADED - 4
      REAL ONE6TH,PT5LN2,PI,PADEA,PADEB,PADEC,PADED,PADECB,PADED4
      DOUBLE PRECISION P5OLN2,REC_PI
      PARAMETER(ONE6TH=1./6.,PT5LN2=.34657359,P5OLN2=.721347520444482D0,&
     &  PI=3.141592654,REC_PI=.318309886183791D0,                       &
     &  PADED=PI*(95*PI-208)/((32*PI-124)*PI+96),PADEB=4*PADED-2*PI,    &
     &  PADEC=8*PI*PADED-3*PADEB-32*PI,PADEA=2*PI*(PADEC-PADEB),        &
     &  PADECB=PADEC-PADEB,PADED4=PADED-4)
      INCLUDE 'PARAMS.h'

!     COMMONS:
!     /BMHEAD/
!       BNDWID   SPECTRAL BIN WIDTH [CM-1].
!       EDGENR   LINE CENTER TO SPECTRAL BIN NEARER EDGE
!                DISTANCE OVER SPECTRAL BIN WIDTH.
!       EDGEFR   LINE CENTER TO SPECTRAL BIN FURTHER EDGE
!                DISTANCE OVER SPECTRAL BIN WIDTH.
!       EDGEDF   EDGEFR MINUS EDGENR.
!                |                            |
!                |<-------------------------->|                 d
!                |<------>|   BNDWID          |      EDGENR = ------
!                |    d   |                   |               BNDWID
!                |      _/ \_                 |
!                |_____/     \________________|    EDGEFR = 1 - EDGENR
      INCLUDE 'BMHEAD.h'

!     LOCAL VARIABLES:
!       BLINES   EFFECTIVE NUMBER OF LINES IN THE SPECTRAL BIN.
!       BLINNR   EFFECTIVE NUMBER OF LINES FROM LINE CENTER
!                TO THE NEAREST SPECTRAL BIN EDGE.
!       HLFWDW   HALF THE WEAK-LINE EQUIVALENT WIDTH OVER THE BIN WIDTH.
!       SUBYGC   LINE STRENGTH S TIMES COLUMN DENSITY U OVER
!                LORENTZ (COLLISION) HALF WIDTH.
!       WDDENL   LORENTZ EQUIVALENT WIDTH PADE APPROXIMATE DENOMINATOR.
!       DOPTRM   TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH.
!       WDRATD   DOPPLER TO WEAK-LINE EQUIVALENT WIDTH RATIO SQUARED.
!       WDRATL   LORENTZ TO WEAK-LINE EQUIVALENT WIDTH RATIO SQUARED.
!       ONEML2   ONE MINUS THE LORENTZ TO WEAK-LINE EQUIVALENT
!                WIDTH RATIO SQUARED (1 - WDRATL).
!       VGTTRM   TERM IN EXPRESSION FOR VOIGT EQUIVALENT WIDTH.
!       WDRATD   DOPPLER TO WEAK-LINE EQUIVALENT WIDTH RATIO SQUARED.
!       HLFWDV   HALF TOTAL VOIGT EQUIVALENT WIDTH OVER BIN WIDTH.
!       GAMC2    SQUARE OF THE RATIO OF THE LORENTZ (COLLISION)
!                HALF-WIDTH TO LINE SPACING PARAMETER.
!       GAMD2    SQUARE OF THE RATIO OF THE DOPPLER HALF-WIDTH TO
!                LINE SPACING PARAMETER, ALL DIVIDED BY (2 LN2).
!       SUGCPI   PRODUCT OF [1] ABSORPTION COEFFICIENT (S/d, LINE
!                STRENGTH OVER LINE SPACING), [2] COLUMN AMOUNT U
!                AND [3] GAMC, I.E., THE RATIO OF LORENTZ HALF-WIDTH
!                TO LINE SPACING, ALL DIVIDED BY PI.
!       VGTRAT   VOIGT EQUIVALENT WIDTH IN A FINITE SPECTRAL INTERVAL
!                NORMALIZED BY THE INTERVAL WIDTH.
!       CVGTRT   ONE MINUS VGTRAT.
      REAL BLINNR,SUBYGC,WDDENL,DOPTRM,WDRATD,WDRATL,                   &
     &  ONEML2,VGTTRM,CVGTRT
      DOUBLE PRECISION BLINES,HLFWDW,HLFWDV,GAMC2,GAMD2,SUGCPI,VGTRAT

!     FUNCTIONS:
!       BMWID    FINITE BIN VOIGT EQUIVALENT WIDTH OVER FULL BIN WIDTH.
!       BMWID2   SUB-BIN VOIGT EQUIVALENT WIDTH OVER FULL BIN WIDTH.
!       CVGTWD   ONE MINUS THE QUOTIENT OF THE VOIGT EQUIVALENT WIDTH IN
!                A FINITE SPECTRAL INTERVAL OVER THE INTERVAL WIDTH.
!2BIN REAL BMWID,BMWID2,CVGTWD
      DOUBLE PRECISION BMWID2
      REAL CVGTWD

!     THE RODGERS AND WILLIAMS APPROXIMATION IS USED TO
!     CALCULATE THE TOTAL VOIGT EQUIVALENT WIDTH (WDV):

!                2            2            2            2          2
!       (WDV/WDW)  = (WDL/WDW)  + (WDD/WDW)  - (WDL/WDW)  (WDD/WDW)

!     WHERE WDW IS THE WEAK-LINE EQUIVALENT WIDTH,
!           WDL IS THE LORENTZ EQUIVALENT WIDTH, AND
!           WDD IS THE DOPPLER EQUIVALENT WIDTH.

!     THE APPROXIMATIONS FOR THE EQUIVALENT WIDTHS FOLLOW:
!                                       _
!          WDW   =   <SU>         =  [ >_ SU/D ] / <1/D>

!                             2        A + X (B + 4 X)
!       WDRATL   =   (WDL/WDW)   =  ---------------------
!                                   A + X [C + X (D + X)]

!                             2
!       WDRATD   =   (WDD/WDW)   =  LN (1 + DOPTRM) / DOPTRM

!     WHERE

!                    X  =  WDW/<ALFCOL>
!                                                   2
!                    D  =  PI (95 PI - 208) / (32 PI  - 124 PI + 96)

!                    B  =  4 D - 2 PI

!                    C  =  8 PI D - 3 B - 32 PI

!                    A  =  2 PI (C - B)

!                                                   2
!               DOPTRM  =  [LN(2)/2]  (WDW/<ALFDOP>)

!     THE HALF-WIDTHS ARE CALCULATED FROM THE FOLLOWING EXPRESSIONS:
!                                                                 _
!       <ALFCOL>/WDW = [<ALFCOL/D> / <1/D>] / WDW = <ALFCOL/D> / >_ SU/D
!                                                                 _
!       <ALFDOP>/WDW = [<ALFDOP/D> / <1/D>] / WDW = <ALFDOP/D> / >_ SU/D

!     EFFECTIVE NUMBER OF LINES:
      BLINES=DBLE(BNDWID)*DBLE(ODBAR)
      IF(LBMWID)THEN

!         HALF THE WEAK-LINE EQUIVALENT WIDTH OVER THE BIN WIDTH:
          HLFWDW=DBLE(DEPTH)/(2*BLINES)

!         DENOMINATOR IN THE LORENTZ EQUIVALENT WIDTH EXPRESSION:
          SUBYGC=DEPTH/ACBAR
          WDDENL=PADEA+SUBYGC*(PADEC+SUBYGC*(PADED+SUBYGC))

!         CALCULATE TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH.
          DOPTRM=PT5LN2*(DEPTH/ADBAR)**2

!         VOIGT TOTAL EQUIVALENT WIDTH (RODGERS-WILLIAMS APPROXIMATION):
          IF(DOPTRM.GT..01)THEN

!             STANDARD EXPRESSIONS
              WDRATD=LOG(1+DOPTRM)/DOPTRM
              WDRATL=(PADEA+SUBYGC*(PADEB+4*SUBYGC))/WDDENL
              HLFWDV=HLFWDW*SQRT(DBLE(WDRATD+WDRATL-WDRATD*WDRATL))
          ELSE

!             WDRATD IS NEAR ONE.  CALCULATE ONE MINUS WDRATL.
              ONEML2=SUBYGC*(PADECB+SUBYGC*(PADED4+SUBYGC))/WDDENL
              IF(DOPTRM.GT..0001)THEN

!                 IF DOPTRM IS SMALL, WDRATD IS NEAR ONE.  REPLACE
!                 THE LOG AND THE SQRT WITH POWER SERIES EXPANSIONS.
                  VGTTRM=DOPTRM*(.25-DOPTRM*(ONE6TH-DOPTRM/8))*ONEML2
                  HLFWDV=HLFWDW*DBLE(1-VGTTRM*(1+VGTTRM*(1+VGTTRM)/2))
              ELSE

!                 IF DOPTRM IS VERY SMALL, TRUNCATED EXPANSION SUFFICES.
                  HLFWDV=HLFWDW*DBLE(1-DOPTRM*ONEML2/4)
              ENDIF
          ENDIF

!         SUBTRACT LINE TAIL CONTRIBUTIONS ASSUMING HALF-WIDTHS
!         SMALL COMPARED TO LINE CENTER TO LINE TAIL DISTANCE:
          GAMC2=DBLE(ACBAR)**2
          GAMD2=P5OLN2*DBLE(ADBAR)**2
          SUGCPI=DBLE(DEPTH)*DBLE(ACBAR)*REC_PI
          VGTRAT=BMWID2(DBLE(EDGEFR),HLFWDV,SUGCPI,                     &
     &                  DBLE(ACBAR),GAMC2,GAMD2,BLINES,LBMWID)          &
     &          +BMWID2(DBLE(EDGENR),HLFWDV,SUGCPI,                     &
     &                  DBLE(ACBAR),GAMC2,GAMD2,BLINES,LBMWID)

!         CHECK IF LBMWID IS STILL TRUE AND VGTRAT IS NOT TOO LARGE:
          IF(VGTRAT.LT..9D0 .AND. LBMWID)THEN

!             CALCULATE TRANSMITTANCE USING POWER LAW EXPRESSION:
              BMTRN=SNGL((1-VGTRAT)**BLINES)
              RETURN
          ENDIF
          LBMWID=.FALSE.
      ENDIF

!     SINGLE LINE TRANSMITTANCE IS LESS THAN 10%.  INCREASE ACCURACY
!     BY NUMERICALLY INTEGRATING THE VOIGT LINE SHAPE:
      BLINNR=EDGENR*SNGL(BLINES)
      CVGTRT=2*EDGENR*CVGTWD(DEPTH,ACBAR,ADBAR,BLINNR,0.)               &
     &  +EDGEDF*CVGTWD(DEPTH,ACBAR,ADBAR,EDGEDF*SNGL(BLINES),BLINNR)
      IF(CVGTRT.LE.0.)THEN
          BMTRN=0.
      ELSE
          BMTRN=CVGTRT**SNGL(BLINES)
      ENDIF

!     RETURN TO BMOD:
      RETURN
      END

!SAV  REAL FUNCTION BMWID(HLFWDV,SUGCPI,GAMC,GAMC2,GAMD2,BLINES)

!     BMWID RETURNS THE FINITE BIN VOIGT EQUIVALENT WIDTH
!     DIVIDED BY THE WIDTH OF THE FINITE BIN.

!                    1
!                    /    /                       2   \
!          BMWID =   |    |  1 - EXP[ -T(R BLINES)  ] |  d R
!                    /    \                           /
!                    0

!     WHERE
!                                                               2
!                                     1/2  INF        - (LN 2) Z
!                  2          / LN 2 \     /         E            DZ
!              T(N)  = SUGCPI | ---- |     |     -----------------------
!                             \  PI  /     /         2                 2
!                                         -INF   GAMC  + ( N - GAMD Z )

!                        2
!              AND   GAMD  = GAMD2

!     THE INPUT HLFWDV IS HALF THE TOTAL VOIGT EQUIVALENT WIDTH OVER
!     THE SPECTRAL BIN WIDTH.  UNLESS T(BLINES) IS LARGE (>16),
!     BMWID IS CALCULATED AS THE DIFFERENCE BETWEEN HLFWDV AND THE
!     TAIL CONTRIBUTION (BMTAIL).  THE BMTAIL CALCULATION ASSUMES
!     THAT THE LORENTZ AND DOPPLER HALF-WIDTHS ARE SMALL COMPARED TO
!     THE LINE TAIL OFFSET, I.E., GAMC << BLINES AND GAMD << BLINES.

!                          /                                       \
!              2   SUGCPI  |      GAMD2 /     15 GAMD2 - 4 GAMC \  |
!          T(N)  = ------  |  1 + ----- | 3 + ----------------- |  |
!                  DENOM   |      DENOM \           DENOM       /  |
!                          \                                       /

!                                2
!          DENOM = GAMC2 + BLINES

!                     /                     2   \
!          BMTAIL = - | 1 - EXP [ -T(BLINES)  ] |
!                     \                         /

!                          T(BLINES)        /
!            SQRT(SUGCPI)      /         2  |     / 3 GAMD2 - GAMC2 \  2
!          + ------------      |  EXP (-t ) | 2 + | --------------- | t
!               BLINES         /            |     \     SUGCPI      /
!                              0            \

!                             2                         2      \
!                   / 15 GAMD2  - 10 GAMD2 GAMC2 - GAMC2  \  4 |
!                 + | ----------------------------------- | t  |  d t
!                   \                      2              /    |
!                                  4 SUGCPI                    /

!     THE ERROR FUNCTION INTEGRALS ARE EVALUATED DIFFERENTLY
!     DEPENDING ON THE MAGNITUDE OF T(BLINES)**2 [VARIABLE "OPTDEP"].

!     INPUT ARGUMENTS (NOTE THAT THE ARGUMENTS CAN ALL BE SCALED BY
!     AN ARBITRARY POSITIVE CONSTANT WITHOUT CHANGING THE FUNCTION).
!       HLFWDV   HALF TOTAL VOIGT EQUIVALENT WIDTH OVER BIN WIDTH.
!       SUGCPI   PRODUCT OF [1] ABSORPTION COEFFICIENT (S/d, LINE
!                STRENGTH OVER LINE SPACING), [2] COLUMN AMOUNT U
!                AND [3] GAMC, I.E., THE RATIO OF LORENTZ HALF-WIDTH
!                TO LINE SPACING, ALL DIVIDED BY PI.
!       GAMC     LORENTZ (COLLISION) HALF-WIDTH TO LINE SPACING RATIO.
!       GAMC2    SQUARE OF GAMC.
!       GAMD2    SQUARE OF THE RATIO OF THE DOPPLER HALF-WIDTH TO
!                LINE SPACING PARAMETER, ALL DIVIDED BY (2 LN2).
!       BLINES   NUMBER OF LINES IN FULL SPECTRAL BIN.
!SAV  REAL HLFWDV,SUGCPI,GAMC,GAMC2,GAMD2,BLINES

!     PARAMETERS:
!       PCON     CONSTANT USED IN RATIONAL APPROXIMATION TO THE ERROR
!                FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A1       LINEAR COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A2       QUADRATIC COEFFICIENT USED IN RATIONAL APPROXIMATION
!                TO THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A3       CUBIC COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A4       QUARTIC COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A5       5TH-ORDER COEFFICIENT USED IN RATIONAL APPROXIMATION
!                TO THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       RTPI     SQUARE ROOT OF PI.
!       RTPI4    SQUARE ROOT OF PI DIVIDED BY 4.
!       RTPIC    SQUARE ROOT OF PI TIMES 3/32.
!SAV  REAL PCON,A1,A2,A3,A4,A5,RTPI,RTPI4,RTPIC
!SAV  PARAMETER(PCON=.3275911,A1=.451673692,A2=-.504257335,             &
!SAV &  A3=2.519390259,A4=-2.575644906,A5=1.881292140,                  &
!SAV &  RTPI=A1+A2+A3+A4+A5,RTPI4=RTPI/4,RTPIC=3*RTPI/32)

!     LOCAL VARIABLES:
!       COEF2    COEFFICIENT OF VOIGT LINE TAIL SERIES QUADRATIC TERM.
!       COEF4    COEFFICIENT OF VOIGT LINE TAIL SERIES QUARTIC TERM.
!       BLINSQ   BLINES SQUARED.
!       DENLOR   LORENTZ LINE SHAPE FUNCTION DENOMINATOR EVALUATED AT
!                THE BIN EDGE AND DIVIDED BY THE LINE SPACING SQUARED.
!       STOR0    CONSTANT EQUAL TO 5 TIMES COEF2 PLUS GAMC2.
!       OPTDEP   OPTICAL DEPTH EVALUATED AT THE SPECTRAL BIN EDGE.
!       TAIL     5 LORENTZ LINE TAIL INTEGRALS.
!       STOR1    CONSTANT EQUAL TO 9 TIMES GAMD2 PLUS STOR0.
!       BMTAIL   TAIL CONTRIBUTION TO NORMALIZED VOIGT EQUIVALENT WIDTH.
!       C4RAT    COEF4 DIVIDED BY SUGCPI.
!       RTDEP    SQUARE ROOT OF THE OPTICAL DEPTH.
!       PHI      PRODUCT OF RTPI, THE EXPONENTIAL FUNCTION EVALUATED AT
!                OPTDEP, AND THE COMPLEMENTARY ERROR FUNCTION OF RTDEP.
!       RTRAT    RATIO OF THE SQUARE ROOT OF SUGCPI TO BLINES.
!       RAPP     RATIO USED IN RATIONAL APPROXIMATION TO ERROR FUNCTION.
!SAV  REAL COEF2,COEF4,BLINSQ,DENLOR,STOR0,OPTDEP,                      &
!SAV &  TAIL(5),STOR1,BMTAIL,C4RAT,RTDEP,PHI,RTRAT,RAPP

!     DEFINE REQUIRED LOCAL VARIABLES:
!SAV  COEF2=3*GAMD2-GAMC2
!SAV  COEF4=5*(COEF2-GAMC2)*GAMD2-GAMC2**2
!SAV  BLINSQ=BLINES**2
!SAV  DENLOR=GAMC2+BLINSQ
!SAV  STOR0=5*COEF2+GAMC2
!SAV  OPTDEP=SUGCPI*(1+GAMD2*(3.+STOR0/DENLOR)/DENLOR)/DENLOR

!     BRANCH BASED ON MAGNITUDE OF THE OPTICAL DEPTH:
!SAV  IF(OPTDEP.LE..01)THEN

!         SMALL OPTICAL DEPTH:
!SAV      CALL LOREN5(BLINES,BLINSQ,GAMC,GAMC2,TAIL)
!SAV      STOR1=9*GAMD2+STOR0
!SAV      BMTAIL=SUGCPI*(TAIL(1)+GAMD2*(3*TAIL(2)+STOR0*TAIL(3))        &
!SAV &      -(SUGCPI/2)*(TAIL(2)+GAMD2*(6*TAIL(3)+(STOR0+STOR1)*TAIL(4))&
!SAV &      -SUGCPI*(TAIL(3)/3+GAMD2*(3*TAIL(4)+STOR1*TAIL(5)))))/BLINES
!SAV  ELSEIF(OPTDEP.GT.16.)THEN

!         SPECTRAL BIN IS BLACK:
!SAV      BMWID=1.
!SAV      RETURN
!SAV  ELSE

!         LARGE OPTICAL DEPTH CONTRIBUTION:
!SAV      C4RAT=COEF4/SUGCPI
!SAV      RTRAT=SQRT(SUGCPI)/BLINES
!SAV      BMTAIL=RTRAT*(RTPI+(RTPI4*COEF2+RTPIC*C4RAT)/SUGCPI)-1
!SAV      IF(OPTDEP.LT.11.)THEN

!             INTERMEDIATE OPTICAL DEPTH:
!SAV          RTDEP=SQRT(OPTDEP)
!SAV          IF(RTDEP.LT.2.33)THEN

!                 RATIONAL APPROXIMATION TO ERROR FUNCTION:
!SAV              RAPP=1/(1+PCON*RTDEP)
!SAV              PHI=RAPP*(A1+RAPP*(A2+RAPP*(A3+RAPP*(A4+RAPP*A5))))
!SAV          ELSE

!                 CONTINUED FRACTION REPRESENTATION FOR ERROR FUNCTION:
!SAV              PHI=(2.+OPTDEP*(4.5+OPTDEP))                          &
!SAV &              /(RTDEP*(3.75+OPTDEP*(5+OPTDEP)))
!SAV          ENDIF
!SAV          BMTAIL=BMTAIL+EXP(-OPTDEP)                                &
!SAV &          *(1-RTRAT*(PHI+(2*COEF2*(2*RTDEP+PHI)                  &
!SAV &          +C4RAT*(RTDEP*(OPTDEP+1.5)+.75*PHI))/(8*SUGCPI)))
!SAV      ENDIF
!SAV  ENDIF
!SAV  BMWID=HLFWDV-BMTAIL
!SAV  RETURN
!SAV  END

      REAL FUNCTION CVGTWD(DEPTH,ACBAR,ADBAR,BLINES,FMN)

!     CVGTWD RETURNS ONE MINUS THE QUOTIENT OF THE VOIGT EQUIVALENT
!     WIDTH IN A FINITE SPECTRAL INTERVAL OVER THE INTERVAL WIDTH.

!     ARGUMENTS:
!       DEPTH    PRODUCT OF LINE STRENGTH AND COLUMN AMOUNT OVER
!                LINE SPACING.
!       ACBAR    LORENTZ HALF-WIDTH OVER THE LINE SPACING.
!       ADBAR    HALF-WIDTH OVER THE LINE SPACING.
!       BLINES   EFFECTIVE NUMBER OF LINES IN SPECTRAL BIN (EQUAL
!                TO THE SPECTRAL BIN WIDTH OVER THE LINE SPACING).
!       FMN      SPECTRAL BIN MINIMUM FREQUENCY OVER THE LINE SPACING.
        REAL DEPTH,ACBAR,ADBAR,BLINES,FMN

!     LOCAL VARIABLES:
!       NPOINT   ONE LESS THAN THE NUMBER OF INTEGRATION POINTS.
!       IPOINT   INDEX FOR INTEGRATION POINTS.
!       COEF     COEFFICIENTS USED IN EQUIVALENT WIDTH CALCULATIONS:
!       X        REAL PART OF THE COMPLEX PROBABILITY FUNCTION
!                ARGUMENT [EQUAL TO SQRT(LN 2) TIMES THE FREQUENCY
!                DISPLACEMENT OVER THE DOPPLER HALF-WIDTH].
!       DELX     THE INTERGRATION STEP SIZE.
!       Y1P5     IMAGINARY PART OF THE COMPLEX PROBABILITY FUNCTION
!                ARGUMENT [EQUAL TO SQRT(LN 2) TIMES THE LORENTZ
!                HALF-WIDTH OVER THE DOPPLER HALF-WIDTH] PLUS 1.5.
!       Y1P5SQ   Y1P5 SQUARED.
      INTEGER NPOINT,IPOINT
      REAL COEF,X,DELX,Y1P5,Y1P5SQ

!     FUNCTIONS:
!       RCPF12   REAL PART OF THE COMPLEX PROBABILITY FUNCTION.
      REAL RCPF12

!     DATA:
!       RTPI     SQUARE ROOT OF PI.
!       RTLN2    SQUARE ROOT OF LN(2).
      REAL RTPI,RTLN2
      DATA RTPI,RTLN2/1.7724539,.83255461/

!     INITIALIZATIONS:
      NPOINT=4*(INT(BLINES/MAX(ACBAR,ADBAR))+4)
      COEF=RTLN2/ADBAR
      X=COEF*FMN
      DELX=COEF*BLINES/NPOINT
      Y1P5=COEF*ACBAR+1.5
      Y1P5SQ=Y1P5**2
      COEF=-DEPTH*COEF/RTPI
      CVGTWD=EXP(COEF*RCPF12(X,Y1P5,Y1P5SQ))/2

!     SET Y1P5 TO ZERO TO INDICATE THAT ITS VALUE
!     WAS PREVIOUSLY PASSED TO FUNCTION RCPF12.
      Y1P5=0.

!     LOOP OVER INTERMEDIATE INTEGRATION POINTS:
      DO IPOINT=2,NPOINT
          X=X+DELX
          CVGTWD=CVGTWD+EXP(COEF*RCPF12(X,Y1P5,Y1P5SQ))
      ENDDO

!     ADD FINAL INTEGRATION POINT
      X=X+DELX
      CVGTWD=CVGTWD+EXP(COEF*RCPF12(X,Y1P5,Y1P5SQ))/2
      CVGTWD=CVGTWD/NPOINT
      RETURN
      END

      SUBROUTINE LOREN5(ELINES,ELINSQ,GAMC,GAMC2,TAIL)

!     LOREN5 COMPUTES INTEGRAL OVER THE LORENTZ LINE TAIL
!     RAISED TO POWERS ONE THROUGH FIVE:

!                     INF
!                      /           dV
!         TAIL(M)  =   |    -----------------   ;  M = 1, 2, 3, 4 AND 5.
!                      /           2     2  M
!                   ELINES   ( GAMC  +  V  )

!     INPUT ARGUMENTS:
!       ELINES   NUMBER OF LINES BETWEEN LINE CENTER AND BIN EDGE.
!       ELINSQ   ELINES SQUARED.
!       GAMC     LORENTZ (COLLISION) HALF-WIDTH TO LINE SPACING RATIO.
!       GAMC2    SQUARE OF GAMC.
      DOUBLE PRECISION ELINES,ELINSQ,GAMC,GAMC2

!     OUTPUT ARGUMENTS:
!       TAIL     THE M LINE TAIL INTEGRALS [CM**(2M-1)].
      DOUBLE PRECISION TAIL(5)

!     PARAMETERS:
!       MX       THE POWER VALUE.
!       XMN      THE N-TH SERIES EXPANSION COEFFICIENT FOR TAIL(M).
      DOUBLE PRECISION M1,X10,X11,X12,X13,X14,M2,X20,X21,X22,X23,X24,   &
     &                 M3,X30,X31,X32,X33,X34,M4,X40,X41,X42,X43,X44,   &
     &                 M5,X50,X51,X52,X53,X54
      PARAMETER(M1=1.D0,X10=1/(2*M1-1),X11=M1/(2*M1+1),                 &
     &  X12=M1*(M1+1)/(4*M1+6),X13=M1*(M1+1)*(M1+2)/(12*M1+30),         &
     &  X14=M1*(M1+1)*(M1+2)*(M1+3)/(48*M1+168),                        &
     &  M2=2.D0,X20=1/(2*M2-1),X21=M2/(2*M2+1),                         &
     &  X22=M2*(M2+1)/(4*M2+6),X23=M2*(M2+1)*(M2+2)/(12*M2+30),         &
     &  X24=M2*(M2+1)*(M2+2)*(M2+3)/(48*M2+168),                        &
     &  M3=3.D0,X30=1/(2*M3-1),X31=M3/(2*M3+1),                         &
     &  X32=M3*(M3+1)/(4*M3+6),X33=M3*(M3+1)*(M3+2)/(12*M3+30),         &
     &  X34=M3*(M3+1)*(M3+2)*(M3+3)/(48*M3+168),                        &
     &  M4=4.D0,X40=1/(2*M4-1),X41=M4/(2*M4+1),                         &
     &  X42=M4*(M4+1)/(4*M4+6),X43=M4*(M4+1)*(M4+2)/(12*M4+30),         &
     &  X44=M4*(M4+1)*(M4+2)*(M4+3)/(48*M4+168),                        &
     &  M5=5.D0,X50=1/(2*M5-1),X51=M5/(2*M5+1),                         &
     &  X52=M5*(M5+1)/(4*M5+6),X53=M5*(M5+1)*(M5+2)/(12*M5+30),         &
     &  X54=M5*(M5+1)*(M5+2)*(M5+3)/(48*M5+168))

!     LOCAL VARIABLES:
!       DENLOR   LORENTZ LINE SHAPE FUNCTION DENOMINATOR EVALUATED AT
!                THE BIN EDGE AND DIVIDED BY THE LINE SPACING SQUARED.
!       TWOGM2   TWICE GAMC2.
!       DENOM    STORED DENOMINATOR].
!       RAT      STORED RATIO].
      DOUBLE PRECISION DENLOR,TWOGM2,DENOM,RAT

!     BRANCH BASED ON RELATIVE VALUE OF ELINES AND GAMC:
      IF(ELINSQ.LT.9*GAMC2)THEN

!         RECURSION FORMULA:
          TAIL(1)=ATAN2(GAMC,ELINES)/GAMC
          DENLOR=GAMC2+ELINSQ
          TWOGM2=GAMC2+GAMC2
          DENOM=TWOGM2
          RAT=ELINES/DENLOR
          TAIL(2)=(TAIL(1)-RAT)/DENOM
          DENOM=DENOM+TWOGM2
          RAT=RAT/DENLOR
          TAIL(3)=(3*TAIL(2)-RAT)/DENOM
          DENOM=DENOM+TWOGM2
          RAT=RAT/DENLOR
          TAIL(4)=(5*TAIL(3)-RAT)/DENOM
          DENOM=DENOM+TWOGM2
          RAT=RAT/DENLOR
          TAIL(5)=(7*TAIL(4)-RAT)/DENOM
      ELSE

!         POWER SERIES EXPANSION:
          RAT=GAMC2/ELINSQ
          DENOM=ELINES
          IF(ELINSQ.GE.81*GAMC2)THEN

!             3 TERMS:
              TAIL(1)=(X10-RAT*(X11-RAT*X12))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(2)=(X20-RAT*(X21-RAT*X22))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(3)=(X30-RAT*(X31-RAT*X32))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(4)=(X40-RAT*(X41-RAT*X42))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(5)=(X50-RAT*(X51-RAT*X52))/DENOM
          ELSE

!             5 TERMS:
              TAIL(1)=(X10-RAT*(X11-RAT*(X12-RAT*(X13-RAT*X14))))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(2)=(X20-RAT*(X21-RAT*(X22-RAT*(X23-RAT*X24))))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(3)=(X30-RAT*(X31-RAT*(X32-RAT*(X33-RAT*X34))))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(4)=(X40-RAT*(X41-RAT*(X42-RAT*(X43-RAT*X44))))/DENOM
              DENOM=DENOM*ELINSQ
              TAIL(5)=(X50-RAT*(X51-RAT*(X52-RAT*(X53-RAT*X54))))/DENOM
          ENDIF
      ENDIF
      RETURN
      END

      REAL FUNCTION RCPF12(X,Y1P5,Y1P5SQ)

!     RCPF12 RETURNS THE REAL PART OF THE COMPLEX PROBABILITY
!     FUNCTION [ W(Z)=EXP(-Z**2)*ERFC(-I*Z);   Z = X + IY ]
!     IN THE UPPER HALF-PLANE (I.E., FOR Y >= 0.).  THE MAXIMUM
!     RELATIVE ERROR OF RCPF12 IS <2.E-6.  THIS ROUTINE WAS
!     DEVELOPED BY J.HUMLICEK, JQSRT, VOL 21, P309 (1980)

!     ARGUMENTS:
!       X        REAL PART OF ARGUMENT.
!       Y1P5     IMAGINARY PART OF ARGUMENT PLUS 1.5.
!       Y1P5SQ   Y1P5 SQUARED.
      REAL X,Y1P5,Y1P5SQ

!     LOCAL VARIABLES:
!       XMH      X MINUS THE ROOT OF A HERMITE POLYNOMIAL.
!       XPH      X PLUS THE ROOT OF A HERMITE POLYNOMIAL.
!       CY1P5    COTANGENT DATA TIMES Y1P5.
      REAL XMH,XPH,CY1P5(6)

!     DATA:
!       HERM12   ZEROES OF THE 12TH ORDER HERMITE POLYNOMIAL
!                (ABRAMOWITZ AND STEGUN, TABLE 25.10).

!                  WHRM12   9/4
!       SINDAT  =  ------  E     SIN( 3 HERM12 )  ,  WHERE WHRM12 IS
!                    PI

!                THE WEIGHT OF THE 12TH ORDER HERMITE POLYNOMIAL.
!       CTNDAT   COTANGENT DATA (DERIVED FROM COSDAT AND SINDAT).
      REAL HERM12(6),CTNDAT(6),SINDAT(6)
      SAVE HERM12,CTNDAT,SINDAT,CY1P5
      DATA HERM12/ 0.314240376,     0.947788391,     1.59768264,        &
     &             2.27950708,      3.02063703,      3.8897249     /,   &
     &     CTNDAT/ 0.72617082,     -3.25314144,     -0.08083430,        &
     &             1.61167866,     -2.63380063,     -0.798131224   /,   &
     &     SINDAT/ 1.393237,        0.231152406,    -0.155351466,       &
     &             6.21836624E-03,  9.19082986E-05, -6.27525958E-07/
!OLD
!OLD  DATA COSDAT/ 1.01172805,     -0.75197147,      1.2557727E-02,
!OLD 1             1.00220082E-02, -2.42068135E-04,  5.00848061E-07/,

!     Y1P5 IS INPUT AS ZERO TO INDICATE THAT ARRAY CY1P5
!     HAS ALREADY BEEN DEFINED:
      IF(Y1P5.GT.0.)THEN

!         LOOP OVER DATA:
          CY1P5(1)=CTNDAT(1)*Y1P5
          XMH=X-HERM12(1)
          XPH=X+HERM12(1)
          RCPF12=SINDAT(1)*((CY1P5(1)-XMH)/(XMH**2+Y1P5SQ)              &
     &                     +(CY1P5(1)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(2)=CTNDAT(2)*Y1P5
          XMH=X-HERM12(2)
          XPH=X+HERM12(2)
          RCPF12=RCPF12+SINDAT(2)*((CY1P5(2)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(2)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(3)=CTNDAT(3)*Y1P5
          XMH=X-HERM12(3)
          XPH=X+HERM12(3)
          RCPF12=RCPF12+SINDAT(3)*((CY1P5(3)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(3)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(4)=CTNDAT(4)*Y1P5
          XMH=X-HERM12(4)
          XPH=X+HERM12(4)
          RCPF12=RCPF12+SINDAT(4)*((CY1P5(4)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(4)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(5)=CTNDAT(5)*Y1P5
          XMH=X-HERM12(5)
          XPH=X+HERM12(5)
          RCPF12=RCPF12+SINDAT(5)*((CY1P5(5)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(5)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(6)=CTNDAT(6)*Y1P5
          XMH=X-HERM12(6)
          XPH=X+HERM12(6)
          RCPF12=RCPF12+SINDAT(6)*((CY1P5(6)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(6)+XPH)/(XPH**2+Y1P5SQ))
      ELSE

!         LOOP OVER DATA:
          XMH=X-HERM12(1)
          XPH=X+HERM12(1)
          RCPF12=SINDAT(1)*((CY1P5(1)-XMH)/(XMH**2+Y1P5SQ)              &
     &                     +(CY1P5(1)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(2)
          XPH=X+HERM12(2)
          RCPF12=RCPF12+SINDAT(2)*((CY1P5(2)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(2)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(3)
          XPH=X+HERM12(3)
          RCPF12=RCPF12+SINDAT(3)*((CY1P5(3)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(3)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(4)
          XPH=X+HERM12(4)
          RCPF12=RCPF12+SINDAT(4)*((CY1P5(4)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(4)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(5)
          XPH=X+HERM12(5)
          RCPF12=RCPF12+SINDAT(5)*((CY1P5(5)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(5)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(6)
          XPH=X+HERM12(6)
          RCPF12=RCPF12+SINDAT(6)*((CY1P5(6)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(6)+XPH)/(XPH**2+Y1P5SQ))
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION                                         &
     &  BMWID2(EDGE,HLFWDV,SUGCPI,GAMC,GAMC2,GAMD2,BLINES,LBMWID)

!     BMWID2 RETURNS A NORMALIZED FINITE BIN VOIGT EQUIVALENT WIDTH.

!                          EDGE
!                           /    /                       2   \
!          BMWID2(EDGE) =   |    |  1 - EXP[ -T(R BLINES)  ] |  d R
!                           /    \                           /
!                           0

!     WHERE
!                                                               2
!                                     1/2  INF        - (LN 2) Z
!                  2          / LN 2 \     /         E            DZ
!              T(N)  = SUGCPI | ---- |     |     -----------------------
!                             \  PI  /     /         2                 2
!                                         -INF   GAMC  + ( N - GAMD Z )

!                        2
!              AND   GAMD  = GAMD2

!     THE INPUT HLFWDV IS HALF THE TOTAL VOIGT EQUIVALENT WIDTH, I.E.,
!     HLFWDV=BMWID2(INF).  UNLESS T(EDGE BLINES) IS LARGE (>16),
!     BMWID2(EDGE) IS CALCULATED AS THE DIFFERENCE BETWEEN HLFWDV AND
!     THE TAIL CONTRIBUTION (BMTAIL).  THE BMTAIL CALCULATION ASSUMES
!     THAT THE LORENTZ AND DOPPLER HALF-WIDTHS ARE SMALL COMPARED TO THE
!     LINE TAIL OFFSET, I.E., GAMC << EDGE*BLINES & GAMD << EDGE*BLINES.

!                          /                                       \
!              2   SUGCPI  |      GAMD2 /     15 GAMD2 - 4 GAMC \  |
!          T(N)  = ------  |  1 + ----- | 3 + ----------------- |  |
!                  DENOM   |      DENOM \           DENOM       /  |
!                          \                                       /

!                                2
!          DENOM = GAMC2 + ELINES

!                     /                     2   \
!          BMTAIL = - | 1 - EXP [ -T(ELINES)  ] |
!                     \                         /

!                          T(ELINES)        /
!            SQRT(SUGCPI)      /         2  |     / 3 GAMD2 - GAMC2 \  2
!          + ------------      |  EXP (-t ) | 2 + | --------------- | t
!               ELINES         /            |     \     SUGCPI      /
!                              0            \

!                             2                         2      \
!                   / 15 GAMD2  - 10 GAMD2 GAMC2 - GAMC2  \  4 |
!                 + | ----------------------------------- | t  |  d t
!                   \                      2              /    |
!                                  4 SUGCPI                    /

!     THE ERROR FUNCTION INTEGRALS ARE EVALUATED DIFFERENTLY
!     DEPENDING ON THE MAGNITUDE OF T(ELINES)**2 [VARIABLE "OPTDEP"].

!     INPUT ARGUMENTS:
!       EDGE     LINE CENTER TO SPECTRAL BIN EDGE
!                DISTANCE OVER SPECTRAL BIN WIDTH.
!       HLFWDV   HALF TOTAL VOIGT EQUIVALENT WIDTH OVER BIN WIDTH.
!       SUGCPI   PRODUCT OF [1] ABSORPTION COEFFICIENT (S/d, LINE
!                STRENGTH OVER LINE SPACING), [2] COLUMN AMOUNT U
!                AND [3] GAMC, I.E., THE RATIO OF LORENTZ HALF-WIDTH
!                TO LINE SPACING, ALL DIVIDED BY PI.
!       GAMC     LORENTZ (COLLISION) HALF-WIDTH TO LINE SPACING RATIO.
!       GAMC2    SQUARE OF GAMC.
!       GAMD2    SQUARE OF THE RATIO OF THE DOPPLER HALF-WIDTH TO
!                LINE SPACING PARAMETER, ALL DIVIDED BY (2 LN2).
!       BLINES   NUMBER OF LINES IN FULL SPECTRAL BIN.
      DOUBLE PRECISION EDGE,HLFWDV,SUGCPI,GAMC,GAMC2,GAMD2,BLINES
      LOGICAL LBMWID

!     PARAMETERS:
!       PCON     CONSTANT USED IN RATIONAL APPROXIMATION TO THE ERROR
!                FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A1       LINEAR COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A2       QUADRATIC COEFFICIENT USED IN RATIONAL APPROXIMATION
!                TO THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A3       CUBIC COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A4       QUARTIC COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A5       5TH-ORDER COEFFICIENT USED IN RATIONAL APPROXIMATION
!                TO THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       RTPI     SQUARE ROOT OF PI.
!       RTPI4    SQUARE ROOT OF PI DIVIDED BY 4.
!       RTPIC    SQUARE ROOT OF PI TIMES 3/32.
      DOUBLE PRECISION PCON,A1,A2,A3,A4,A5,RTPI,RTPI4,RTPIC
      PARAMETER(PCON=.3275911D0,A1=.451673692D0,A2=-.504257335D0,       &
     &  A3=2.519390259D0,A4=-2.575644906D0,A5=1.881292140D0,            &
     &  RTPI=A1+A2+A3+A4+A5,RTPI4=RTPI/4,RTPIC=3*RTPI/32)

!     LOCAL VARIABLES:
!       COEF2    COEFFICIENT OF VOIGT LINE TAIL SERIES QUADRATIC TERM.
!       COEF4    COEFFICIENT OF VOIGT LINE TAIL SERIES QUARTIC TERM.
!       ELINES   NUMBER OF LINES BETWEEN LINE CENTER AND BIN EDGE.
!       ELINSQ   ELINES SQUARED.
!       DENLOR   LORENTZ LINE SHAPE FUNCTION DENOMINATOR EVALUATED AT
!                THE BIN EDGE AND DIVIDED BY THE LINE SPACING SQUARED.
!       STOR0    CONSTANT EQUAL TO 5 TIMES COEF2 PLUS GAMC2.
!       OPTDEP   OPTICAL DEPTH EVALUATED AT THE SPECTRAL BIN EDGE.
!       TAIL     5 LORENTZ LINE TAIL INTEGRALS.
!       STOR1    CONSTANT EQUAL TO 9 TIMES GAMD2 PLUS STOR0.
!       BMTAIL   TAIL CONTRIBUTION TO NORMALIZED VOIGT EQUIVALENT WIDTH.
!       C4RAT    COEF4 DIVIDED BY SUGCPI.
!       RTDEP    SQUARE ROOT OF THE OPTICAL DEPTH.
!       PHI      PRODUCT OF RTPI, THE EXPONENTIAL FUNCTION EVALUATED AT
!                OPTDEP, AND THE COMPLEMENTARY ERROR FUNCTION OF RTDEP.
!       RTRAT    RATIO OF THE SQUARE ROOT OF SUGCPI TO BLINES.
!       RAPP     RATIO USED IN RATIONAL APPROXIMATION TO ERROR FUNCTION.
      DOUBLE PRECISION COEF2,COEF4,ELINES,ELINSQ,DENLOR,STOR0,OPTDEP,   &
     &  STOR1,BMTAIL,C4RAT,RTDEP,PHI,RTRAT,RAPP,TAIL(5)

!     DEFINE REQUIRED LOCAL VARIABLES:
      COEF2=3*GAMD2-GAMC2
      COEF4=5*(COEF2-GAMC2)*GAMD2-GAMC2**2
      ELINES=BLINES*EDGE
      ELINSQ=ELINES**2
      DENLOR=GAMC2+ELINSQ
      STOR0=5*COEF2+GAMC2
      OPTDEP=SUGCPI*(1+GAMD2*(3+STOR0/DENLOR)/DENLOR)/DENLOR

!     BRANCH BASED ON MAGNITUDE OF THE OPTICAL DEPTH:
      IF(OPTDEP.LE..01D0)THEN

!         SMALL OPTICAL DEPTH:
          CALL LOREN5(ELINES,ELINSQ,GAMC,GAMC2,TAIL)
          STOR1=9*GAMD2+STOR0
          BMTAIL=SUGCPI*(TAIL(1)+GAMD2*(3*TAIL(2)+STOR0*TAIL(3))        &
     &      -(SUGCPI/2)*(TAIL(2)+GAMD2*(6*TAIL(3)+(STOR0+STOR1)*TAIL(4))&
     &      -SUGCPI*(TAIL(3)/3+GAMD2*(3*TAIL(4)+STOR1*TAIL(5)))))/BLINES
      ELSEIF(OPTDEP.GT.16.D0)THEN

!         SPECTRAL BIN IS BLACK:
          BMWID2=EDGE
          RETURN
      ELSE

!         LARGE OPTICAL DEPTH CONTRIBUTION:
          C4RAT=COEF4/SUGCPI
          RTRAT=SQRT(SUGCPI)/BLINES
          BMTAIL=RTRAT*(RTPI+(RTPI4*COEF2+RTPIC*C4RAT)/SUGCPI)-EDGE
          IF(OPTDEP.LT.11.D0)THEN

!             INTERMEDIATE OPTICAL DEPTH:
              RTDEP=SQRT(OPTDEP)
              IF(RTDEP.LT.2.33D0)THEN

!                 RATIONAL APPROXIMATION TO ERROR FUNCTION:
                  RAPP=1/(1+PCON*RTDEP)
                  PHI=RAPP*(A1+RAPP*(A2+RAPP*(A3+RAPP*(A4+RAPP*A5))))
              ELSE

!                 CONTINUED FRACTION REPRESENTATION FOR ERROR FUNCTION:
                  PHI=(2+OPTDEP*(4.5D0+OPTDEP))                         &
     &              /(RTDEP*(3.75D0+OPTDEP*(5+OPTDEP)))
              ENDIF
              BMTAIL=BMTAIL+EXP(-OPTDEP)                                &
     &          *(EDGE-RTRAT*(PHI+(2*COEF2*(2*RTDEP+PHI)+C4RAT          &
     &          *(RTDEP*(OPTDEP+1.5D0)+.75D0*PHI))/(8*SUGCPI)))
          ENDIF
      ENDIF
      BMWID2=HLFWDV-BMTAIL
      IF(20*BMWID2.LE.HLFWDV)LBMWID=.FALSE.
      RETURN
      END
