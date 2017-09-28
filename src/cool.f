      SUBROUTINE COOL(VCEN)

!     THIS ROUTINE CALCULATES AND WRITES OUT SPECTRAL COOLING RATES.
!     COOLING RATES, THE CHANGE IN TEMPERATURE T WITH TIME t, EQUAL
!     THE GRAVITATIONAL CONSTANT g OVER THE SPECIFIC HEAT OF AIR AT
!     CONSTANT PRESSURE Cp TIMES THE CHANGE IN NET UPWARD FLUX F
!     WITH PRESSURE P:

!            dT     g   dF     g   1      dF
!            --  =  --  --  =  --  -  ----------
!            dt     Cp  dP     Cp  P  d(ln P/Po)

!     Po IS STANDARD PRESSURE (Po = 1 ATM).  THE DERIVATIVE OF
!     THE NET FLUX WITH LOG PRESSURE IS DEFINED AT EACH LEVEL BY
!     FITTING THE VALUES OF F AND ln P/Po FROM TWO LEVELS BELOW
!     TO TWO LEVELS ABOVE THE CURRENT LEVEL TO A POLYNOMIAL
!     (FOURTH DEGREE EXCEPT AT THE TOP AND BOTTOM OF THE
!     ATMOSPHERE), AND ANALYTICALLY COMPUTING THE DERIVATIVE.
      IMPLICIT NONE

!     ARGUMENTS
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
      REAL VCEN

!     PARAMETERS:
!       EPSILN   A SMALL POSITIVE NUMBER.
      INCLUDE 'PARAMS.h'
      REAL EPSILN
      PARAMETER(EPSILN=1.E-30)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /COOLRT/
!       ISPCCR  UNIT NUMBER FOR SPECTRAL COOLING RATE FILE.
!       P3COEF  PRESSURE DEPENDENT COEFFICIENTS FOR CALCULATING COOLING
!               RATES WITH 3 POINT SPLINE FIT [(K/DAY) / (W/CM2)].
!       P5COEF  PRESSURE DEPENDENT COEFFICIENTS FOR CALCULATING COOLING
!               RATES WITH 5 POINT SPLINE FIT [(K/DAY) / (W/CM2)].
!       CLRTSM  LEVEL IN-BAND COOLING RATES [K/DAY].
!       CLRTS0  LEVEL IN-BAND COOLING RATES EXCLUDING
!               DIRECT SOLAR CONTRIBUTION [K/DAY].
!       NTFSM   LEVEL IN-BAND NET (THERMAL PLUS SCATTERED
!               SOLAR PLUS DIRECT SOLAR) UPWARD FLUX [W/CM2].
!       UPFSM   LEVEL IN-BAND UPWARD THERMAL (PLUS SCATTERED
!               SOLAR IF DISORT & NO AZIMUTH DEPENDENCE) FLUX [W/CM2].
!       DNFSM   LEVEL IN-BAND DOWNWARD THERMAL (PLUS SCATTERED
!               SOLAR IF DISORT & NO AZIMUTH DEPENDENCE) FLUX [W/CM2].
!       UPFSSM  LEVEL IN-BAND UPWARD SCATTERED SOLAR FLUX
!               (USED WITH DISORT ONLY IF AZIMUTH DEPENDENT) [W/CM2].
!       DNFSSM  LEVEL IN-BAND DOWNWARD SCATTERED SOLAR FLUX
!               (USED WITH DISORT ONLY IF AZIMUTH DEPENDENT) [W/CM2].
      DOUBLE PRECISION CLRTSM,CLRTS0,NTFSM,UPFSM,DNFSM,UPFSSM,DNFSSM
      REAL P3COEF,P5COEF
      INTEGER ISPCCR
      COMMON/COOLRT/CLRTSM(LAYDIM),CLRTS0(LAYDIM),NTFSM(LAYDIM),        &
     &  UPFSM(LAYDIM),DNFSM(LAYDIM),UPFSSM(LAYDIM),DNFSSM(LAYDIM),      &
     &  P3COEF(3,LAYDIM),P5COEF(5,LAYDIM),ISPCCR

!     /NETFLX/
!       SUBINT   SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.
!       UPFLX    LEVEL UPWARD THERMAL SPECTRAL
!                FLUX (INCLUDES SCATTERED SOLAR IF DISORT
!                & NO AZIMUTH DEPENDENCE) [W CM-2 / CM-1].
!       DNFLX    LEVEL DOWNWARD THERMAL SPECTRAL
!                FLUX (INCLUDES SCATTERED SOLAR IF DISORT
!                & NO AZIMUTH DEPENDENCE) [W CM-2 / CM-1].
!       UPFLXS   LEVEL UPWARD SCATTERED SOLAR
!                SPECTRAL FLUX (USED WITH DISORT ONLY
!                IF AZIMUTH DEPENDENT) [W CM-2 / CM-1].
!       DNFLXS   LEVEL DOWNWARD SCATTERED SOLAR
!                SPECTRAL FLUX (USED WITH DISORT ONLY
!                IF AZIMUTH DEPENDENT) [W CM-2 / CM-1].
!       NTFLX    LEVEL NET (THERMAL PLUS SCATTERED SOLAR PLUS
!                DIRECT SOLAR) UPWARD SPECTRAL FLUX [W CM-2 / CM-1].
      REAL SUBINT,UPFLX,DNFLX,UPFLXS,DNFLXS,NTFLX
      COMMON/NETFLX/SUBINT(MXKSUB),UPFLX(LAYDIM),DNFLX(LAYDIM),         &
     &  UPFLXS(LAYDIM),DNFLXS(LAYDIM),NTFLX(LAYDIM)

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

!     BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     VARIABLES:
!       LM2,LM1,L,LP1,LP2   LEVEL INDICES.
!       CLRATE   SPECTRAL COOLING RATES [K DAY-1 / CM-1].
!       UPDFUS   LEVEL UPWARD DIFFUSE FLUX [W CM-2 / CM-1].
!       DNDFUS   LEVEL DOWNWARD DIFFUSE FLUX [W CM-2 / CM-1].
!       DNTOTL   LEVEL DOWNWARD TOTAL FLUX [W CM-2 / CM-1].
!       LUDFUS   NATURAL LOGARITHM OF LEVEL UPWARD DIFFUSE FLUX.
!       LDDFUS   NATURAL LOGARITHM OF LEVEL DOWNWARD DIFFUSE FLUX.
!       LDTOTL   NATURAL LOGARITHM OF LEVEL DOWNWARD TOTAL FLUX.
      INTEGER LM2,LM1,L,LP1,LP2
      REAL CLRATE(LAYDIM),UPDFUS(LAYDIM),DNDFUS(LAYDIM),DNTOTL(LAYDIM), &
     &  LUDFUS(LAYDIM),LDDFUS(LAYDIM),LDTOTL(LAYDIM)

!     INITIALIZE FIRST THREE LEVELS OF LOCAL ARRAYS:
      UPDFUS(1)=UPFLX(1)+UPFLXS(1)
      DNDFUS(1)=DNFLX(1)+DNFLXS(1)
      DNTOTL(1)=UPDFUS(1)-NTFLX(1)
      IF(UPDFUS(1).GT.EPSILN)THEN
          LUDFUS(1)=LOG(UPDFUS(1))
      ELSE
          LUDFUS(1)=-99.
      ENDIF
      IF(DNDFUS(1).GT.EPSILN)THEN
          LDDFUS(1)=LOG(DNDFUS(1))
      ELSE
          LDDFUS(1)=-99.
      ENDIF
      IF(DNTOTL(1).GT.EPSILN)THEN
          LDTOTL(1)=LOG(DNTOTL(1))
      ELSE
          LDTOTL(1)=-99.
      ENDIF
      UPDFUS(2)=UPFLX(2)+UPFLXS(2)
      DNDFUS(2)=DNFLX(2)+DNFLXS(2)
      DNTOTL(2)=UPDFUS(2)-NTFLX(2)
      IF(UPDFUS(2).GT.EPSILN)THEN
          LUDFUS(2)=LOG(UPDFUS(2))
      ELSE
          LUDFUS(2)=-99.
      ENDIF
      IF(DNDFUS(2).GT.EPSILN)THEN
          LDDFUS(2)=LOG(DNDFUS(2))
      ELSE
          LDDFUS(2)=-99.
      ENDIF
      IF(DNTOTL(2).GT.EPSILN)THEN
          LDTOTL(2)=LOG(DNTOTL(2))
      ELSE
          LDTOTL(2)=-99.
      ENDIF
      UPDFUS(3)=UPFLX(3)+UPFLXS(3)
      DNDFUS(3)=DNFLX(3)+DNFLXS(3)
      DNTOTL(3)=UPDFUS(3)-NTFLX(3)
      IF(UPDFUS(3).GT.EPSILN)THEN
          LUDFUS(3)=LOG(UPDFUS(3))
      ELSE
          LUDFUS(3)=-99.
      ENDIF
      IF(DNDFUS(3).GT.EPSILN)THEN
          LDDFUS(3)=LOG(DNDFUS(3))
      ELSE
          LDDFUS(3)=-99.
      ENDIF
      IF(DNTOTL(3).GT.EPSILN)THEN
          LDTOTL(3)=LOG(DNTOTL(3))
      ELSE
          LDTOTL(3)=-99.
      ENDIF

!     GROUND SPECTRAL COOLING RATE UPWARD COMPONENT.
      IF(UPDFUS(1).GT.EPSILN .AND. UPDFUS(3).GT.EPSILN)THEN

!         FOUR POINT SPLINE FIT CENTERED ON LEVEL 2:
          CLRATE(1)=UPDFUS(1)*    (P5COEF(3,1)*LUDFUS(1)                &
     &      +P5COEF(4,1)*LUDFUS(2)+P5COEF(5,1)*LUDFUS(3))
          CLRTS0(1)=CLRTS0(1)+DBLE(CLRATE(1))
      ELSE

!         NOT ENOUGH INFORMATION TO CALCULATE COOLING RATE:
          CLRATE(1)=0.
      ENDIF

!     GROUND SPECTRAL COOLING RATE DOWNWARD COMPONENT.
      IF(DNTOTL(1).GT.EPSILN .AND. DNTOTL(3).GT.EPSILN)                 &
     &  CLRATE(1)=CLRATE(1)-DNTOTL(1)*(P5COEF(3,1)*LDTOTL(1)            &
     &          +P5COEF(4,1)*LDTOTL(2)+P5COEF(5,1)*LDTOTL(3))
      CLRTSM(1)=CLRTSM(1)+DBLE(CLRATE(1))

!     GROUND SPECTRAL COOLING RATE DIFFUSE DOWNWARD COMPONENT.
      IF(DNDFUS(1).GT.EPSILN .AND. DNDFUS(3).GT.EPSILN)                 &
     &  CLRTS0(1)=CLRTS0(1)-DBLE(DNDFUS(1)*(P5COEF(3,1)*LDDFUS(1)       &
     &          +P5COEF(4,1)*LDDFUS(2)+P5COEF(5,1)*LDDFUS(3)))

!     GROUND SPECTRAL COOLING RATE:
      UPFSM(1)=UPFSM(1)+DBLE(UPFLX(1))
      DNFSM(1)=DNFSM(1)+DBLE(DNFLX(1))
      UPFSSM(1)=UPFSSM(1)+DBLE(UPFLXS(1))
      DNFSSM(1)=DNFSSM(1)+DBLE(DNFLXS(1))
      NTFSM(1)=NTFSM(1)+DBLE(NTFLX(1))

!     INITIALIZE FOURTH LEVEL OF LOCAL ARRAYS:
      UPDFUS(4)=UPFLX(4)+UPFLXS(4)
      DNDFUS(4)=DNFLX(4)+DNFLXS(4)
      DNTOTL(4)=UPDFUS(4)-NTFLX(4)
      IF(UPDFUS(4).GT.EPSILN)THEN
          LUDFUS(4)=LOG(UPDFUS(4))
      ELSE
          LUDFUS(4)=-99.
      ENDIF
      IF(DNDFUS(4).GT.EPSILN)THEN
          LDDFUS(4)=LOG(DNDFUS(4))
      ELSE
          LDDFUS(4)=-99.
      ENDIF
      IF(DNTOTL(4).GT.EPSILN)THEN
          LDTOTL(4)=LOG(DNTOTL(4))
      ELSE
          LDTOTL(4)=-99.
      ENDIF

!     SECOND LEVEL SPECTRAL COOLING RATE UPWARD COMPONENT.
      IF(UPDFUS(1).GT.EPSILN .AND. UPDFUS(4).GT.EPSILN)THEN

!         FOUR POINT SPLINE FIT CENTERED ON LEVEL 2:
          CLRATE(2)=UPDFUS(2)*                                          &
     &      (P5COEF(2,2)*LUDFUS(1)+P5COEF(3,2)*LUDFUS(2)                &
     &      +P5COEF(4,2)*LUDFUS(3)+P5COEF(5,2)*LUDFUS(4))
          CLRTS0(2)=CLRTS0(2)+DBLE(CLRATE(2))
      ELSEIF(UPDFUS(1).GT.EPSILN .AND. UPDFUS(3).GT.EPSILN)THEN

!         THREE POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(2)=UPDFUS(2)*    (P3COEF(1,2)*LUDFUS(1)                &
     &      +P3COEF(2,2)*LUDFUS(2)+P3COEF(3,2)*LUDFUS(3))
          CLRTS0(2)=CLRTS0(2)+DBLE(CLRATE(2))
      ELSE

!         NOT ENOUGH INFORMATION TO CALCULATE COOLING RATE:
          CLRATE(2)=0.
      ENDIF

!     SECOND LEVEL SPECTRAL COOLING RATE DOWNWARD COMPONENT.
      IF(DNTOTL(1).GT.EPSILN .AND. DNTOTL(4).GT.EPSILN)THEN

!         FOUR POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(2)=CLRATE(2)-DNTOTL(2)*                                &
     &      (P5COEF(2,2)*LDTOTL(1)+P5COEF(3,2)*LDTOTL(2)                &
     &      +P5COEF(4,2)*LDTOTL(3)+P5COEF(5,2)*LDTOTL(4))
      ELSEIF(DNTOTL(1).GT.EPSILN .AND. DNTOTL(3).GT.EPSILN)THEN

!         THREE POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(2)=CLRATE(2)-DNTOTL(2)*(P3COEF(1,2)*LDTOTL(1)          &
     &            +P3COEF(2,2)*LDTOTL(2)+P3COEF(3,2)*LDTOTL(3))
      ENDIF
      CLRTSM(2)=CLRTSM(2)+DBLE(CLRATE(2))

!     SECOND LEVEL SPECTRAL COOLING RATE DIFFUSE DOWNWARD COMPONENT.
      IF(DNDFUS(1).GT.EPSILN .AND. DNDFUS(4).GT.EPSILN)THEN

!         FOUR POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRTS0(2)=CLRTS0(2)-DBLE(DNDFUS(2)*                           &
     &      (P5COEF(2,2)*LDDFUS(1)+P5COEF(3,2)*LDDFUS(2)                &
     &      +P5COEF(4,2)*LDDFUS(3)+P5COEF(5,2)*LDDFUS(4)))
      ELSEIF(DNDFUS(1).GT.EPSILN .AND. DNDFUS(3).GT.EPSILN)THEN

!         THREE POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRTS0(2)=CLRTS0(2)-DBLE(DNDFUS(2)*(P3COEF(1,2)*LDDFUS(1)     &
     &       +P3COEF(2,2)*LDDFUS(2)+P3COEF(3,2)*LDDFUS(3)))
      ENDIF

!     SECOND LEVEL SPECTRAL COOLING RATE:
      UPFSM(2)=UPFSM(2)+DBLE(UPFLX(2))
      DNFSM(2)=DNFSM(2)+DBLE(DNFLX(2))
      UPFSSM(2)=UPFSSM(2)+DBLE(UPFLXS(2))
      DNFSSM(2)=DNFSSM(2)+DBLE(DNFLXS(2))
      NTFSM(2)=NTFSM(2)+DBLE(NTFLX(2))

!     LOOP OVER INTERMEDIATE LAYERS:
      LM2=1
      LM1=2
      L=3
      LP1=4
      DO LP2=5,ML

!         INITIALIZE THE LP2 LEVEL OF LOCAL ARRAYS:
          UPDFUS(LP2)=UPFLX(LP2)+UPFLXS(LP2)
          DNDFUS(LP2)=DNFLX(LP2)+DNFLXS(LP2)
          DNTOTL(LP2)=UPDFUS(LP2)-NTFLX(LP2)
          IF(UPDFUS(LP2).GT.EPSILN)THEN
              LUDFUS(LP2)=LOG(UPDFUS(LP2))
          ELSE
              LUDFUS(LP2)=-99.
          ENDIF
          IF(DNDFUS(LP2).GT.EPSILN)THEN
              LDDFUS(LP2)=LOG(DNDFUS(LP2))
          ELSE
              LDDFUS(LP2)=-99.
          ENDIF
          IF(DNTOTL(LP2).GT.EPSILN)THEN
              LDTOTL(LP2)=LOG(DNTOTL(LP2))
          ELSE
              LDTOTL(LP2)=-99.
          ENDIF

!         INCREMENT SPECTRAL SUMS:
          UPFSM(L)=UPFSM(L)+DBLE(UPFLX(L))
          DNFSM(L)=DNFSM(L)+DBLE(DNFLX(L))
          UPFSSM(L)=UPFSSM(L)+DBLE(UPFLXS(L))
          DNFSSM(L)=DNFSSM(L)+DBLE(DNFLXS(L))
          NTFSM(L)=NTFSM(L)+DBLE(NTFLX(L))

!         INTERMEDIATE LEVEL SPECTRAL COOLING RATE UPWARD COMPONENT.
          IF(UPDFUS(LM2).GT.EPSILN .AND. UPDFUS(LP2).GT.EPSILN)THEN

!             FIVE POINT SPLINE FIT CENTERED ON LEVEL L:
              CLRATE(L)=UPDFUS(L)*      (P5COEF(1,L)*LUDFUS(LM2)        &
     &          +P5COEF(2,L)*LUDFUS(LM1)+P5COEF(3,L)*LUDFUS(L  )        &
     &          +P5COEF(4,L)*LUDFUS(LP1)+P5COEF(5,L)*LUDFUS(LP2))
              CLRTS0(L)=CLRTS0(L)+DBLE(CLRATE(L))
          ELSEIF(UPDFUS(LM1).GT.EPSILN .AND. UPDFUS(LP1).GT.EPSILN)THEN

!             THREE POINT SPLINE FIT CENTERED ON LEVEL L:
              CLRATE(L)=UPDFUS(L)*    (P3COEF(1,L)*LUDFUS(LM1)          &
     &          +P3COEF(2,L)*LUDFUS(L)+P3COEF(3,L)*LUDFUS(LP1))
              CLRTS0(L)=CLRTS0(L)+DBLE(CLRATE(L))
          ELSE

!             NOT ENOUGH INFORMATION TO CALCULATE COOLING RATE:
              CLRATE(L)=0.
          ENDIF

!         INTERMEDIATE LEVEL SPECTRAL COOLING RATE DOWNWARD COMPONENT.
          IF(DNTOTL(LM2).GT.EPSILN .AND. DNTOTL(LP2).GT.EPSILN)THEN

!             FIVE POINT SPLINE FIT CENTERED ON LEVEL L:
              CLRATE(L)=CLRATE(L)-DNTOTL(L)*(P5COEF(1,L)*LDTOTL(LM2)    &
     &              +P5COEF(2,L)*LDTOTL(LM1)+P5COEF(3,L)*LDTOTL(L  )    &
     &              +P5COEF(4,L)*LDTOTL(LP1)+P5COEF(5,L)*LDTOTL(LP2))
          ELSEIF(DNTOTL(LM1).GT.EPSILN .AND. DNTOTL(LP1).GT.EPSILN)THEN

!             THREE POINT SPLINE FIT CENTERED ON LEVEL L:
              CLRATE(L)=CLRATE(L)-DNTOTL(L)*(P3COEF(1,L)*LDTOTL(LM1)    &
     &                +P3COEF(2,L)*LDTOTL(L)+P3COEF(3,L)*LDTOTL(LP1))
          ENDIF
          CLRTSM(L)=CLRTSM(L)+DBLE(CLRATE(L))

!         INTERMEDIATE LEVEL SPECTRAL COOLING RATE DIFFUSE DOWNWARD
!         COMPONENT.
          IF(DNDFUS(LM2).GT.EPSILN .AND. DNDFUS(LP2).GT.EPSILN)THEN

!             FIVE POINT SPLINE FIT CENTERED ON LEVEL L:
              CLRTS0(L)=CLRTS0(L)-DBLE(DNDFUS(L)*(P5COEF(3,L)*LDDFUS(L) &
     &          +P5COEF(1,L)*LDDFUS(LM2)+P5COEF(2,L)*LDDFUS(LM1)        &
     &          +P5COEF(4,L)*LDDFUS(LP1)+P5COEF(5,L)*LDDFUS(LP2)))
          ELSEIF(DNDFUS(LM1).GT.EPSILN .AND. DNDFUS(LP1).GT.EPSILN)THEN

!             THREE POINT SPLINE FIT CENTERED ON LEVEL L:
              CLRTS0(L)=CLRTS0(L)-DBLE(DNDFUS(L)*(P3COEF(2,L)*LDDFUS(L) &
     &          +P3COEF(1,L)*LDDFUS(LM1)+P3COEF(3,L)*LDDFUS(LP1)))
          ENDIF
          LM2=LM1
          LM1=L
          L=LP1
          LP1=LP2
      ENDDO

!     SECOND-TO-LAST LEVEL SPECTRAL COOLING RATE UPWARD COMPONENT.
      IF(UPDFUS(LM2).GT.EPSILN .AND. UPDFUS(ML).GT.EPSILN)THEN

!         FOUR POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(L)=UPDFUS(L)*                                          &
     &      (P5COEF(1,L)*LUDFUS(LM2)+P5COEF(2,L)*LUDFUS(LM1)            &
     &      +P5COEF(3,L)*LUDFUS(L  )+P5COEF(4,L)*LUDFUS(ML ))
          CLRTS0(L)=CLRTS0(L)+DBLE(CLRATE(L))
      ELSEIF(UPDFUS(LM1).GT.EPSILN .AND. UPDFUS(ML).GT.EPSILN)THEN

!         THREE POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(L)=UPDFUS(L)*    (P3COEF(1,L)*LUDFUS(LM1)              &
     &      +P3COEF(2,L)*LUDFUS(L)+P3COEF(3,L)*LUDFUS(ML ))
          CLRTS0(L)=CLRTS0(L)+DBLE(CLRATE(L))
      ELSE

!         NOT ENOUGH INFORMATION TO CALCULATE COOLING RATE:
          CLRATE(L)=0.
      ENDIF

!     SECOND-TO-LAST LEVEL SPECTRAL COOLING RATE DOWNWARD COMPONENT.
      IF(DNTOTL(LM2).GT.EPSILN .AND. DNTOTL(ML).GT.EPSILN)THEN

!         FOUR POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(L)=CLRATE(L)-DNTOTL(L)*                                &
     &      (P5COEF(1,L)*LDTOTL(LM2)+P5COEF(2,L)*LDTOTL(LM1)            &
     &      +P5COEF(3,L)*LDTOTL(L  )+P5COEF(4,L)*LDTOTL(ML ))
      ELSEIF(DNTOTL(LM1).GT.EPSILN .AND. DNTOTL(ML).GT.EPSILN)THEN

!         THREE POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(L)=CLRATE(L)-DNTOTL(L)*(P3COEF(1,L)*LDTOTL(LM1)        &
     &            +P3COEF(2,L)*LDTOTL(L)+P3COEF(3,L)*LDTOTL(ML ))
      ENDIF

!     SECOND-TO-LAST LEVEL SPECTRAL COOLING RATE DIFFUSE DOWNWARD
!     COMPONENT.
      IF(DNDFUS(LM2).GT.EPSILN .AND. DNDFUS(ML ).GT.EPSILN)THEN

!         FOUR POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRTS0(L)=CLRTS0(L)-DBLE(DNDFUS(L)*                           &
     &      (P5COEF(1,L)*LDDFUS(LM2)+P5COEF(2,L)*LDDFUS(LM1)            &
     &      +P5COEF(3,L)*LDDFUS(L  )+P5COEF(4,L)*LDDFUS(ML )))
      ELSEIF(DNDFUS(LM1).GT.EPSILN .AND. DNDFUS(ML ).GT.EPSILN)THEN

!         THREE POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRTS0(L)=CLRTS0(L)-DBLE(DNDFUS(L)*(P3COEF(1,L)*LDDFUS(LM1)   &
     &                 +P3COEF(2,L)*LDDFUS(L)+P3COEF(3,L)*LDDFUS(ML )))
      ENDIF

!     SECOND-TO-LAST LEVEL SPECTRAL COOLING RATE.
      CLRTSM(L)=CLRTSM(L)+DBLE(CLRATE(L))
      UPFSM(L)=UPFSM(L)+DBLE(UPFLX(L))
      DNFSM(L)=DNFSM(L)+DBLE(DNFLX(L))
      UPFSSM(L)=UPFSSM(L)+DBLE(UPFLXS(L))
      DNFSSM(L)=DNFSSM(L)+DBLE(DNFLXS(L))
      NTFSM(L)=NTFSM(L)+DBLE(NTFLX(L))

!     TOP-OF-ATMOSPHERE SPECTRAL COOLING RATE UPWARD COMPONENT.
      IF(UPDFUS(LM1).GT.EPSILN .AND. UPDFUS(ML ).GT.EPSILN)THEN

!         THREE POINT SPLINE FIT CENTERED ON LEVEL L:
          CLRATE(ML)=UPDFUS(ML)*   (P5COEF(1,ML)*LUDFUS(LM1)            &
     &      +P5COEF(2,ML)*LUDFUS(L)+P5COEF(3,ML)*LUDFUS(ML ))
          CLRTS0(ML)=CLRTS0(ML)+DBLE(CLRATE(ML))
      ELSE

!         NOT ENOUGH INFORMATION TO CALCULATE COOLING RATE:
          CLRATE(L)=0.
      ENDIF

!     TOP-OF-ATMOSPHERE SPECTRAL COOLING RATE DOWNWARD COMPONENT.
      IF(DNTOTL(LM1).GT.EPSILN .AND. DNTOTL(ML).GT.EPSILN)              &
     &  CLRATE(ML)=CLRATE(ML)-DNTOTL(ML)*(P5COEF(1,ML)*LDTOTL(LM1)      &
     &            +P5COEF(2,ML)*LDTOTL(L)+P5COEF(3,ML)*LDTOTL(ML ))

!     TOP-OF-ATMOSPHERE SPECTRAL COOLING RATE DIFFUSE DOWNWARD COMPONENT
      IF(DNDFUS(LM1).GT.EPSILN .AND. DNDFUS(ML).GT.EPSILN)              &
     &  CLRTS0(ML)=CLRTS0(ML)-DBLE(DNDFUS(ML)*(P5COEF(1,ML)*LDDFUS(LM1) &
     &    +P5COEF(2,ML)*LDDFUS(L)+P5COEF(3,ML)*LDDFUS(ML )))

!     TOP OF ATMOSPHERE SPECTRAL COOLING RATE.
      CLRTSM(ML)=CLRTSM(ML)+DBLE(CLRATE(ML))
      UPFSM(ML)=UPFSM(ML)+DBLE(UPFLX(ML))
      DNFSM(ML)=DNFSM(ML)+DBLE(DNFLX(ML))
      UPFSSM(ML)=UPFSSM(ML)+DBLE(UPFLXS(ML))
      DNFSSM(ML)=DNFSSM(ML)+DBLE(DNFLXS(ML))
      NTFSM(ML)=NTFSM(ML)+DBLE(NTFLX(ML))

!     WRITE OUT SPECTRAL COOLING RATE DATA.
      IF(NPR.LE.-2)WRITE(ISPCCR,'(F8.2,1P,11E11.3:,/(8X,1P,11E11.3))')  &
     &  VCEN,(-CLRATE(L),L=1,ML)
      RETURN
      END