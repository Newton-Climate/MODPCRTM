      SUBROUTINE VSANSM(K,AHAZE,IHA1,ZNEW)
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
      INTEGER K

!     OUTPUT ARGUMENTS:
      INTEGER IHA1
      REAL AHAZE
      DOUBLE PRECISION ZNEW

!     COMMONS:

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

!     /CARD1B/
      INTEGER JUNITP,JUNITT,JUNIT,JLOW
      REAL WMOL
      COMMON/CARD1B/JUNITP,JUNITT,JUNIT(13),WMOL(12),JLOW

!     /ZVSALY/
      DOUBLE PRECISION ZVSA
      REAL RHVSA,AHVSA
      INTEGER IHVSA
      COMMON/ZVSALY/ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)

!     /NSINP/
      DOUBLE PRECISION ZNSINP
      REAL PNSINP,TNSINP,WNSINP
      COMMON/NSINP/ZNSINP(40),PNSINP(40),TNSINP(40),WNSINP(40,12)

!     /M_PTWO/
!       PPROF    PRESSURE PROFILE [MB].
!       TPROF    TEMPERATURE PROFILE [K].
!       WH2O     H2O VOLUME MIXING RATIO PROFILE [PPMV].
!       WO3      O3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL PPROF,TPROF,WH2O,WO3
      COMMON/M_PTWO/PPROF(LAYDIM),TPROF(LAYDIM),WH2O(LAYDIM),WO3(LAYDIM)

!     /M_UMIX/
!       WCO2     CO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WN2O     N2O VOLUME MIXING RATIO PROFILE [PPMV].
!       WCO      CO VOLUME MIXING RATIO PROFILE [PPMV].
!       WCH4     CH4 VOLUME MIXING RATIO PROFILE [PPMV].
!       WO2      O2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WNO      NO VOLUME MIXING RATIO PROFILE [PPMV].
!       WSO2     SO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WNO2     NO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WHNO3    HNO3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL WCO2,WN2O,WCO,WCH4,WO2,WNO,WSO2,WNO2,WNH3,WHNO3
      COMMON/M_UMIX/WCO2(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),              &
     &  WCH4(LAYDIM),WO2(LAYDIM),WNO(LAYDIM),WSO2(LAYDIM),              &
     &  WNO2(LAYDIM),WNH3(LAYDIM),WHNO3(LAYDIM)

!     LOCAL VARIABLES:
      DOUBLE PRECISION DIF
      INTEGER JML,J,KN,JL,JP,JLS,KM
      REAL DLIN,FAC,WWMOL(12)

!     MODEL 7 CODING
!     OLD LAYERS  AEROSOL RETURNED
!     NEW LAYERS P,T,DP,AEROSOL
      JML=ML
      J=1
      KN=K
      IF(KN.GT.10)GO TO 140
  110 CONTINUE
      JL=J-1
      IF(JL.LT.1)JL=1
      JP=JL+1
      JLS = JL
      IF(ZVSA(KN).EQ.ZNSINP(JL))GO TO 140
      JLS = JP
      IF(ZVSA(KN).EQ.ZNSINP(JP))GO TO 140
      IF(ZVSA(KN).GT.ZNSINP(JL).AND.ZVSA(KN).LT.ZNSINP(JP))GO TO 115
      IF(J.GE.JML) GO TO 115
      J = J + 1
      GO TO 110
115   ZNEW=ZVSA(KN)
      DIF=ZNSINP(JP)-ZNSINP(JL)
      DLIN=SNGL((ZVSA(KN)-ZNSINP(JL))/DIF)
      PPROF(K)  = (PNSINP(JP)-PNSINP(JL))*DLIN+PNSINP(JL)
      TPROF(K)   =(TNSINP(JP)-TNSINP(JL))*DLIN+TNSINP(JL)
      DO KM = 1,12
          WWMOL(KM)=(WNSINP(JP,KM)-WNSINP(JL,KM))*DLIN+WNSINP(JL,KM)
      ENDDO
      IHA1  =IHVSA(KN)
      AHAZE  =AHVSA(KN)
      FAC=SNGL((ZVSA(KN)-ZNSINP(JL))/DIF)
      IF(PNSINP(JP).GT.0..AND.PNSINP(JL).GT.0.) THEN
           PPROF(K)  =PNSINP(JL)*(PNSINP(JP)/PNSINP(JL))**FAC
      ENDIF
      IF(TNSINP(JP).GT.0..AND.TNSINP(JL).GT.0.) THEN
           TPROF(K)   =TNSINP(JL)*(TNSINP(JP)/TNSINP(JL))**FAC
      ENDIF
      DO KM = 1,12
      IF(WNSINP(JP,KM).GT.0. .AND. WNSINP(JL,KM).GT.0.)                 &
     &  WWMOL(KM)=(WNSINP(JL,KM)*WNSINP(JP,KM)/WNSINP(JL,KM))**FAC
      ENDDO
       WH2O(K)    = WWMOL(1)
       WCO2(K)  = WWMOL(2)
       WO3(K)    = WWMOL(3)
       WN2O(K)  = WWMOL(4)
       WCO(K)   = WWMOL(5)
       WCH4(K)  = WWMOL(6)
       WO2(K)   = WWMOL(7)
       WNO(K)   = WWMOL(8)
       WSO2(K)  = WWMOL(9)
       WNO2(K)  = WWMOL(10)
       WNH3(K)  = WWMOL(11)
       WHNO3(K) = WWMOL(12)
      RETURN
140   CONTINUE
      IF(K.GT.10) THEN
         J = K - 10 + JLOW
         IHA1  =0
         AHAZE  =0.
      ELSE
         J = JLS
      ENDIF
      ZNEW=ZNSINP(J)
      PPROF(K)  =PNSINP(J)
      TPROF(K) = TNSINP(J)
      DO KM = 1,12
      WWMOL(KM)= WNSINP(J,KM)
      ENDDO
       WH2O(K)    = WWMOL(1)
       WCO2(K)  = WWMOL(2)
       WO3(K)    = WWMOL(3)
       WN2O(K)  = WWMOL(4)
       WCO(K)   = WWMOL(5)
       WCH4(K)  = WWMOL(6)
       WO2(K)   = WWMOL(7)
       WNO(K)   = WWMOL(8)
       WSO2(K)  = WWMOL(9)
       WNO2(K)  = WWMOL(10)
       WNH3(K)  = WWMOL(11)
       WHNO3(K) = WWMOL(12)
      IF(KN.LE.9) IHA1  =IHVSA(KN)
      IF(KN.LE.9)AHAZE  =AHVSA(KN)
      RETURN
      END
