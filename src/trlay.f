      SUBROUTINE TRLAY(TAU,TSCAT,G,CSZEN,S0DEP,DEPRAT,EX,TDFS,REFS)

!     CALCULATE PARAMETERS FOR SOLAR HYBRID MODIFIED DELTA EDDINGTON
!     2-STREAM APPROXIMATION.  THE EQUATIONS ARE DERIVED FROM
!     W. E. MEADOR AND W. R. WEAVER, J. ATMOS. SCI. 37, 630-642 (1980).

!     INPUTS
!       TAU     VERTICAL EXTINCTION OPTICAL DEPTH.
!       TSCAT   VERTICAL SCATTERING OPTICAL DEPTH.
!       G       HENYEY-GREENSTEIN ASYMMETRY FACTOR.
!       CSZEN   COSINE OF THE SOLAR ZENITH ANGLE.
!       S0DEP   EXTINCTION OPTICAL DEPTH FROM LAYER BOTTOM TO SUN/MOON.
!       DEPRAT  FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN
!               ACROSS THE CURRENT LAYER.

!     OUTPUTS
!       EX      SOLAR ATTENUATION ACROSS LAYER
!       TDFS    LAYER TRANSMITTANCE MINUS EX
!       REFS    LAYER REFLECTANCE
      REAL TAU,TSCAT,G,S0DEP,DEPRAT,CSZEN,EX,TDFS,REFS

!     PARAMETERS:
      REAL CUTOFF,PT3CUT
      PARAMETER(CUTOFF=0.01,PT3CUT=.3*CUTOFF)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /LSTART/
!       TOLABS   OUTPUT WARNING THRESHOLD FOR EXCESS ABSORPTION.
!       LLORMN   OUTPUT WARNING FLAG FOR SMALL LORENTZ HALF WIDTH.
!       LLORMX   OUTPUT WARNING FLAG FOR LARGE LORENTZ HALF WIDTH.
!       LDOPMN   OUTPUT WARNING FLAG FOR SMALL DOPPLER HALF WIDTH.
!       LDOPMX   OUTPUT WARNING FLAG FOR LARGE DOPPLER HALF WIDTH.
!       LLINMN   OUTPUT WARNING FLAG FOR SMALL EFFECTIVE LINE NUMBER.
!       LSSA     OUTPUT WARNING FLAG FOR SCATTERING ALBEDO NEAR 1.
!       LAEREX   OUTPUT WARNING FLAG FOR ROUTINE AEREXT.
!       LBMCRK   OUTPUT WARNING FLAG FOR ROUTINE BMCRKS.
!       LFLUXS   OUTPUT WARNING FLAG FOR ROUTINE FLUXES.
!       LSSRAD   OUTPUT WARNING FLAG FOR ROUTINE SSRAD.
!       LTRLAY   OUTPUT WARNING FLAG FOR ROUTINE TRLAY.
!       LNVERS   OUTPUT WARNING FLAG FOR ROUTINE DENFAC.
!       LGEOM    OUTPUT WARNING FLAG FOR INDEX OF REFRACTION GRADIENT.
!       LO3TRN   OUTPUT WARNING FLAG FOR OZONE CURTIS-GODSON PROBLEM.
      REAL TOLABS
      LOGICAL LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,                  &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      COMMON/LSTART/TOLABS,LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,     &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      SAVE /LSTART/

!     /TP8CTL/
!       LTAPE8   .TRUE. IF TAPE8 IS BEING WRITTEN.
      LOGICAL LTAPE8
      COMMON/TP8CTL/LTAPE8

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      REAL BETABS

!     LOCAL VARIABLES:
      REAL OMEGA,MU0,GAMMA1,G1MG2,GAMMA3,GAMMA4,ONEMG2,DENOM,           &
     &  A1,A2,K,KT,KM,X,ONEPKM,ONEMKM,G3M,G4M,XHALF,TERM1,              &
     &  TERM2,TERM3,COEF,EKT,TWOEKT,RATP,RATM,XMKT,G3K,G4K


      IF(TAU.LE.0.)THEN

!         CASE 1:  TAU = 0.
          EX=1.
          REFS=0.
          TDFS=0.
      ELSE
          OMEGA=TSCAT/TAU
          GAMMA3=BETABS(CSZEN,G)
          GAMMA4=1-GAMMA3
          X=S0DEP*DEPRAT
          IF(TAU.GT.X)THEN

!             ONLY WARN USER IF NOT PREVIOUSLY WARNED
!             AND TAU EXCEEDS X BY MORE THAN 5%.
              IF(LTRLAY .AND. TAU.GT.1.05*X)THEN
                  IF(BINOUT)THEN
                      CALL BNWT3F(X,TAU,CSZEN)
                  ELSEIF(LTAPE8)THEN
                      WRITE(IPR1,'(/A,1P,2(E12.4,A,/22X,A),/(22X,A,0P,  &
     &                  F7.4,A))')' Warning from TRLAY:  The estimated' &
     &                  //' layer optical depth towards the sun is',X,  &
     &                  ', but','the vertical layer optical depth is',  &
     &                  TAU,'.  The solar path','layer optical depth'// &
     &                  ' was reset to the vertical layer extinction',  &
     &                  'optical depth over the cosine (=',CSZEN,       &
     &                  ') of the solar zenith.',                       &
     &                  '***  THIS WARNING WILL NOT BE REPEATED  ***'
                  ENDIF
                  LTRLAY=.FALSE.
              ENDIF
              MU0=ABS(CSZEN)+1.E-10
              X=TAU/MU0
          ELSE
              MU0=TAU/X
          ENDIF
          ONEMG2=(1-G)*(1+G)
          DENOM=ONEMG2*(1-MU0)+MU0
          GAMMA1=((1-OMEGA)+.75*ONEMG2*(1-G*OMEGA)                      &
     &      +GAMMA3*OMEGA*(1-ONEMG2))/DENOM
          G1MG2=(1+ONEMG2)*(1-OMEGA)/DENOM
          A1=GAMMA1-GAMMA3*G1MG2
          A2=GAMMA1-GAMMA4*G1MG2
          K=SQRT(G1MG2*(2*GAMMA1-G1MG2))
          KT=K*TAU
          KM=K*MU0
          ONEPKM=1+KM
          IF(X+KT.LT.CUTOFF)THEN

!             CASE 2:  (K TAU) AND (TAU/MU0) SMALL
!                      [(K + 1/MU0) TAU < CUTOFF]
              G3M=GAMMA3/MU0
              G4M=GAMMA4/MU0
              XHALF=X/2
              TERM1=ONEPKM+KM
              EX=1-X*(1-XHALF)
              COEF=TSCAT/(1-TAU*(K-GAMMA1))
              REFS=MAX(0.,COEF*(G3M-XHALF*(G3M*TERM1-A2)))
              TDFS=COEF*(G4M-XHALF*(G4M*TERM1-A1))
!NEXT
!NEXT         FOR GREATER ACCURACY, ADD THE NEXT HIGHER
!NEXT         ORDER TERM BY REPLACING THE EXPRESSIONS ABOVE
!NEXT         FOR EX, COEF, REFS AND TDFS WITH THOSE BELOW.
!NEXT         XTHIRD=X/3
!NEXT         TERM2=1+KM*(3+4*KM)
!NEXT         TERM3=TERM1+KM
!NEXT         EX=1-X*(1-XHALF*(1-XTHIRD))
!NEXT         COEF=TSCAT/(1-TAU*(K-GAMMA1)*(1-KT))
!NEXT         REFS=COEF*(G3M-XHALF*(G3M*TERM1-A2                        &
!NEXT&          -XTHIRD*(G3M*TERM2-A2*TERM3)))
!NEXT         TDFS=COEF*(G4M-XHALF*(G4M*TERM1-A1                        &
!NEXT&          -XTHIRD*(G4M*TERM2-A1*(1+TERM3))))
          ELSEIF(KT.LT.PT3CUT)THEN

!             CASE 3:  (K TAU) SMALL.
!                      [K TAU < .3 CUTOFF]
              COEF=OMEGA/(2-2*TAU*(K-GAMMA1)*(1-KT))
              EX=EXP(-X)
              TERM1=TAU*(1-KT*(3-KT)/6)
              EKT=1-K*TERM1
              TERM2=(1+EKT)*TERM1
              TERM3=1+EKT*EKT
              TWOEKT=2*EKT
              RATP=(1-EKT*EX)/ONEPKM
              ONEMKM=1-KM
              RATM=(EKT-EX)/ONEMKM
              DENOM=ONEPKM*ONEMKM
              REFS=MAX(0.,COEF*(GAMMA3*(RATP+EKT*RATM)                  &
     &          +A2*(MU0*(TWOEKT*EX-TERM3)+TERM2)/DENOM))
              TDFS=COEF*(GAMMA4*(RATM+EKT*RATP)                         &
     &          +A1*(MU0*(TWOEKT-EX*TERM3)-EX*TERM2)/DENOM)
          ELSE

!             CASE 4:  (K TAU) NOT SMALL.
!                      [K TAU > .3 CUTOFF]
              EKT=EXP(-KT)
              EX=EXP(-X)
              RATP=(1-EKT*EX)/ONEPKM
              XMKT=X-KT
              IF(ABS(XMKT).LT..004)THEN

!                 CASE 4A:  (K - 1/MU0) TAU NEAR ZERO.
                  RATM=X*EX*(1+XMKT*(3+XMKT)/6)
              ELSE

!                 CASE 4B:  (K - 1/MU0) TAU NOT NEAR ZERO.
                  RATM=(EKT-EX)/(1-KM)
              ENDIF
              COEF=OMEGA/(K+GAMMA1+(K-GAMMA1)*EKT**2)
              G3K=GAMMA3*K
              REFS=MAX(0.,COEF*((G3K+A2)*RATP+(G3K-A2)*RATM*EKT))
              G4K=GAMMA4*K
              TDFS=COEF*((G4K+A1)*RATM+(G4K-A1)*RATP*EKT)
          ENDIF
      ENDIF
      RETURN
      END
