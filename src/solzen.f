      SUBROUTINE SOLZEN(L,IKMAX1,REE,DEG,ANGLE)

!     ROUTINE SOLZEN STORES THE COSINE OF THE SOLAR/LUNAR ZENITH
!     ANGLE DATA FOR THE MULTIPLE SCATTERING CALCULATIONS:
!       CSZEN0  COSINE OF SOLAR/LUNAR ZENITH AT THE GROUND.
!       CSZEN   LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH FOR VERTICAL
!               PATH (USED IN CALCULATION OF BACKSCATTER FRACTION).
!       CSZENX  AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING
!               (AWAY FROM EARTH) THE CURRENT LAYER.
      IMPLICIT NONE

!     DEFINE INPUTS
!       L       CURRENT LAYER BOUNDARY.
!       IKMAX1  NUMBER OF LAYER BOUNDARIES (NUMBER OF LAYERS PLUS ONE).
!       REE     EARTH RADIUS [KM]
!       DEG     NUMBER OF DEGREES IN ONE RADIAN
!       ANGLE   SOLAR ZENITH ANGLE AT BOUNDARY L [DEG]
      INTEGER L,IKMAX1
      DOUBLE PRECISION REE
      REAL DEG,ANGLE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:

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

!     /MSRD/
!       CSSCAT   COSINE OF THE SCATTERING ANGLE.
!                (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       SLEGEN   Nth LEGENDRE POLYNOMIAL EVALUATED AT THE COSINE OF THE
!                SCATTERING ANGLE TIMES (2N+1)/4pi (N=0 TO NSTR-1).
!       CSZEN0   LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.
!       CSZEN    LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.
!       CSZENX   AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING
!                (AWAY FROM EARTH) THE CURRENT LAYER.
!       BBGRND   THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].
!       BBNDRY   LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1].
!       TCONT    LAYER CONTINUUM OPTICAL DEPTH.
!       TAUT     LAYER TOTAL VERTICAL EXTINCTION OPTICAL DEPTH.
!       GTSCAT   SUM OVER SCATTERING SOURCES OF SCATTERING OPTICAL DEPTH
!                AND PHASE FUNCTION LEGENDRE COEFFICIENT PRODUCTS.
!       COSBAR   LAYER EFFECTIVE SCATTERING ASYMMETRY FACTOR.
!       DEPRAT   FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.
!       S0DEP    OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.
!       S0TRN    TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]
!       UPF      LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].
!       DNF      LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].
!       UPFS     LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].
!       DNFS     LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].
!       CO_LIN   TRUE IF LOS AND SOLAR PATHS ARE NEARLY IDENTICAL.
      REAL CSSCAT,SLEGEN,CSZEN0,CSZEN,CSZENX,TCONT,TAUT,GTSCAT,COSBAR,  &
     &  BBGRND,BBNDRY,S0DEP,S0TRN,DEPRAT,UPF,DNF,UPFS,DNFS
      LOGICAL CO_LIN
      COMMON/MSRD/CSSCAT(MLOS),SLEGEN(0:MAZ,MLOS),CSZEN0(LAYDIM),       &
     &  CSZEN(LAYDIM),CSZENX(LAYDIM),TCONT(LAYDIM),TAUT(MXKSUB,LAYDIM), &
     &  GTSCAT(0:MXCMU,1:LAYDIM),COSBAR(LAYDIM),BBGRND,BBNDRY(LAYDIM),  &
     &  S0DEP(MXKSUB,LAYTWO),S0TRN(MXKSUB,LAYTWO),DEPRAT(MXKSUB,LAYDIM),&
     &  UPF(MXKSUB,LAYDIM),DNF(MXKSUB,LAYDIM),UPFS(MXKSUB,LAYDIM),      &
     &  DNFS(MXKSUB,LAYDIM),CO_LIN(MLOS)

!     DECLARE LOCAL VARIABLES
      INTEGER LP1
      REAL RIDIF,DSTDIF,DSTRAT,PROD

!     COSINE OF THE SOLAR ZENITH AT UPPER LAYER BOUNDARY
      CSZEN0(L)=COS(ANGLE/DEG)

!     COSINE OF THE SOLAR ZENITH AVERAGED FOR LAYER.
      IF(L.GT.1)CSZEN(L-1)=(CSZEN0(L-1)+CSZEN0(L))/2

!     COSINE OF THE SOLAR ZENITH EXITING THE LAYER ABOVE.
      IF(L.GE.IKMAX1)RETURN
      LP1=L+1
      RIDIF=(RFNDX(L)-RFNDX(LP1))/(1+RFNDX(LP1))
      DSTDIF=SNGL((ZM(LP1)-ZM(L))/(REE+ZM(LP1)))
      DSTRAT=1-DSTDIF
      PROD=(1+RIDIF)*DSTRAT
      CSZENX(L)=(ABS(CSZEN0(L))+                                        &
     &  SQRT((CSZEN0(L)*PROD)**2+(1+PROD)*(DSTDIF-RIDIF*DSTRAT)))/2
      RETURN
      END