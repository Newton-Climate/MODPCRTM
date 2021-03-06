      DOUBLE PRECISION FUNCTION DPANDX(ECD,HEIGHT,RSCLHT,GAMMA0)

!     DPANDX COMPUTES THE PRODUCT OF THE EARTH CENTER DISTANCE AND
!     THE INDEX OF REFRACTION.
      IMPLICIT NONE

!     ARGUMENTS:
!       ECD      EARTH CENTER DISTANCE [KM].
!       HEIGHT   HEIGHT ABOVE THE SPHERICAL EARTH [KM].
!       RSCLHT   RECIPROCAL OF THE REFRACTIVITY SCALE HEIGHT [KM-1].
!       GAMMA0   REFRACTIVITY, EQUAL TO INDEX OF REFRACTION MINUS ONE.
      DOUBLE PRECISION ECD,HEIGHT,RSCLHT,GAMMA0

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

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

!     LOCAL VARIABLES:
!       ARG      ARGUMENT OF THE EXPONENTIAL.
!       RFRCTV   REFRACTIVITY AT HEIGHT.
!       REFIND   REFRACTIVE INDEX AT HEIGHT.
!       NRGRAD   GRADIENT OF THE REFRACTIVE INDEX AND ECD PRODUCT.
      DOUBLE PRECISION ARG,RFRCTV,REFIND,NRGRAD
      IF(HEIGHT.LT.GNDALT)THEN
          ARG=RSCLHT*GNDALT
      ELSE
          ARG=RSCLHT*HEIGHT
      ENDIF
      IF(ARG.GT.BIGEXP)THEN
          DPANDX=ECD
      ELSE
          RFRCTV=GAMMA0*EXP(-ARG)
          REFIND=1+RFRCTV
          DPANDX=ECD*REFIND
          NRGRAD=REFIND-RSCLHT*ECD*RFRCTV
          IF(LGEOM .AND. NRGRAD.LT.0.)THEN
              WRITE(IPR,'(/A,/(22X,2(A,F10.5)))')                       &
     &          ' Warning from DPANDX:  The product of the Earth'//     &
     &          ' center distance R+H and the index of refraction n(H)',&
     &          ' decreases with altitude H at H =',HEIGHT,             &
     &          ' km, with d[(R+H) * n(H)] / dH =',NRGRAD,              &
     &          ' The spherical refractive geometry routines may'//     &
     &          ' fail for path elevation angles near zero degrees.'
              LGEOM=.FALSE.
          ENDIF
      ENDIF
      RETURN
      END
