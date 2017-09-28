      REAL FUNCTION BT_WN(CHNRAD,MOM1ST,M2DIFF)

!     COMPUTES SPECTRAL CHANNEL BRIGHTNESS TEMPERATURE [K] FROM
!     RADIANCE SPECTRALLY CONVOLVED OVER A FREQUENCY (wavenumbers)
!     DEPENDENT CHANNEL RESPONSE FUNCTION.
      IMPLICIT NONE

!     PARAMETERS:
!       C1BYPI   FIRST RADIATION CONSTANT DIVIDED PI [W CM^2 / SR].
!       C2       SECOND RADIATION CONSTANT [CM K].
      REAL C1BYPI,C2
      PARAMETER(C1BYPI=1.190956E-12,C2=1.438786)

!     INPUT ARGUMENTS:
!       CHNRAD   CHANNEL SPECTRAL RADIANCE [ W CM-2 SR-1 / CM-1 ].
!       MOM1ST   CHANNEL FIRST SPECTRAL FREQUENCY (v) MOMENT [CM-1]:

!                     /                   where RSR is the normalized
!                V  = |  v RSR(v) dv ,    spectral frequency dependent
!                 1   /                   Relative Spectral Response
!                                         function.

!       M2DIFF   SPECTRAL FREQUENCY 2ND MOMENT DIFFERENCE INTEGRAL.

!                       /   2
!                = V   /  V     -   1
!                   2 /    1

      REAL CHNRAD,MOM1ST,M2DIFF

!     LOCAL VARIABLES:
!       RHO      DIMENSIONLESS CHANNEL SPECTRAL RADIANCE.
!       LN       DENOMINATOR IN ZERO ORDER APPROXIMATION.
      REAL RHO,LN

!     SQUARE SLIT APPROXIMATION WITH ( BNDWID / VCEN ) ^ 2 SMALL:
      RHO=CHNRAD/(C1BYPI*MOM1ST**3)
      LN=LOG(1+1/RHO)
      BT_WN=C2*MOM1ST/(LN+M2DIFF*(3/(1+RHO)-LN*(3-(.5+RHO)*LN)))
      RETURN
      END

      REAL FUNCTION BT_UM(CHNRAD,MOM1ST,M2DIFF)

!     COMPUTES SPECTRAL CHANNEL BRIGHTNESS TEMPERATURE [K] FROM
!     RADIANCE SPECTRALLY CONVOLVED OVER A WAVELENGTH (microns)
!     DEPENDENT CHANNEL RESPONSE FUNCTION.
      IMPLICIT NONE

!     PARAMETERS:
!       C1BYPI   FIRST RADIATION CONSTANT DIVIDED PI [W CM-2 uM^4 / SR].
!       C2       SECOND RADIATION CONSTANT [MICRONS K].
      REAL C1BYPI,C2
      PARAMETER(C1BYPI=11909.56,C2=14387.86)

!     INPUT ARGUMENTS:
!       CHNRAD   CHANNEL SPECTRAL RADIANCE [ W CM-2 SR-1 / MICRON ].
!       MOM1ST   CHANNEL FIRST SPECTRAL WAVELENGTH (w) MOMENT [MICRONS]:

!                     /                   where RSR is the normalized
!                W  = |  w RSR(w) dw ,    spectral wavelength dependent
!                 1   /                   Relative Spectral Response
!                                         function.

!       M2DIFF   SPECTRAL WAVELENGTH 2ND MOMENT DIFFERENCE INTEGRAL.

!                       /   2
!                = W   /  W     -   1
!                   2 /    1

      REAL CHNRAD,MOM1ST,M2DIFF

!     LOCAL VARIABLES:
!       SIG      DIMENSIONLESS CHANNEL SPECTRAL RADIANCE.
!       LN       DENOMINATOR IN ZERO ORDER APPROXIMATION.
      REAL SIG,LN

!     SQUARE SLIT APPROXIMATION WITH ( BNDWID / VCEN ) ^ 2 SMALL:
      SIG=CHNRAD*MOM1ST**5/C1BYPI
      LN=LOG(1+1/SIG)
      BT_UM=C2/(MOM1ST*(LN+M2DIFF*(15/(1+SIG)-LN*(6-(.5+SIG)*LN))))
      RETURN
      END

      REAL FUNCTION BT_NM(CHNRAD,MOM1ST,M2DIFF)

!     COMPUTES SPECTRAL CHANNEL BRIGHTNESS TEMPERATURE [K] FROM
!     RADIANCE SPECTRALLY CONVOLVED OVER A WAVELENGTH (nanometers)
!     DEPENDENT CHANNEL RESPONSE FUNCTION.
      IMPLICIT NONE

!     PARAMETERS:
!       C1BYPI   FIRST RADIATION CONSTANT DIVIDED PI [W CM-2 NM^4 / SR].
!       C2       SECOND RADIATION CONSTANT [NM K].
      REAL C1BYPI,C2
      PARAMETER(C1BYPI=1.190956E+16,C2=1.438786E+07)

!     INPUT ARGUMENTS:
!       CHNRAD   CHANNEL SPECTRAL RADIANCE [ W CM-2 SR-1 / NM ].
!       MOM1ST   CHANNEL FIRST SPECTRAL WAVELENGTH (w) MOMENT [NM]:

!                     /                   where RSR is the normalized
!                W  = |  w RSR(w) dw ,    spectral wavelength dependent
!                 1   /                   Relative Spectral Response
!                                         function.

!       M2DIFF   SPECTRAL WAVELENGTH 2ND MOMENT DIFFERENCE INTEGRAL.

!                       /   2
!                = W   /  W     -   1
!                   2 /    1

      REAL CHNRAD,MOM1ST,M2DIFF

!     LOCAL VARIABLES:
!       SIG      DIMENSIONLESS CHANNEL SPECTRAL RADIANCE.
!       LN       DENOMINATOR IN ZERO ORDER APPROXIMATION.
      REAL SIG,LN

!     SQUARE SLIT APPROXIMATION WITH ( BNDWID / VCEN ) ^ 2 SMALL:
      SIG=CHNRAD*MOM1ST**5/C1BYPI
      LN=LOG(1+1/SIG)
      BT_NM=C2/(MOM1ST*(LN+M2DIFF*(15/(1+SIG)-LN*(6-(.5+SIG)*LN))))
      RETURN
      END
