      DOUBLE PRECISION FUNCTION WTRGT(DIFF,COEF,AVG,WRIGHT,NORM)

!     WTRGT RETURNS A RIGHT SIDE WAVELENGTH TRIANGULAR SLIT INTEGRATION:

!               C / v1
!            /                                      2
!     WT  =  |            [ (w + D) - w' ] dw'  /  D
!            /
!               C / v2

!                         C    (               C    v1 + v2 )      2
!         =  (v2 - v1)  -----  ( (w + D)  -  -----  ------- )  /  D
!                       v1 v2  (             v1 v2     2    )

!     WHERE

!          v1 = minimum frequency [CM-1]
!          v2 = maximum frequency [CM-1]
!          C  = conversion factor [WAVELENGTHS / CM]
!          w  = central wavelength [MICRONS OR NANOMETERS]
!          D  = full-width-at-half-maximum [MICRONS OR NANOMETERS]

!     ARGUMENTS (WAVELENGTH UNIT CAN BE MICRONS OR NANOMETERS):
!       DIFF     FREQUENCY DIFFERENCE, v2 - v1 [CM-1].
!       COEF     COEFFICIENT, C / (v1 v2) [WAVELENGTH CM].
!       AVG      AVERAGE FREQUENCY, (v1 + v2) / 2 [CM-1].
!       WRIGHT   MAXIMUM SLIT WAVELENGTH, w - D [WAVELENGTH].
!       NORM     SLIT NORMALIZATION, D*D [WAVELENGTH**2].
      DOUBLE PRECISION DIFF,COEF,AVG,WRIGHT,NORM

!     DEFINE WEIGHT:
      WTRGT=DIFF*COEF*(WRIGHT-COEF*AVG)/NORM
      RETURN
      END
