      BLOCK DATA ATMCON

!     THIS SUBROUTINE INITIALIZES COMMONS /CONSTN/.

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:

!     /CONSTN/
!       AMWT     MOLECULAR WEIGHTS [GM/MOLE].
!                  1  H2O   2  CO2   3  O3    4  N2O   5  CO    6  CH4
!                  7  O2    8  NO    9  SO2   10 NO2   11 NH3   12 HNO3
      REAL AMWT
      COMMON/CONSTN/AMWT(12)
      SAVE /CONSTN/
      DATA AMWT/H2OMWT, 44.010, 47.998, 44.010, 28.011, 16.043,         &
     &  31.999, 30.010, 64.060, 46.010, 17.030, 63.010/
      END
