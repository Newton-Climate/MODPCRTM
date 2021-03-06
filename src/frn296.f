      SUBROUTINE FRN296(V1C,FH2O)

!     LOADS FOREIGN CONTINUUM 296K

!     INPUT ARGUMENTS:
!       V1C      SPECTRAL FREQUENCY [CM-1].
!       FH2O     FOREIGN BROADENED H2O CONTINUUM DATA [1E-20*CM**3/MOL].
      REAL V1C,FH2O

!     COMMONS:

!     /FH2ODT/
!       V1FH2O   INITIAL SPECTRAL FREQUENCY [CM-1].
!       DVFH2O   SPECTRAL FREQUENCY STEP SIZE [CM-1].
!       N_FH2O   NUMBER OF SPECTRAL POINTS.
!       F296     FOREIGN BROADENED H2O CONTINUUM DATA [1E-20*CM**3/MOL].
      REAL V1FH2O,DVFH2O,F296
      INTEGER N_FH2O
      COMMON/FH2ODT/V1FH2O,DVFH2O,N_FH2O,F296(2003)
      SAVE /FH2ODT/

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL BFH2O
      CALL SINT(V1FH2O,V1C,DVFH2O,N_FH2O,F296,FH2O)
      RETURN
      END
