      LOGICAL FUNCTION TSTBAD(VARNAM,RELERR)

!     Write name (-VARNAM-) of variable failing self-test and its
!     percent error from the correct value;  return  'FALSE'.

      CHARACTER VARNAM*(*)
      REAL RELERR

      TSTBAD=.FALSE.
      WRITE(*, '(/,3A,1P,E11.2,A)')                                     &
     &       ' OUTPUT VARIABLE  ', VARNAM,'  DIFFERED BY', 100.*RELERR, &
     &       '  PER CENT FROM CORRECT VALUE.  SELF-TEST FAILED.'

      RETURN
      END
