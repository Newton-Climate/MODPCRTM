      LOGICAL FUNCTION WRTBAD(VARNAM)

!                   Write names of erroneous variables and return 'TRUE'
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!      INPUT :   VARNAM=name of erroneous variable to be written
!                         (character of any length)

      CHARACTER      VARNAM*(*)
      INTEGER        NUMMSG
      SAVE           NUMMSG
      DATA           NUMMSG / 0 /

      WRTBAD=.TRUE.
      NUMMSG=NUMMSG + 1
      WRITE (*, '(3A)')  ' ****  INPUT VARIABLE  ', VARNAM,             &
     &                     '  IN ERROR  ****'
      IF (NUMMSG.EQ.MAXMSG)                                             &
     &   CALL  ERRMSG ('TOO MANY INPUT ERRORS.  ABORTING...$', .TRUE.)

      RETURN
      END
